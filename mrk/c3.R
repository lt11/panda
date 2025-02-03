## header ---------------------------------------------------------------------

options(scipen = 999)
options(stringsAsFactors = F)
rm(list = ls())
library(data.table)
library(this.path)
library(scriptName)
library(GenomicRanges)
library(tictoc)

## function(s) ----------------------------------------------------------------

SplitSubCol <- function(x, n, s) {
  y <- sapply(strsplit(x, split = s, fixed = T), "[[", n)
  return(y)
}
### It splits a character column of a data-table and collects one sub-column.
### 
### Arguments:
### (1) x: a character column of a data-table
### (2) n: the number of sub-column to print in output
### (3) s: the string to apply the split
###
### Returns:
### (1) y: a character vector

## settings -------------------------------------------------------------------

### fixed settings
dirBase <- dirname(this.dir())
### dev dirBase <- "/Users/Lorenzo/dev/panda"
dirAnnoBed <- file.path(dirBase, "anno", "bed")
dirOut <- file.path(dirBase, "png")
unlink(dirOut, recursive = T)
dir.create(dirOut, showWarnings = F)
pathGen <- file.path(dirOut, "generators.txt")
cat("class:feature", "\n", file = pathGen)
idRef <- "SGDref"
hdImpg <- c("Query_id", "Query_start", "Query_end",
            "Target_id", "Target_start", "Target_end",
            "Target_clsfeat", "N_score", "Query_alndir", "Target_alndir")

# vtClassSrt <- c("gene",
#                 "pseudogene",
#                 "intron",
#                 "five_prime_UTR_intron",
#                 "noncoding_exon",
#                 "plus_1_translational_frameshift",
#                 "blocked_reading_frame",
#                 "external_transcribed_spacer_region",
#                 "internal_transcribed_spacer_region",
#                 "ARS",
#                 "ARS_consensus_sequence",
#                 "TY1",
#                 "TY1/TY2_soloLTR",
#                 "TY1_truncated",
#                 "TY2",
#                 "TY3_soloLTR",
#                 "TY4_soloLTR",
#                 "TSU4_soloLTR",
#                 "TY4_truncated",
#                 "TY5",
#                 "TY5_soloLTR",
#                 "LTR_retrotransposon",
#                 "long_terminal_repeat",
#                 "W_region",
#                 "Z1_region",
#                 "Z2_region",
#                 "centromere",
#                 "centromere_DNA_Element_I",
#                 "centromere_DNA_Element_II",
#                 "centromere_DNA_Element_III",
#                 "matrix_attachment_site",
#                 "X_element",
#                 "X_element_partial",
#                 "X_element_combinatorial_repeat",
#                 "X_region",
#                 "Y_prime_element",
#                 "Y_region",
#                 "recombination_enhancer",
#                 "silent_mating_type_cassette_array",
#                 "mating_type_region",
#                 "tRNA",
#                 "pseudogenic_transcript",
#                 "ncRNA",
#                 "ncRNA_gene",
#                 "non_transcribed_region",
#                 "rRNA",
#                 "rRNA_gene",
#                 "snRNA",
#                 "snRNA_gene",
#                 "snoRNA",
#                 "snoRNA_gene",
#                 "intein_encoding_region")

### the path to the paf file
dirAlnPaf <- file.path(dirBase, "aln")
pathAlnPaf <- list.files(path = dirAlnPaf, pattern = "paf$",
                         recursive = T, full.names = T)
if (length(pathAlnPaf) > 1) {
  stop("Multiple paf files found.")
} else if (length(pathAlnPaf) == 0) {
  stop("No paf file found.")
}

## clmnt ----------------------------------------------------------------------

### script name
myName <- current_filename()
cat("[", myName, "] ",
    "Making the pangenome ",
    "\n", sep = "")

### read strain-haplotypes ids from file (the order is maintained in the output)
pathIds <- file.path(dirBase, "ids", "ids-ps.txt")
vtStrainHaplo <- as.character(fread(file = pathIds, header = F)[[1]])
vtStrainHaploHash <- gsub(pattern = "-", replacement = "#", x = vtStrainHaplo)
### initialise the output features table
nHaplos <- length(vtStrainHaploHash)
dtPanFeats <- data.table(matrix(character(length = 0), ncol = nHaplos))
setnames(dtPanFeats, vtStrainHaploHash)
### dev 
# vtStrainHaploHash <- c("SGDref#0", "AFI#0", "S288C#0",
#                        "XXX-h3", "stoca-zo", "DBVPG6765#0")

### dev indS <- vtStrainHaplo[1] 
indL <- 1
lsImpg <- list()
### for each strain-haplotypes (hyphen-separated)
for (indS in vtStrainHaplo) {
  
  ### pick the strain-haplotypes bed file (query)
  pathAnnoBed <- list.files(path = dirAnnoBed, pattern = indS, full.names = T)
  ### string for bash
  strBashImpg <- paste0("impg -I -p ", pathAlnPaf, " -b ", pathAnnoBed)
  ### run "$ impg -p file.paf -b file.bed"  which is 
  ### (50x faster than going line-by-line)
  ### with direct load of the output
  strOut <- system(strBashImpg, intern = T)
  ### make a single-file data-table
  lsImpgOne <- strsplit(strOut, split = "\t")
  dtImpgOne <- rbindlist(lapply(lsImpgOne,
                                function(x) as.data.table(as.list(x))))
  ### append the single-file data-table to a list
  lsImpg[[indL]] <- dtImpgOne
  indL <- indL + 1
}

### get the redundant data-table
dtImpgAll <- rbindlist(lsImpg)
colnames(dtImpgAll) <- hdImpg
### add supplementary columns
dtImpgAll[, c("Query_start", "Query_end") := lapply(.SD, as.numeric),
          .SDcols = c("Query_start", "Query_end")]
dtImpgAll[, Target_cls := SplitSubCol(Target_clsfeat, 1, ":")]
dtImpgAll[, c("Target_start", "Target_end") := lapply(.SD, as.numeric),
          .SDcols = c("Target_start", "Target_end")]
dtImpgAll[, Target_len := Target_end - Target_start]

### size filters
dtImpgAllSzFlt <- dtImpgAll[Target_len > 5 & Query_end - Query_start > 5]

### sorting ascending or descending is not equivalent since the 
### overlap is done in two rounds and 
### is strand-aware (if two features have overlapping 
### coordinates but one is in the + strand and the other is in 
### the - strand they do not overlap), thus
### we sort in increasing (ascending) order to avoid spurious
### overlaps that may occur when multiple classes are analysed
### (a spurious overlap will produce a big block that cannot be easily
### decomposed, except for the different classes);
### the YFL066C locus (where we have the Y' region TEL06L) is a good
### example of a spurious overlap
setorder(dtImpgAllSzFlt, Target_len)

### transfer the strand data from Target_clsfeat to Query_id and Target_id
strStrand <- SplitSubCol(dtImpgAllSzFlt[, Target_clsfeat], 2, "#")
dtImpgAllSzFlt[, Query_id := paste0(Query_id, "#", strStrand)]
dtImpgAllSzFlt[, Target_id := paste0(Target_id, "#", strStrand)]

nBlocks <- 1
vtUnq <- unique(dtImpgAllSzFlt[, Target_clsfeat])

### dev
# indTarClsFeat <- "Y_prime_element:TEL06L#-"
# indTarClsFeat <- "gene:YAR014C#-"
# indTarClsFeat <- vtUnq[1]
# indTarClsFeat <- "gene:YBL037W#+"
# indTarClsFeat <- "X_element:TEL01L#-"
# indTarClsFeat <- "gene:YBR140C#-"
# indTarClsFeat <- "ARS:ARS102#+"
# indTarClsFeat <- "gene:YCR012W#+"
### quello con molti haplo-id doppioni, indSb <- "TY1/TY2_soloLTR"
# indTarClsFeat <- "TY1:MN#-" 
# indTarClsFeat <- "gene:YGR296W#+"
# which(vtUnq == "gene:YGR296W#+")
# indTarClsFeat <- "gene:S288C_G0022800#+"
# which(vtUnq == "gene:S288C_G0022800#+")
# which(vtUnq == "ARS+ARS102#+")
# which(vtUnq == "X_element+TEL01L#+")
# indTarClsFeat <- unique(dtImpgAllSzFlt[, Target_clsfeat])[3]
# indTarClsFeat <- vtUnq[1]
for (indTarClsFeat in vtUnq) {
  
  ## block calculation --------------------------------------------------------
  
  ### dev print(nBlocks)
  
  ### get the target class-features that will generate the block
  dtTarClsFeat <- dtImpgAllSzFlt[indTarClsFeat, on = "Target_clsfeat"]
  ### check if the corresponding interval is still in dtImpgAllSzFlt
  ### (if it was ovelapping with a former indTarClsFeat
  ### it has already been removed, but dtTarClsFeat will still have a line
  ### e.g. the following):
  ### Query_id Query_start Query_end Target_id Target_start Target_end
  ###   <char>       <num>     <num>    <char>        <num>      <num>
  ###     <NA>         NA        NA      <NA>           NA         NA
  ### Target_clsfeat N_score Query_alndir Target_alndir Target_cls Target_len
  ###         <char>  <char>       <char>        <char>     <char>      <num>
  ###    buciodeculo    <NA>         <NA>          <NA>       <NA>         NA
  
  if (any(is.na(dtTarClsFeat[, Query_id]))) {
    next
  }
  
  ### first round of overlap (on query coordinates)
  setkey(dtTarClsFeat, Query_id, Query_start, Query_end)
  dtIndOverOne <- foverlaps(dtImpgAllSzFlt, dtTarClsFeat,
                            nomatch = NULL, type = "any", which = T)
  dtBlock <- dtImpgAllSzFlt[dtIndOverOne[, xid], ]
  print(format(object.size(dtBlock), units = "auto"))
  
  ### second round of overlap (on target coordinates)
  setkey(dtBlock, Target_id, Target_start, Target_end)
  dtIndOverTwo <- foverlaps(dtImpgAllSzFlt, dtBlock, 
                            nomatch = NULL, type = "any", which = T)
  dtBlock <- rbindlist(list(dtBlock, dtImpgAllSzFlt[dtIndOverTwo[, xid], ]))
  print(format(object.size(dtBlock), units = "auto"))
  
  ### block classes
  allClassesInBlock <- unique(dtBlock[, Target_cls])
  
  ### dev
  # cat(indTarClsFeat, "\t",
  #     nBlocks, "\t", length(allClassesInBlock), "\n",
  #     file = "~/Desktop/check-nsb.txt", append = T)
  
  ## sub-blocks reduction (class-by-class, strand-aware) ----------------------
  
  ### dev
  # indSb <- "gene"
  # indSb <- "TY1/TY2_soloLTR"
  # indSb <- "Y_prime_element"
  for (indSb in allClassesInBlock) {
    ### reduce (same as bedtools merge) the rows for each class
    ### on the basis of the query coordinates
    grSblock <- GRanges(seqnames = dtBlock[Target_cls == indSb, Query_id],
                        ranges = IRanges(start = dtBlock[Target_cls == indSb,
                                                         Query_start],
                                         end = dtBlock[Target_cls == indSb,
                                                       Query_end]),
                        mcols = dtBlock[Target_cls == indSb, Target_clsfeat])
    grSblockRed <- reduce(grSblock)
    
    ### vector with the features names and class name
    ### cleaned of the "#+" or "#-" at the end of the string
    ### e.g. "gene:S288C_G0022800" and "gene:YGR296W" 
    vtClsFeat <- sub(pattern = "#[\\+\\-]$",
                     replacement = "",
                     x = unique(grSblock$mcols))
    ### features string
    strFeats <- paste(gsub(pattern = paste0(indSb, ":"), replacement = "",
                           x = vtClsFeat),
                      collapse = ",")
    
    dtSblockRed <- data.table(as.character(grSblockRed@seqnames),
                              as.character(grSblockRed@ranges))
    
    ## transformation of the reduced sub-block --------------------------------
    
    ### e.g. this dtSblockRed:
    ###                 V1              V2
    ###             <char>          <char>
    ###  SGDref#0#chrVII#+ 1084864-1090591
    ###   S288C#0#chrVII#+ 1085116-1090843
    ###
    ### becomes a two-column data-table with columns:
    ###                 SGDref#0                     S288C#0
    ###                   <char>                      <char>
    ### chrVII:1084864-1090591#+    chrVII:1085116-1090843#+
    
    ### dev
    # vtHaplo <- sub("#chr.*", "", dtSblockRed[, V1])
    # tbHaplo <- table(vtHaplo)
    # if (length(which(tbHaplo > 1)) != 0) {
    #   stop("found it!")
    # }
    
    ### transposition and collapsing explained with
    ### another example, this dtSblockRed:
    ###                      V1            V2
    ###                  <char>        <char>
    ###     DBVPG6765#0#chrIV#- 960288-966619
    ###         AFI#0#chrXIII#- 349539-349680
    ###         AFI#0#chrXIII#- 349713-349989
    ###   DBVPG6765#0#chrXIII#- 354903-355179
    ###           AFI#0#chrIV#- 960661-961002
    ###       S288C#0#chrXIII#- 384508-384784
    ###      SGDref#0#chrXIII#- 378732-379008
    ###
    ### must produce this column:
    ### DBVPG6765#0
    ### chrIV:960288-966619#-,chrXIII:354903-355179#-
    
    ### formatting dtSblockRed
    dtSblockRed[, Haplo_id := sub("#chr.*", "", dtSblockRed[, V1])]
    vtCoord <- paste(SplitSubCol(x = dtSblockRed[, V1], n = 3, s = "#"),
                     dtSblockRed[, V2],
                     sep = ":")
    dtSblockRed[, Info_str :=  paste(vtCoord,
                                     SplitSubCol(dtSblockRed[, V1], 4, "#"),
                                     sep = "#")]
    
    ### TODO: retrieve and store the sequences with samtools
    
    ### collapsing with ";" all the Haplo_id elements of a Haplo_id
    dtSblockRedCo <- dtSblockRed[, paste(Info_str, collapse = ";"),
                                 by = Haplo_id]
    
    ### transpose, producing a data-table
    ### we keep it as data-table even if 
    ### qui
    dtTra <- transpose(dtSblockRedCo)
    ### set column names and format
    setnames(dtTra, as.character(dtTra[1, ]))
    dtTra <- dtTra[-1]
    ### add class and features columns
    dtTra[, ':='(class = rep(indSb, .N), features = rep(strFeats, .N))]
    
    ### add the missing columns: not needed 
    ### since rbindlist makes it by default using fill = TRUE
    # newCols <- setdiff(vtStrainHaploHash, colnames(dtTra))
    # dtTra[, (newCols) := lapply(newCols, function(x) NA)]
    # setcolorder(dtTra, vtStrainHaploHash)
    
    ### append with rbindlist matching (default operation) column names
    dtPanFeats <- rbindlist(list(dtPanFeats, dtTra), fill = T)
    
    ### append the generator feature to the output
    cat(rep(indTarClsFeat, nrow(dtTra)), sep = "\n",
        append = T, file = pathGen)
    ### dev print(nrow(dtTra))
  }
  
  ### delete dtIndOverOne and dtIndOverTwo rows in dtImpgAllSzFlt,
  ### the (any(is.na(dtTarClsFeat[, Query_id]))) will avoid 
  ### the loop to break
  indOut <- unique(c(dtIndOverOne[, xid], dtIndOverTwo[, xid]))
  dtImpgAllSzFlt <- dtImpgAllSzFlt[-indOut]
  
  nBlocks <- nBlocks + 1
  
  ### dev if (nBlocks == 20) stop("We did 20 iterations!")
}

### format and write dtPanFeats
leftCols <- c("class", "features")
colOrder <- c(leftCols, setdiff(names(dtPanFeats), leftCols))
setcolorder(dtPanFeats, colOrder)
pathOutPanFeat <- file.path(dirOut, "pan-features.txt")
# fwrite(dtPanFeats, file = pathOutPanFeat, append = F, quote = F, sep = "\t",
#        col.names = T, na = "MA", nThread = 8, buffMB = 1024)
#>Error in fwrite(dtPanFeats, file = pathOutPanFeat, append = F, quote = F,  : 
#> Bad address: '/home/ltattini/prog/graphs/panda/run-11/png/pan-features.txt'
write.table(x = dtPanFeats, file = pathOutPanFeat, append = F, quote = F,
            sep = "\t", col.names = T, row.names = F, na = "MA")
save(dtPanFeats, file = file.path(dirOut, "pan-features.RData"))
### TODO: write dtPanSeqs
