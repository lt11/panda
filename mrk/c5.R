## header ---------------------------------------------------------------------

### This script requires "samtools faidx". Given: a list of "essential" genes, a folder with the genomes used to 
### build the graph, and the corresponding pangenome, here we extract 
### the sequence encoded in each element or semicolon delimited sub-element 
### of the pangenome table (dtPanFeatsGns) from the corresponding FASTA file.
### We verify that: (1) the sequence starts with a START codon, 
### (2) the sequence ends wit han END codon, and (3) that it does not contain 
### other END codons. If the three conditions are valid, we append the 
### sequence in a row-specific file in FASTA format (specifying e.g. 
### ">SGDref#chrIV:30882-32244" by combining the name of the column from 
### which the element or the sub-element has been extracted with the 
### element or the sub-element itself). If one the conditions is not met, 
### calculate the reverse-complement of the sequence and re-apply 
### the three conditions to the reverse-complement sequence. 
### If the three conditions are valid append the reverse-complement 
### of the sequence in a row-specific file in FASTA format 
### (specifying e.g. ">SGDref#chrIV#RC:30882-32244").
### Otherwise, mark the element or the sub-element as "invalid" 
### by modifying it in the pangenome table e.g. as inv#chrXII:490405-490594 
### and counting the number of invalid elements per row.

options(scipen = 999)
options(stringsAsFactors = F)
rm(list = ls())
library(data.table)
library(this.path)
library(scriptName)
library(stringr)
library(Biostrings)

## function(s) ----------------------------------------------------------------

GetSequence <- function(fa_path, region) {
  cmd <- sprintf("samtools faidx %s '%s'", shQuote(fa_path), region)
  res <- system(cmd, intern = TRUE)
  if (length(res) <= 1) return(NA)
  paste(res[-1], collapse = "")
}

IsValidCoding <- function(seq) {
  strStartCod <- c("ATG")
  strStopCods <- c("TAA", "TAG", "TGA")
  if (nchar(seq) %% 3 != 0) return(FALSE)
  codons <- str_sub(seq, seq(1, nchar(seq), 3), seq(3, nchar(seq), 3))
  if (!(codons[1] %in% strStartCod)) return(FALSE)
  if (!(codons[length(codons)] %in% strStopCods)) return(FALSE)
  if (any(codons[-length(codons)] %in% strStopCods)) return(FALSE)
  return(TRUE)
}

## settings -------------------------------------------------------------------

### reference genome
idRefHap <- "SGDref-0"
### minimal frequency of a gene to be labelled as core
freqCore <- 0.95
### genome path base
dirGenomes <- file.path(Sys.getenv("HOME"), "data",
                        "nano-assemblies-pansn-2024", "genomes")

### fixed settings
dirBase <- dirname(this.dir())
### the type of sub-blocks to process
strSblock <- "gene"
### input
dirBlocks <- file.path(dirBase, "sts", strSblock)
pathPanGenesHaplo <- file.path(dirBlocks, "sts-by-haplos.txt")
pathEsse <- file.path(Sys.getenv("HOME"), "data",
                      "SGD-essential-genes", "essentiality.txt")
### pattern of systematic genes
ptnSys <- "Y[A-P][L,R][0-9]{3}[W,C]"
### pattern of random id genes
ptnRid <- "_G[0-9]{7}"
### skipped columns
skippedCols <- c("Class_id", "Features_id", "Ν_pres", "F_pres", "N_feats",
                 "N_feats_sys", "N_feats_rid", "Sblock_type" )
### output folders
dirSeqGenes <- file.path(dirBase, "seqs", strSblock)
unlink(dirSeqGenes, recursive = T)
dir.create(dirSeqGenes, recursive = T)

## clmnt ----------------------------------------------------------------------

### script name
myName <- current_filename()
cat("[", myName, "] ",
    "Checking selection patterns. ",
    "\n", sep = "")

### read sub-blocks
if (!exists("dtPanFeatsGns")) {
  dtPanFeatsGns <- fread(pathPanGenesHaplo)
}

### dev dtPanFeatsGns <- dtPanFeatsGns[1:100, ]

### read essential non-evolvable (esne) genes
dtEsse <- fread(pathEsse, sep = "\t", col.names = c("Gene_id", "Gene_type"))
### get the list of non-evolvable genes
vtEsne <- dtEsse[Gene_type == "non-evolvable", unique(Gene_id)]

## data preparation -----------------------------------------------------------

### annotate sub-blocks as non-evolvable, core, or dispensable
dtPanFeatsGns[, Sblock_type := fifelse(
  grepl(paste(vtEsne, collapse = "|"), Features_id),
  "non-evolvable",
  fifelse(F_pres > freqCore, "core", "dispensable")
)]

## data processing ------------------------------------------------------------

### process each row
dtPanFeatsGns[, c("N_total", "N_invalid") := .(0, 0)]

### dev
indR <- 5
for (indR in seq(1:nrow(dtPanFeatsGns))) {
  lsRow <- dtPanFeatsGns[indR]
  firstFeat <- sub(",.*", "", lsRow$Features_id)
  sblock_file <- file.path(dirSeqGenes, paste0(lsRow$Sblock_type, "-", firstFeat, ".fa"))
  bonCols <- setdiff(names(lsRow), skippedCols)
  
  ### dev
  col_name <- bonCols[3]
  for (col_name in bonCols) {
    if (!grepl("#", col_name)) next
    elements <- unlist(strsplit(lsRow[[col_name]], ";"))
    if (length(elements) == 0 || all(is.na(elements))) next
    
    ### dev
    el <- elements[1]
    for (el in elements) {
      if (is.na(el) || el == "") next
      lsRow$N_total <- lsRow$N_total + 1
      strGenome <- paste0(sub("#", "-", col_name), "-genome.fa.gz$")
      pathGenome <- list.files(file.path(dirGenomes),
                               recursive = T,
                               full.names = T,
                               pattern = strGenome)
      if (length(pathGenome) != 1) {
        cat("[", myName, "] ",
            "Check the folder with the genomes.",
            "\n", sep = "")
        ### stop prints "Error: " by default
        stop("multiple genomes were found.")
        
      }
      seq <- GetSequence(pathGenome, el)
      hdFasta <- paste0(col_name, "#", el)
      
      if (!is.na(seq) && IsValidCoding(seq)) {
        cat(sprintf(">%s\n%s\n", hdFasta, seq), file = sblock_file, append = T)
      } else {
        rc_seq <- as.character(reverseComplement(DNAString(seq)))
        hdFasta_rc <- paste0(col_name, "#RC#", el)
        if (!is.na(rc_seq) && IsValidCoding(rc_seq)) {
          cat(sprintf(">%s\n%s\n", hdFasta_rc, rc_seq), file = sblock_file, append = T)
        } else {
          el_hdFastaged <- paste0("inv#", el)
          lsRow[[col_name]] <- gsub(el, el_hdFastaged, lsRow[[col_name]], fixed = T)
          lsRow$N_invalid <- lsRow$N_invalid + 1
        }
      }
    }
  }
  dtPanFeatsGns[indR, c("N_total", "N_invalid") := list(lsRow$N_total, lsRow$N_invalid)]
}