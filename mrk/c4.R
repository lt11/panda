options(scipen = 999)
options(stringsAsFactors = F)
rm(list = ls())
library(data.table)
library(this.path)
library(scriptName)

## settings -------------------------------------------------------------------

### pattern of systematic genes
ptnSys <- "Y[A-P][L,R][0-9]{3}[W,C]"
### pattern of random id genes
ptnRid <- "_G[0-9]{7}"

### fixed settings
dirBase <- dirname(this.dir())

### output
dirOut <- file.path(dirBase, "sts")
unlink(dirOut, recursive = T)
dir.create(dirOut, recursive = T)
pathPanGenes <- file.path(dirOut, "sts-genes.txt")

## clmnt ----------------------------------------------------------------------

### script name
myName <- current_filename()
cat("[", myName, "] ",
    "Making the pangenome. ",
    "\n", sep = "")

### load SGDref features
pathBedRef <- "/Users/Lorenzo/prog/graphs/panda/sc-phhd/run-1/anno/bed/SGDref-0.bed"
dtBedRef <- fread(file = pathBedRef, sep = "\t", header = F)
names(dtBedRef) <- c("Strain_hap_chr", "S_coord", "E_coord", "Class_feat_strand")
dtBedRefGns <- dtBedRef[grep("^gene", Class_feat_strand)]

### load the generators
pathInGrt <- "/Users/Lorenzo/prog/graphs/panda/sc-phhd/run-1/png/generators.txt"
dtGrt <- fread(file = pathInGrt,
               sep = ":", header = T)
dtGrt[, c("Feature_id", "Strand_id") := 
        tstrsplit(.SD[[which(colnames(dtGrt) == "feature#strand")]], "#")]
setnames(dtGrt, old = which(colnames(dtGrt) == "class"), new = "Class_id")

### load the pangenome (dtPanFeats)
pathInPan <- "/Users/Lorenzo/prog/graphs/panda/sc-phhd/run-1/png/pan-features.RData"
load(pathInPan)
# setnames(dtPanFeats, old = names(dtPanFeats)[1:2],
#          new = c("Class_id", "Features_id"))

### strain column boundaries
nHaplos <- ncol(dtPanFeats) - 2
nSkip <- 2
firstStrainCol <- 1 + nSkip
lastStrainCol <- ncol(dtPanFeats)

### sub-blocks statistics
dtPanFeatsGns <- dtPanFeats[Class_id == "gene"]
nSys <- nrow(dtPanFeatsGns[grep(ptnSys, Features_id)])
nUnl <- nrow(dtPanFeatsGns[grep(ptnSys, Features_id, invert = T)])

### sub-blocks with at least one labellable gene
nLabellable <- nrow(dtPanFeatsGns[grepl(ptnSys, Features_id)
                                  & grepl(ptnRid, Features_id)])

### fraction of presence: how many strains have
### at least on path though a sub-block
# dtPanFeatsGns[, Ν_notna := Reduce(`+`, lapply(.SD, Negate(is.na))),
#               .SDcols = 3:nHaplos]


# cucco <- dtPanFeatsGns[5, as.character(unlist(.SD)), .SDcols = 3:109]
# sum(!is.na(cucco))

### number of strains with at least one region in the sub-block
dtPanFeatsGns[, Ν_pres := rowSums(!is.na(as.matrix(.SD))),
              .SDcols = firstStrainCol:lastStrainCol]
### fraction of strains with at least one region in the sub-block
dtPanFeatsGns[, F_pres := c(Ν_pres / nHaplos) ]

## core and accessory sub-blocks statistics -----------------------------------

### all
dtPanFeatsGns[, sum(F_pres == 1)]
### all but one
dtPanFeatsGns[, sum(Ν_pres >= c(nHaplos - 1))]
### all but two
dtPanFeatsGns[, sum(Ν_pres >= c(nHaplos - 2))]
### private
dtPanFeatsGns[, sum(Ν_pres == 1)]

# the numbers of panda: pangenome annotations data analysis
# we have 6579 systematic genes in SGD
# the gene pangenome (107 haplotypes) has 9433 entries
# an entry is an homology sub-block aka a unique gene ("gene" for short in the next lines)
# 5549 entries have at least one systematic name (some of the 6589 SGD genes are overlapping or copies at different loci, that in the panda output are reported in the same sub-block)
# among these 5549 entries with a systematic name, 1178 have at least one unlabelled gene (so it can get a systematic name)
# 3884 genes do not have a systematic label
# 
# 1625 genes are private
# 4470 are core (for N-2 haplotypes)

## sub-blocks decomposition ---------------------------------------------------

### count how many features are present in a sub-block with gregexpr
### and using integers only (e.g. 1L), altough it does not work with empty 
### strings (but this cannot happen in dtPanFeatsGns[, Features_id])
### gregexpr(",", c("dst", "Sto,ca"))
### [[1]]
### [1] -1
### attr(,"match.length")
### [1] -1
### attr(,"index.type")
### [1] "chars"
### attr(,"useBytes")
### [1] TRUE
### [[2]]
### [1] 4
### attr(,"match.length")
### [1] 1
### attr(,"index.type")
### [1] "chars"
### attr(,"useBytes")
### [1] TRUE
dtPanFeatsGns[, N_feats := sapply(gregexpr(",", Features_id),
                                  function(x) ifelse(x[1] == -1, 1L,
                                                     length(x) + 1L))]
dtPanFeatsGns[, N_feats_sys := sapply(gregexpr(ptnSys,Features_id,
                                               perl = T),
                                      function(x) sum(x > 0))]
dtPanFeatsGns[, N_feats_rid := sapply(gregexpr(ptnRid, Features_id,
                                               perl = T),
                                      function(x) sum(x > 0))]
### save dtPanFeatsGns with counts
fwrite(file = pathPanGenes, x = dtPanFeatsGns, sep = "\t",
       row.names = F, col.names = T)

### check if N_feats[i] = N_feats_sys[i] + N_feats_rid[i]
dtPanFeatsGns[, Good_sum := ifelse(N_feats_sys + N_feats_rid == N_feats,
                                   T, F)]
nTrue <- dtPanFeatsGns[, sum(Good_sum)]
if (nTrue != nrow(dtPanFeatsGns)) {
  cat("[", myName, "] ",
      "N_feats[i] != N_feats_sys[i] + N_feats_rid[i]\n", "\n", sep = "")
  cat("[", myName, "] ",
      "Check i = ", paste(dtPanFeatsGns[Good_sum == F, .I], collapse = ", "),
      "\n", sep = "")
  ### stop prints "Error: " by default
  stop("features patterns probably did not work.")
}

