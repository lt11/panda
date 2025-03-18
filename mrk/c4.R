options(scipen = 999)
options(stringsAsFactors = F)
rm(list = ls())
library(data.table)
library(this.path)
library(scriptName)

## function(s) ----------------------------------------------------------------

#' Count Non-NA Values in a Vector
#'
#' This function takes a numeric, character, or logical vector and counts 
#' the number of non-NA (non-missing) values in it.
#'
#' @param x A vector of any type (numeric, character, or logical) 
#'          from which non-NA values will be counted.
#' @return An integer representing the number of non-NA values in `x`.
#' @examples
#' # Example usage:
#' CountNoNa(c(1, 2, NA, 4, NA))   # Returns 3
#' CountNoNa(c("A", "B", NA, "D")) # Returns 3
#' CountNoNa(c(TRUE, NA, FALSE, TRUE)) # Returns 3
#' @export
CountNoNa <- function(x) {
  sum(!is.na(x))
}

#' Count Semicolon-Separated Sub-strings in a Vector
#'
#' This function counts the total number of sub-strings 
#' in a string vector where 
#' elements are made of sub-strings separated by semicolons (`;`). 
#' It first removes `NA` values, 
#' then splits each remaining string by `;`, 
#' and finally sums up the total 
#' number of separated elements.
#'
#' @param x A character vector where elements may contain 
#'          semicolon-separated values.
#' @return An integer representing the total number of elements after splitting 
#'         non-NA values by semicolons.
#' @examples
#' # Example usage:
#' CountSemicSep(c("A;B;C", "D;E", NA, "F")) # Returns 6 (3+2+1)
#' CountSemicSep(c(NA, "X;Y;Z", "P;Q", ""))  # Returns 5 (3+2+0)
#' CountSemicSep(c(NA, NA))                  # Returns 0
#' @export
CountSemicSep <- function(x) {
  sum(lengths(strsplit(x[!is.na(x)], ";")))
}

#' Count Features in a Comma-Separated String
#'
#' This function counts the number of features in a character vector where 
#' values are separated by commas (`,`). If a value does not contain any commas, 
#' it is assumed to be a single feature. 
#' The function correctly handles `NA` values.
#'
#' @param x A character vector where each element may contain 
#'          comma-separated features.
#' @return An integer vector of the same length as `x`, where each value 
#'         represents the number of features 
#'         in the corresponding element of `x`.
#' @examples
#' # Example usage:
#' CountAnyFeat(c("A,B,C", "D,E", NA, "F")) # Returns c(3, 2, NA, 1)
#' CountAnyFeat(c("X,Y,Z", "P,Q", "R", "")) # Returns c(3, 2, 1, 1)
#' CountAnyFeat(c(NA, NA, "A,B"))           # Returns c(NA, NA, 2)
#' @export
CountAnyFeat <- function(x) {
  z <- sapply(gregexpr(",", x),
              function(y) ifelse(y[1] == -1, 1L, length(y) + 1L))
  return(z)
}

#' Count Systematic Gene Features in a Character Vector
#'
#' This function counts the number of occurrences of systematic gene features 
#' in a character vector based on a predefined regular expression pattern.
#' The pattern follows the format: `Y[A-P][L,R][0-9]{3}[W,C]`, 
#' which is commonly used for systematic yeast gene names.
#'
#' @param x A character vector where each element 
#'          may contain systematic gene names.
#' @return An integer vector of the same length as `x`, 
#'         where each value represents 
#'         the number of systematic gene matches 
#'         in the corresponding element of `x`.
#' @examples
#' # Example usage:
#' CountSysFeat(c("YAL001W,YBR102C", "YDL113C,YGR098W", NA, "YPR202W")) 
#' # Returns c(2, 2, NA, 1)
#' 
#' CountSysFeat(c("YAL003W", "YCR008W,YHR209C", "NotAGene", "Blinda")) 
#' # Returns c(1, 2, 0, 0)
#' 
#' CountSysFeat(c(NA, NA, "YML075C,YPL197W,YJL112W")) 
#' # Returns c(NA, NA, 3)
#' @export
CountSysFeat <- function(x) {
  ### the pattern for systematic yeast gene names
  ptnSys <- "Y[A-P][L,R][0-9]{3}[W,C]"
  
  ### apply the regular expression and count occurrences
  y <- sapply(gregexpr(ptnSys, x, perl = T), function(x) sum(x > 0))
  
  return(y)
}

#' Count Random Identifier Gene Features in a Character Vector
#'
#' This function counts the number of occurrences of genes with a random id 
#' in a character vector based on a predefined regular expression pattern.
#' The pattern follows the format: `"_G[0-9]{7}"`, which identifies gene ids 
#' starting with the stain name, followed by `_G`, 
#' followed by exactly seven digits.
#'
#' @param x A character vector where each element 
#'          may contain random ID gene names.
#' @return An integer vector of the same length as `x`, 
#'         where each value represents 
#'         the number of random ID gene matches 
#'         in the corresponding element of `x`.
#' @examples
#' # Example usage:
#' CountRidFeat(c("ABA_G0053280,CMF_HP1_G0054150", "AAB_G7654321"))
#' # Returns c(2, 1)
#' 
#' CountRidFeat(c("DBVPG6765_G0027560,YFL051C,YJL216C", "", "NoMatch"))
#' # Returns c(1, 0, 0)
#' 
#' CountRidFeat(c(NA, NA, "AAB_G5555555,AIL_G6666666"))
#' # Returns c(NA, NA, 2)
#' @export
CountRidFeat <- function(x) {
  ### the pattern for random ID genes
  ptnRid <- "_G[0-9]{7}"
  
  ### apply the regular expression and count occurrences
  y <- sapply(gregexpr(ptnRid, x, perl = T), function(x) sum(x > 0))
  
  return(y)
}

## settings -------------------------------------------------------------------

### reference genome
idRefHap <- "SGDref-0"
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
pathCountSblocsRegs <- file.path(dirOut, "n-sblocs-regs.txt")

## clmnt ----------------------------------------------------------------------

### script name
myName <- current_filename()
cat("[", myName, "] ",
    "Making the genes statistics. ",
    "\n", sep = "")

### load SGDref features
pathBedRef <- file.path(dirBase, "anno", "bed", paste0(idRefHap, ".bed"))
dtBedRef <- fread(file = pathBedRef, sep = "\t", header = F)
names(dtBedRef) <- c("Strain_hap_chr", "S_coord", "E_coord", "Class_feat_strand")
dtBedRefGns <- dtBedRef[grep("^gene", Class_feat_strand)]

# ### load the generators
# pathInGrt <- file.path(dirBase, "png", "generators.txt")
# dtGrt <- fread(file = pathInGrt,
#                sep = ":", header = T)
# dtGrt[, c("Feature_id", "Strand_id") := 
#         tstrsplit(.SD[[which(colnames(dtGrt) == "feature#strand")]], "#")]
# setnames(dtGrt, old = which(colnames(dtGrt) == "class"), new = "Class_id")

### load the pangenome (dtPanFeats)
pathInPan <- file.path(dirBase, "png", "pan-features.RData")
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

## core and dispensable sub-blocks statistics (haplotype-based) ---------------

### all core
cat("[", myName, "] ",
    "Number of all-strains core genes: ",
    dtPanFeatsGns[, sum(F_pres == 1)],
    "\n", sep = "")
### all but one core
cat("[", myName, "] ",
    "Number of N-1 core genes: ",
    dtPanFeatsGns[, sum(Ν_pres >= c(nHaplos - 1))],
    "\n", sep = "")
### all but two core
cat("[", myName, "] ",
    "Number of N-2 core genes: ",
    dtPanFeatsGns[, sum(Ν_pres >= c(nHaplos - 2))],
    "\n", sep = "")
### strictly dispensable: 1 ≤ N presence < all
cat("[", myName, "] ",
    "Number of strictly dispensable: ",
    dtPanFeatsGns[, sum(Ν_pres < nHaplos)],
    "\n", sep = "")
### moderately dispensable: 1 < N presence < all
cat("[", myName, "] ",
    "Number of strictly dispensable: ",
    dtPanFeatsGns[, sum(Ν_pres < nHaplos
                        & 1 < Ν_pres)],
    "\n", sep = "")
### 1-leniently dispensable: 1 < N presence < all - 1
cat("[", myName, "] ",
    "Number of 1-leniently dispensable: ",
    dtPanFeatsGns[, sum(Ν_pres < c(nHaplos - 1)
                        & 1 < Ν_pres)],
    "\n", sep = "")
### 2-leniently dispensable: 1 < N presence < all - 2
cat("[", myName, "] ",
    "Number of 1-leniently dispensable: ",
    dtPanFeatsGns[, sum(Ν_pres < c(nHaplos - 2)
                        & 1 < Ν_pres)],
    "\n", sep = "")
### private
cat("[", myName, "] ",
    "Number of private genes: ", dtPanFeatsGns[, sum(Ν_pres == 1)],
    "\n", sep = "")



## sub-blocks decomposition ---------------------------------------------------

### count how many features are present in a sub-block with gregexpr
### and using integers only (e.g. 1L), although it does not work with empty 
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

dtPanFeatsGns[, N_feats := CountAnyFeat(Features_id)]
dtPanFeatsGns[, N_feats_sys := CountSysFeat(Features_id)]
dtPanFeatsGns[, N_feats_rid := CountRidFeat(Features_id)]

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

## count sub-blocks and regions per strain ------------------------------------

dtSb <- dtPanFeatsGns[, lapply(.SD, CountNoNa), 
                      .SDcols = firstStrainCol:lastStrainCol]

dtRe <- dtPanFeatsGns[, lapply(.SD, CountSemicSep),
                      .SDcols = firstStrainCol:lastStrainCol]

dtN <- rbindlist(list(dtSb, dtRe))
dtT <- transpose(dtN)
dtCounts <- data.table(colnames(dtN),
                       dtT)
colnames(dtCounts) <- c("Id_strain", "N_sblocks", "N_regions")

### save dtT with counts
fwrite(file = pathCountSblocsRegs, x = dtCounts, sep = "\t",
       row.names = F, col.names = T)

cvSbloc <- sd(dtCounts[, N_sblocks]) / mean(dtCounts[, N_sblocks])
cvRegis <- sd(dtCounts[, N_regions]) / mean(dtCounts[, N_regions])
