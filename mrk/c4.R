## header ---------------------------------------------------------------------

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
#' # example usage:
#' CountNoNa(c(1, 2, NA, 4, NA))    # returns 3
#' CountNoNa(c("A", "B", NA, "D"))  # returns 3
#' CountNoNa(c(TRUE, NA, NA, TRUE)) # returns 2
#' @export
CountNoNa <- function(x) {
  z <- sum(!is.na(x))
  return(z)
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
#' # example usage:
#' CountSemicSep(c("A;B;C", "D;E", NA, "F")) # returns 6 (3+2+1)
#' CountSemicSep(c(NA, "X;Y;Z", "P;Q", ""))  # returns 5 (3+2+0)
#' CountSemicSep(c(NA, NA))                  # returns 0
#' @export
CountSemicSep <- function(x) {
  z <- sum(lengths(strsplit(x[!is.na(x)], ";")))
  return(z)
}

#' Count Features in a Comma-Separated String
#'
#' This function counts the number of features in a character vector where 
#' values are separated by commas (`,`).
#' If a value does not contain any commas, 
#' it is assumed to be a single feature. 
#' The function correctly handles `NA` values.
#'
#' @param x A character vector where each element may contain 
#'          comma-separated features.
#' @return An integer vector of the same length as `x`, where each value 
#'         represents the number of features 
#'         in the corresponding element of `x`.
#' @examples
#' # example usage:
#' CountAnyFeat(c("A,B,C", "D,E", NA, "F")) # returns c(3, 2, NA, 1)
#' CountAnyFeat(c("X,Y,Z", "P,Q", "R", "")) # returns c(3, 2, 1, 1)
#' CountAnyFeat(c(NA, NA, "A,B"))           # returns c(NA, NA, 2)
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
#' # example usage:
#' CountSysFeat(c("YAL001W,YBR102C", "YDL113C,YGR098W", NA, "YPR202W")) 
#' # returns c(2, 2, NA, 1)
#' 
#' CountSysFeat(c("YAL003W", "YCR008W,YHR209C", "NotAGene", "Blinda")) 
#' # returns c(1, 2, 0, 0)
#' 
#' CountSysFeat(c(NA, NA, "YML075C,YPL197W,YJL112W")) 
#' # returns c(NA, NA, 3)
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
#' # example usage:
#' CountRidFeat(c("ABA_G0053280,CMF_HP1_G0054150", "AAB_G7654321"))
#' # returns c(2, 1)
#' 
#' CountRidFeat(c("DBVPG6765_G0027560,YFL051C,YJL216C", "", "NoMatch"))
#' # returns c(1, 0, 0)
#' 
#' CountRidFeat(c(NA, NA, "AAB_G5555555,AIL_G6666666"))
#' # returns c(NA, NA, 2)
#' @export
CountRidFeat <- function(x) {
  ### the pattern for random ID genes
  ptnRid <- "_G[0-9]{7}"
  
  ### apply the regular expression and count occurrences
  y <- sapply(gregexpr(ptnRid, x, perl = T), function(x) sum(x > 0))
  
  return(y)
}

#' Prepend Column Names to Substrings in Selected Columns of a data.table
#'
#' This function modifies a `data.table` by prepending each value 
#' in the specified columns 
#' with the corresponding column name followed by `#`. 
#' Substrings within a column are assumed 
#' to be separated by semicolons (`;`). `NA` values remain unchanged.
#'
#' @param x A `data.table` containing the dataset to be modified.
#' @param y A numeric vector of column indices specifying 
#'          which columns should be modified.
#'
#' @return The input `data.table` is modified in-place, 
#' with the selected columns updated 
#' to have their column names prefixed to each substring.
#'
#' @examples
#' library(data.table)
#' 
#' # create example data.table
#' dt <- data.table(
#'   ID = 1:3,
#'   SGD#0 = c("chrIV:30-32", "chrXIV:386-541;chrXIV:260-415", NA),
#'   AAB#1 = c("", "chrV:35-80", "chrXII:99-111")
#' )
#' 
#' # define columns to modify 
#' # (assuming SGD#0 is in column 2 and AAB#1 in column 3)
#' indHapCols <- c(2, 3)
#'
#' # apply the function
#' AddColPref(dt, indHapCols)
#'
#' # expected output:
#' #     ID   SGD#0                                        AAB#1
#' # 1:  1    SGD#0#chrIV:30-32                            ""  (unchanged)
#' # 2:  2    SGD#0#chrXIV:386-541;SGD#0chrXIV:260-415     AAB#1#chrV:35-80
#' # 3:  3    NA (unchanged)                               AAB#1#chrXII:99-111
#'
#' @export
AddColPref <- function(x, y) {
  x[, (y) := lapply(y, function(indC) {
    col_name <- names(x)[indC]
    ifelse(!is.na(x[[indC]]),
           gsub("([^;]+)", paste0(col_name, "#\\1"), x[[indC]], perl = T),  
           x[[indC]])})]
  return(x)
}

#' Replace NA Values with Empty Strings in a data.table
#'
#' This function modifies a data.table replacing all NA values 
#' with empty strings ("").
#' You can choose whether to apply this only to character 
#' columns or to all columns
#' (which will coerce all types to character).
#'
#' @param x A `data.table` object to modify in-place.
#' @param y Logical. If `TRUE`, replaces NA in all columns 
#'          (and coerces to character).
#'          If `FALSE` (default), replaces only NA values
#'          in character columns.
#'
#' @return The modified `data.table` with NA values replaced by "".
#' @examples
#' library(data.table)
#' x <- data.table(A = c("apple", NA), B = c(NA, "orange"), C = c(1, NA))
#' ReplaceNAtoEmpty(x)         # only character columns
#' ReplaceNAtoEmpty(x, TRUE)   # all columns, coerced to character
ReplaceNAtoEmpty <- function(x, y = F) {
  if (!data.table::is.data.table(x)) {
    ### stop prints "Error: " by default
    stop("input must be a data.table.")
  }
  colsToMod <- if (y) names(x) else names(x)[sapply(x, is.character)]
  
  x[, (colsToMod) := lapply(.SD, function(x) {
    x[is.na(x)] <- ""
    x
  }), .SDcols = colsToMod]
  
  return(x)
}

#' Replace Empty Strings with NA Values in a data.table
#'
#' This function modifies a data.table by replacing all empty strings ("")
#' with NA values. You can choose whether to apply this only to character 
#' columns or to all columns (which will coerce all types to character).
#'
#' @param x A `data.table` object to modify in-place.
#' @param y Logical. If `TRUE`, replaces "" in all columns 
#'          (and coerces to character).
#'          If `FALSE` (default), replaces only "" values
#'          in character columns.
#'
#' @return The modified `data.table` with "" values replaced by NA.
#' @examples
#' library(data.table)
#' x <- data.table(A = c("apple", ""), B = c("", "orange"), C = c("1", ""))
#' ReplaceEmptyToNA(x)         # only character columns
#' ReplaceEmptyToNA(x, TRUE)   # all columns, coerced to character
ReplaceEmptyToNA <- function(x, y = F) {
  if (!data.table::is.data.table(x)) {
    stop("input must be a data.table.")
  }
  colsToMod <- if (y) names(x) else names(x)[sapply(x, is.character)]
  
  x[, (colsToMod) := lapply(.SD, function(col) {
    col[col == ""] <- NA
    col
  }), .SDcols = colsToMod]
  
  return(x)
}

#' Filter intervals in a semicolon-separated genomic string based on a list
#' of valid ("kept") intervals.
#'
#' @param cell A character string of the form "chr:start-end;chr:start-end;..."
#'             or NA.
#' @param regs A character vector of intervals to retain, in the
#'             same "chr:start-end" format.
#'
#' @return A semicolon-separated string of intervals found in both the cell
#'         and the kept list, or NA if none match.
#' @examples
#' dtRegions <- data.table(
#' "Class_id" = c("gene", "gene", "gene", "gene", "gene"),
#' "Features_id" = c("anello", "del", "conte", "antani", "sblinda"),
#' "SGDref#0" = c(
#'   "chrVII:9161-9695",
#'   "chrVII:9518-9751",
#'   NA,
#'   "chrXIV:1375-2528;chrXIV:36539-36694;chrIV_chrII:36539-37694",
#'   "chrX:1000-1010;chrX:1200-1800"),
#' "BUCIO#1" = c(
#'   "chrV:9161-9295",
#'   NA,
#'   NA,
#'   NA,
#'   "chrX:1000-2010;chrX:1200-1800"),
#' "DECULO#2" = c(
#'   "chrV:9161-9895",
#'   NA,
#'   NA,
#'   NA,
#'   NA)
#' )
FilterCellByKeptIntervals <- function(cell, reg) {
  if (is.na(cell)) return(NA_character_)
  allParts <- unlist(strsplit(cell, ";", fixed = T))
  intvKept <- allParts[allParts %in% reg]
  if (length(intvKept) == 0) return(NA_character_)
  paste(intvKept, collapse = ";")
}

## settings -------------------------------------------------------------------

### reference genome
idRefHap <- "SGDref-0"
### pattern of systematic genes
ptnSys <- "Y[A-P][L,R][0-9]{3}[W,C]"
### pattern of random id genes
ptnRid <- "_G[0-9]{7}"
### the type of sub-blocks to process
strSblock <- "gene"

### fixed settings
dirBase <- dirname(this.dir())

### output
dirOut <- file.path(dirBase, "sts", strSblock)
unlink(dirOut, recursive = T)
dir.create(dirOut, recursive = T)
pathPanGenesHaplo <- file.path(dirOut, "sts-by-haplos.txt")
pathNorefPriv <- file.path(dirOut, "n-rid-private.txt")
pathCountSblocsRegs <- file.path(dirOut, "n-sblocks-regs.txt")
pathPanGenesGenomes <- file.path(dirOut, "sts-by-genomes.txt")

## clmnt ----------------------------------------------------------------------

### script name
myName <- current_filename()
cat("[", myName, "] ",
    "Making the statistics. ",
    "\n", sep = "")

### load the pangenome (dtPanFeats)
pathInPan <- file.path(dirBase, "png", "pan-features.RData")
load(pathInPan)
### dev dtPanFeats <- dtPanFeats[1:20, 1:10, with = F]

### we need the first and last haplotype columns to correctly subset 
### the data-table when we calculate the statistics 
nHaplos <- ncol(dtPanFeats) - 2
indHapCols <- 3:ncol(dtPanFeats)

## filter sub-blocks by regions size ------------------------------------------

dtPanFeatsRaw <- dtPanFeats[Class_id == strSblock]

### columns with genomic intervals
intervalCols <- setdiff(names(dtPanFeatsRaw), c("Class_id", "Features_id"))

### process each row independently
dtPanFeatsGns <- dtPanFeatsRaw[, {
  ### extract all intervals in one row across all genomic columns
  allIntervals <- unlist(tstrsplit(unlist(.SD), ";", fixed = T))
  allIntervals <- allIntervals[!is.na(allIntervals)
                               & grepl(":", allIntervals)]
  
  if (length(allIntervals) == 0) {
    ### if there are no valid intervals
    outCols <- lapply(.SD, function(x) NA_character_)
  } else {
    ### calculate interval lengths
    intervalCoords <- tstrsplit(gsub(".*:", "", allIntervals), "-", fixed = T)
    intervalStarts <- as.numeric(intervalCoords[[1]])
    intervalEnds <- as.numeric(intervalCoords[[2]])
    intervalLengths <- intervalEnds - intervalStarts
    
    ### compute median
    M <- median(intervalLengths, na.rm = T)
    
    ### identify good intervals (length ≥ 0.5 * M)
    indBon <- which(intervalLengths >= 0.5 * M)
    bonIntervals <- allIntervals[indBon]
    
    ### rebuild each genomic column
    outCols <- lapply(.SD, FilterCellByKeptIntervals,
                      reg = bonIntervals)
  }
  
  ### return Class_id and Features_id along with filtered columns
  c(list(Class_id = Class_id, Features_id = Features_id), outCols)
}, by = seq(1:nrow(dtPanFeatsRaw)), .SDcols = intervalCols]

### clean up row index
dtPanFeatsGns[, seq := NULL]

### garbage collection
rm(dtPanFeatsRaw)
invisible(gc())

## sub-blocks statistics ------------------------------------------------------

cat("[", myName, "] ",
    "Counting on sub-blocks labelled as: ",
    strSblock,
    "\n", sep = "")
### dev dtPanFeatsGns[20, 4:10] <- NA

nSBlocks <- nrow(dtPanFeatsGns)
cat("[", myName, "] ",
    "Number of sub-blocks: ",
    nSBlocks,
    "\n", sep = "")

nSys <- nrow(dtPanFeatsGns[grep(ptnSys, Features_id)])
cat("[", myName, "] ",
    "Number of sub-blocks with systematic label(s): ",
    nSys,
    "\n", sep = "")

nLabellable <- nrow(dtPanFeatsGns[grepl(ptnSys, Features_id)
                                  & grepl(ptnRid, Features_id)])
cat("[", myName, "] ",
    "Number of sub-blocks with at least one labellable feature: ",
    nLabellable,
    "\n", sep = "")

nUnl <- nrow(dtPanFeatsGns[grep(ptnSys, Features_id, invert = T)])
cat("[", myName, "] ",
    "Number of sub-blocks with only random label(s): ",
    nUnl,
    "\n", sep = "")

## presence (haplotype-based) -------------------------------------------------

### number of haplotypes with at least one region in the sub-block
dtPanFeatsGns[, Ν_pres := rowSums(!is.na(as.matrix(.SD))),
              .SDcols = indHapCols]
### fraction of haplotypes with at least one region in the sub-block
dtPanFeatsGns[, F_pres := c(Ν_pres / nHaplos) ]

## core and dispensable sub-blocks statistics (haplotype-based) ---------------

### all core
cat("[", myName, "] ",
    "Number of core genes (haplotype-based): ",
    dtPanFeatsGns[, sum(F_pres == 1)],
    "\n", sep = "")
### all but one core
cat("[", myName, "] ",
    "Number of N-1 core genes (haplotype-based): ",
    dtPanFeatsGns[, sum(Ν_pres >= c(nHaplos - 1))],
    "\n", sep = "")
### all but two core
cat("[", myName, "] ",
    "Number of N-2 core genes (haplotype-based): ",
    dtPanFeatsGns[, sum(Ν_pres >= c(nHaplos - 2))],
    "\n", sep = "")
### strictly dispensable: 1 ≤ N presence < all
cat("[", myName, "] ",
    "Number of strictly dispensable genes (haplotype-based): ",
    dtPanFeatsGns[, sum(Ν_pres < nHaplos)],
    "\n", sep = "")
### common dispensable: 1 < N presence < all
cat("[", myName, "] ",
    "Number of common dispensable (haplotype-based): ",
    dtPanFeatsGns[, sum(Ν_pres < nHaplos
                        & 1 < Ν_pres)],
    "\n", sep = "")
### 1-leniently common dispensable: 1 < N presence < all - 1
cat("[", myName, "] ",
    "Number of 1-leniently common dispensable (haplotype-based): ",
    dtPanFeatsGns[, sum(Ν_pres < c(nHaplos - 1)
                        & 1 < Ν_pres)],
    "\n", sep = "")
### 2-leniently common dispensable: 1 < N presence < all - 2
cat("[", myName, "] ",
    "Number of 2-leniently common dispensable (haplotype-based): ",
    dtPanFeatsGns[, sum(Ν_pres < c(nHaplos - 2)
                        & 1 < Ν_pres)],
    "\n", sep = "")
### private
cat("[", myName, "] ",
    "Number of private genes (haplotype-based): ",
    dtPanFeatsGns[, sum(Ν_pres == 1)],
    "\n", sep = "")

## decomposition of the features in each sub-block (haplotype-based) ----------

### count how many features are present in a sub-block with gregexpr
### and using integers only (e.g. 1L), although it does not work with empty 
### strings (but this cannot happen in dtPanFeatsGns[, Features_id])
dtPanFeatsGns[, N_feats := CountAnyFeat(Features_id)]
dtPanFeatsGns[, N_feats_sys := CountSysFeat(Features_id)]
dtPanFeatsGns[, N_feats_rid := CountRidFeat(Features_id)]

### save dtPanFeatsGns with counts
fwrite(file = pathPanGenesHaplo, x = dtPanFeatsGns, sep = "\t",
       row.names = F, col.names = T)

### number of non-reference private features per haplotype
tbNorefPriv <- table(dtPanFeatsGns[Ν_pres == 1 & N_feats_rid == 1,
                                   sub("_(?!.*_).*", "",
                                       Features_id,
                                       perl = T)])
### save
write.table(tbNorefPriv, file = pathNorefPriv, quote = F, sep = "\t",
            row.names = F, col.names = c("Haplo_id", "N_priv"))

### check if N_feats[i] = N_feats_sys[i] + N_feats_rid[i]
dtPanFeatsGns[, Good_sum := ifelse(N_feats_sys + N_feats_rid == N_feats, T, F)]
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

## count sub-blocks and regions (haplotype-based) -----------------------------

dtSb <- dtPanFeatsGns[, lapply(.SD, CountNoNa), 
                      .SDcols = indHapCols]
dtRe <- dtPanFeatsGns[, lapply(.SD, CountSemicSep),
                      .SDcols = indHapCols]
dtN <- rbindlist(list(dtSb, dtRe))
dtT <- transpose(dtN)
dtCounts <- data.table(colnames(dtN),
                       dtT)
colnames(dtCounts) <- c("Haplo_id", "N_sblocks", "N_regions")

### save dtT with counts
fwrite(file = pathCountSblocsRegs, x = dtCounts, sep = "\t",
       row.names = F, col.names = T)

cvSbloc <- sd(dtCounts[, N_sblocks]) / mean(dtCounts[, N_sblocks])
cvRegis <- sd(dtCounts[, N_regions]) / mean(dtCounts[, N_regions])

## collapse haplotypes into genomes -------------------------------------------

dtTmp <- dtPanFeatsGns[, .SD, .SDcols = c(1:2, indHapCols)]
AddColPref(dtTmp, indHapCols)
ReplaceNAtoEmpty(dtTmp)
### remove numeric suffix
vtGnmIds <- gsub("#[0-9]+$", "", names(dtTmp)[indHapCols])
### change column names and identify unique ids
names(dtTmp)[indHapCols] <- vtGnmIds
vtGnmIdsUniq <- unique(vtGnmIds)

### initialize output with genome columns
dtPanFeatsGnsGnm <- dtTmp[, .(Class_id, Features_id)]

### process each unique prefix
for (indP in vtGnmIdsUniq) {
  vtColMatched <- which(names(dtTmp) == indP)
  
  if (length(vtColMatched) > 1) {
    colsBon <- dtTmp[, ..vtColMatched]
    dtPanFeatsGnsGnm[, (indP) := do.call(paste, c(colsBon, sep = ";"))]
    dtPanFeatsGnsGnm[, (indP) := gsub("^;|;$", "", get(indP))]
  } else if (length(vtColMatched) == 1) {
    dtPanFeatsGnsGnm[, (indP) := dtTmp[[vtColMatched]]]
  }
}
ReplaceEmptyToNA(dtPanFeatsGnsGnm)

### garbage collection
rm(dtTmp)
invisible(gc())

## presence (genome-based) ----------------------------------------------------

nGenomes <- ncol(dtPanFeatsGnsGnm) - 2
indGnmCols <- 3:ncol(dtPanFeatsGnsGnm)

### number of genomes with at least one region in the sub-block
dtPanFeatsGnsGnm[, Ν_pres := rowSums(!is.na(as.matrix(.SD))),
                 .SDcols = indGnmCols]
### fraction of genomes with at least one region in the sub-block
dtPanFeatsGnsGnm[, F_pres := c(Ν_pres / nGenomes) ]

## core and dispensable sub-blocks statistics (genome-based) ------------------

### all core
cat("[", myName, "] ",
    "Number of core genes (genome-based): ",
    dtPanFeatsGnsGnm[, sum(F_pres == 1)],
    "\n", sep = "")
### all but one core
cat("[", myName, "] ",
    "Number of N-1 core genes (genome-based): ",
    dtPanFeatsGnsGnm[, sum(Ν_pres >= c(nGenomes - 1))],
    "\n", sep = "")
### all but two core
cat("[", myName, "] ",
    "Number of N-2 core genes (genome-based): ",
    dtPanFeatsGnsGnm[, sum(Ν_pres >= c(nGenomes - 2))],
    "\n", sep = "")
### strictly dispensable: 1 ≤ N presence < all
cat("[", myName, "] ",
    "Number of strictly dispensable genes (genome-based): ",
    dtPanFeatsGnsGnm[, sum(Ν_pres < nGenomes)],
    "\n", sep = "")
### common dispensable: 1 < N presence < all
cat("[", myName, "] ",
    "Number of common dispensable (genome-based): ",
    dtPanFeatsGnsGnm[, sum(Ν_pres < nGenomes
                           & 1 < Ν_pres)],
    "\n", sep = "")
### 1-leniently common dispensable: 1 < N presence < all - 1
cat("[", myName, "] ",
    "Number of 1-leniently common dispensable (genome-based): ",
    dtPanFeatsGnsGnm[, sum(Ν_pres < c(nGenomes - 1)
                           & 1 < Ν_pres)],
    "\n", sep = "")
### 2-leniently common dispensable: 1 < N presence < all - 2
cat("[", myName, "] ",
    "Number of 2-leniently common dispensable (genome-based): ",
    dtPanFeatsGnsGnm[, sum(Ν_pres < c(nGenomes - 2)
                           & 1 < Ν_pres)],
    "\n", sep = "")
### private
cat("[", myName, "] ",
    "Number of private genes (genome-based): ",
    dtPanFeatsGnsGnm[, sum(Ν_pres == 1)],
    "\n", sep = "")

## decomposition of the features in each sub-block (genome-based) -------------

### count how many features are present in a sub-block with gregexpr
### and using integers only (e.g. 1L), although it does not work with empty 
### strings (but this cannot happen in dtPanFeatsGns[, Features_id])
dtPanFeatsGnsGnm[, N_feats := CountAnyFeat(Features_id)]
dtPanFeatsGnsGnm[, N_feats_sys := CountSysFeat(Features_id)]
dtPanFeatsGnsGnm[, N_feats_rid := CountRidFeat(Features_id)]

### all core systematic: must be equal to all core
cat("[", myName, "] ",
    "Number of core systematic sub-blocks (genome-based): ",
    dtPanFeatsGnsGnm[, sum(F_pres == 1
                           & N_feats_sys != 0)],
    "\n", sep = "")
### strictly dispensable systematic: 1 ≤ N presence < all
cat("[", myName, "] ",
    "Number of strictly dispensable systematic sub-blocks (genome-based): ",
    dtPanFeatsGnsGnm[, sum(Ν_pres < nGenomes
                           & N_feats_sys != 0)],
    "\n", sep = "")

### save dtPanFeatsGnsGnm with counts
fwrite(file = pathPanGenesGenomes, x = dtPanFeatsGnsGnm, sep = "\t",
       row.names = F, col.names = T)
