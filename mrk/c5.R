## header ---------------------------------------------------------------------

options(scipen = 999)
options(stringsAsFactors = F)
rm(list = ls())
library(data.table)
library(this.path)
library(scriptName)

## functions ------------------------------------------------------------------

#' Count Systematic Gene Features in a Character Vector
#'
#' This function counts the number of occurrences of systematic gene features 
#' in a character vector based on a predefined regular expression pattern.
#' The pattern follows the format: `Y[A-P][L,R][0-9]{3}[W,C](-[A-Z])?`, 
#' which is used for systematic yeast gene names.
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
  ptnSys <- "Y[A-P][L,R][0-9]{3}[W,C](-[A-Z])?"
  
  ### apply the regular expression and count occurrences
  y <- sapply(gregexpr(ptnSys, x, perl = T), function(x) sum(x > 0))
  
  return(y)
}

#' Extract Systematic Gene Names from a Character Vector
#'
#' This function extracts all occurrences of systematic yeast gene names
#' from each element of a character vector. Systematic gene names follow
#' the pattern \code{"Y[A-P][L,R][0-9]{3}[W,C](-[A-Z])?"}, which corresponds
#' to canonical ORF identifiers such as \code{"YAL001W"} or \code{"YBR102C-A"}.
#'
#' The function returns a list where each element corresponds to the
#' extracted gene names from the corresponding element of the input vector.
#' Elements with no matches return a character vector of length zero.
#'
#' @param x A character vector. Each element may contain one or more 
#'          systematic gene names embedded in arbitrary text.
#'
#' @return A list of character vectors. Each list element contains
#'         the systematic gene names found in the corresponding element
#'         of \code{x}. If a vector element contains no matches, an
#'         empty character vector is returned. If \code{x} contains
#'         \code{NA} values, corresponding outputs are also \code{NA}.
#'
#' @examples
#' ExtractSys(c("gene=YAL001W;other=foo",
#'              "ID=YBR102C,YBR102C-A;desc",
#'              "no_sys_gene_here",
#'              NA,
#'              "YGR098W text YPL160C"))
#'
#' # returns:
#' # [[1]] "YAL001W"
#' # [[2]] "YBR102C" "YBR102C-A"
#' # [[3]] character(0)
#' # [[4]] NA
#' # [[5]] "YGR098W" "YPL160C"
#'
#' @export
ExtractSys <- function(x) {
  ### the pattern for systematic yeast gene names
  ptnSys <- "Y[A-P][L,R][0-9]{3}[W,C](-[A-Z])?"
  
  lsStr <- regmatches(x, gregexpr(ptnSys, x, perl = T))
  return(lsStr)
}

## settings -------------------------------------------------------------------

### fixed settings
dirBase <- dirname(this.dir())
dirAnnoGff <- file.path(dirBase, "anno", "gff")
hdGff <- c("Chr_id", "Strain_id", "Feat_type", "S_coord",
           "E_coord", "S_val", "Strand_id", "Frame_id", "Attribute_str")
ptnSys <- "Y[A-P][L,R][0-9]{3}[W,C](-[A-Z])?"

### input
strSblock <- "gene"
dirBlocks <- file.path(dirBase, "sts", strSblock)
pathPanGenesHaplo <- file.path(dirBlocks, "sts-by-haplos.txt")

### skipped columns
skippedCols <- c("Class_id", "Features_id", "SGDref#0", "Ν_pres", "F_pres",
                 "N_feats", "N_feats_sys", "N_feats_rid")

### output
dirRes <- file.path(dirBase, "check-sys")
unlink(dirRes, recursive = T)
dir.create(dirRes, recursive = T)
pathSysSummary <- file.path(dirRes, "sys-summary.txt")

## clmnt ----------------------------------------------------------------------

### script name
myName <- current_filename()
cat("[", myName, "] ",
    "Checking systematic genes. ",
    "\n", sep = "")

### read sub-blocks
if (!exists("dtPanFeatsGns")) {
  dtPanFeatsGns <- fread(pathPanGenesHaplo)
}

### columns to be checked
idCols <- setdiff(names(dtPanFeatsGns), skippedCols)
### precompute logical index for N_feats_sys > 0
idxSys <- !is.na(dtPanFeatsGns[["N_feats_sys"]]) & dtPanFeatsGns[["N_feats_sys"]] > 0

### initialise summary table
nId <- length(idCols)
dtSummary <- data.table(Genome_id = idCols,
                        N_sys_pan = as.integer(NA),
                        N_sys_gff = as.integer(NA),
                        N_lost = as.integer(NA),
                        N_gained = as.integer(NA))

### build summary of systematic genes
### dev: indC <- "CMF#1"
for (indR in 1:length(idCols)) {
  indC <- idCols[indR]
  ### we make a little magheggio to make
  ### "YAL068C/YBR301W/.../YOL161C" count as 1 (a case we see in the gff)
  ### "YCL042W,YCL040W" count as 1 (a case we see in the pangenome)
  ### e.g.:
  ### x <- c("YAL001W/YBR102C", "YCL042W,YCL040W", NA, "YPR202W")
  ### length(CountSysFeat(x[!is.na(x)]))
  ### > 3
  ### while with:
  ### sum(CountSysFeat(x[!is.na(x)]))
  ### > 5
  
  ### count sys genes in the gff annotation
  indD <- sub(pattern = "#", replacement = "-", x = indC)
  pathAnnoGff <- list.files(dirAnnoGff, pattern = indD, full.names = T)
  if (length(pathAnnoGff) != 1) {
    cat("[", myName, "] ",
        "Found multiple annotations files for ", indD, ".\n",
        "Skipping ", indD, ".\n",
        sep = "")
    next
  }
  dtGff <- fread(file = pathAnnoGff,
                 sep = "\t",
                 header = F,
                 verbose = F,
                 select = c(3L, 9L),
                 col.names = c("Feat_type", "Attribute_str"))
  vtGeneAttribute <- dtGff[Feat_type == "gene", Attribute_str]
  nSysGff <- length(CountSysFeat(vtGeneAttribute))
  
  ### count sys genes in the pangenome annotations
  vtFeatures <- dtPanFeatsGns[get(indC) != "" & idxSys, Features_id]
  nSysPan <- length(CountSysFeat(vtFeatures))

  ### compare the two sets
  strAttributes <- unique(unlist(ExtractSys(vtGeneAttribute), use.names = F))
  strFeatures <- unique(unlist(ExtractSys(vtFeatures)), use.names = F)
  geneLost <- setdiff(strAttributes, strFeatures)
  geneGained <- setdiff(strFeatures, strAttributes)
  nGeneLost <- length(geneLost)
  nGeneGained <- length(geneGained)
  
  ### write gene lists
  fwrite(as.data.table(geneLost),
         file = file.path(dirRes, paste0("gene-lost-", indD, ".txt")),
         row.names = F, col.names = F)
  fwrite(as.data.table(geneGained),
         file = file.path(dirRes, paste0("gene-gained-", indD, ".txt")),
         row.names = F, col.names = F)
  
  ### summary
  dtSummary[indR, `:=`(N_sys_pan = nSysPan,
                       N_sys_gff = nSysGff,
                       N_lost = nGeneLost,
                       N_gained = nGeneGained)]
}

### write summary table
fwrite(dtSummary, pathSysSummary, sep = "\t")
cat("[", myName, "] ",
    "Summary of systematic genes written to: ",
    pathSysSummary, "\n", sep = "")
