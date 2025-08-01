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
### and counting the number of invalid allIntervals per row.

options(scipen = 999)
options(stringsAsFactors = F)
rm(list = ls())
library(data.table)
library(this.path)
library(scriptName)
library(stringr)
library(Biostrings)

## function(s) ----------------------------------------------------------------

#' Extract a genomic sequence from a compressed and indexed FASTA file
#' using samtools
#'
#' This function uses `samtools faidx` to extract a DNA sequence from a 
#' compressed and indexed FASTA file.
#' It returns the sequence as a plain character string, excluding the 
#' FASTA header line.
#'
#' @param x A character string specifying 
#'          the path to the FASTA file (e.g., `"genome.fa.gz"`).
#' @param y A character string specifying the genomic region 
#'          in samtools format (e.g., `"chrIV:0-200"`).
#'
#' @return A character string containing the nucleotide sequence, 
#'         or `NA` if the region is invalid or not found.
#'
#' @details The function assumes that `samtools` is installed and available 
#'          in the system path.
#'          It reads the output of the `samtools faidx` command, 
#'          removes the header line, and concatenates all sequence lines 
#'          into a single string.
#'
#' @examples
#' \dontrun{
#' GetSequence("SGD-genome.fa.gz", "chrIV:0-200")
#' }
#'
#' @export
GetSequence <- function(x, y) {
  cmd <- sprintf("samtools faidx %s '%s'", shQuote(x), y)
  res <- system(cmd, intern = TRUE)
  if (length(res) <= 1) return(NA)
  paste(res[-1], collapse = "")
  ### return the sequence but not the header
  return(res[2])
}

#' Check if a DNA sequence is a valid protein-coding sequence
#'
#' This function verifies whether a given DNA sequence satisfies basic criteria
#' for being a valid open reading frame (ORF) in protein-coding genes.
#'
#' @param seq A character string representing a DNA sequence 
#'            (composed of A, T, G, and C).
#'
#' @return `TRUE` if the sequence:
#' \itemize{
#'   \item Has a length that is a multiple of 3 (i.e., full codons),
#'   \item Starts with a start codon (`ATG`),
#'   \item Ends with a stop codon (`TAA`, `TAG`, or `TGA`),
#'   \item Contains no internal stop codons (i.e., stop codons appear only 
#'         at the end).
#' }
#' Otherwise, it returns `FALSE`.
#'
#' @examples
#' IsValidCoding("ATGGACGACTAA")  # TRUE
#' IsValidCoding("ATGTAGGACTAA")  # FALSE due to internal stop codon
#'
#' @export
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

#' Convert a genomic interval from 0-based to 1-based coordinates
#'
#' This function takes a genomic interval in 0-based format 
#' (as used in BED files) and converts it to 1-based format.
#'
#' @param interval A character string representing a genomic interval 
#'                 in the format `"chr:start-end"`, where `start` is 
#'                 0-based and `end` is exclusive.
#'
#' @return A character string representing the interval in 1-based 
#'         inclusive format, or `NA` if the input is malformed 
#'         or doesn't match the expected pattern.
#'
#' @details The function assumes the input interval matches the format 
#'          `"chr:start-end"`, where `start` and `end` are integers. 
#'          It increments the start coordinate by 1 and keeps the end unchanged.
#'
#' @examples
#' ZeroToOneBased("chrIV:0-100")      # "chrIV:1-100"
#' ZeroToOneBased("chrX:499-1000")    # "chrX:500-1000"
#' ZeroToOneBased("chrIII:abc-def")   # NA
#'
#' @export
ZeroToOneBased <- function(interval) {
  match <- regexec("^([^:]+):(\\d+)-(\\d+)$", interval)
  parts <- regmatches(interval, match)[[1]]
  ### malformed input
  if (length(parts) != 4) return(NA)
  chr <- parts[2]
  start <- as.integer(parts[3]) + 1
  end <- as.integer(parts[4])
  
  return(sprintf("%s:%d-%d", chr, start, end))
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
                 "N_feats_sys", "N_feats_rid", "Sblock_type",
                 "N_total", "N_invalid")
### output folders
dirSeqGenes <- file.path(dirBase, "seqs", strSblock)
unlink(dirSeqGenes, recursive = T)
dir.create(dirSeqGenes, recursive = T)
dirSts <- file.path(dirBase, "sts", strSblock)
pathPanClassValidSeq <- file.path(dirSts, "sts-by-haplos-class-valid.txt")

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

### dev indR <- 1
for (indR in seq(1:nrow(dtPanFeatsGns))) {
  lsRow <- dtPanFeatsGns[indR]
  firstFeat <- sub(",.*", "", lsRow$Features_id)
  pathSblockFa <- file.path(dirSeqGenes,
                            paste0(lsRow$Sblock_type, "-", firstFeat, ".fa"))
  bonCols <- setdiff(names(lsRow), skippedCols)
  
  ### dev oneCol <- bonCols[5]
  for (oneCol in bonCols) {
    if (!grepl("#", oneCol)) next
    allIntervals <- unlist(strsplit(lsRow[[oneCol]], ";"))
    if (length(allIntervals) == 0 || all(is.na(allIntervals))) next
    
    ### dev oneInterval <- allIntervals[1]
    for (oneInterval in allIntervals) {
      if (is.na(oneInterval) || oneInterval == "") next
      oneInterval <- ZeroToOneBased(oneInterval)
      lsRow$N_total <- lsRow$N_total + 1
      strGenome <- paste0(sub("#", "-", oneCol), "-genome.fa.gz$")
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
      strSeq <- GetSequence(pathGenome, oneInterval)
      hdFasta <- paste0(oneCol, "#", oneInterval)
      
      if (!is.na(strSeq) && IsValidCoding(strSeq)) {
        cat(sprintf(">%s\n%s\n", hdFasta, strSeq),
            file = pathSblockFa, append = T)
      } else {
        strSeqRevComp <- as.character(reverseComplement(DNAString(strSeq)))
        hdFastaRevComp <- paste0(oneCol, "#RC#", oneInterval)
        if (!is.na(strSeqRevComp) && IsValidCoding(strSeqRevComp)) {
          cat(sprintf(">%s\n%s\n", hdFastaRevComp, strSeqRevComp),
              file = pathSblockFa, append = T)
        } else {
          invInterval <- paste0("inv#", oneInterval)
          lsRow[[oneCol]] <- gsub(oneInterval, invInterval,
                                  lsRow[[oneCol]], fixed = T)
          lsRow$N_invalid <- lsRow$N_invalid + 1
        }
      }
    }
  }
  dtPanFeatsGns[indR,
                c("N_total", "N_invalid") := list(
                  lsRow$N_total,
                  lsRow$N_invalid)]
}

### save dtPanFeatsGns with class annotation (e.g. core) and invalid counts
fwrite(file = pathPanClassValidSeq, x = dtPanFeatsGns, sep = "\t",
       row.names = F, col.names = T)
