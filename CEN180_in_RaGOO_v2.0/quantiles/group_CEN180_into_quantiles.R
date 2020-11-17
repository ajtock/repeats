#!/applications/R/R-3.5.0/bin/Rscript

# author: Andy Tock
# contact: ajt200@cam.ac.uk
# date: 17.11.2020

# Group CEN180 sequences (identified by Piotr and Ian) in RaGOO v2.0 (assembled by Matt)
# into quantiles according to decreasing:
# 1. weighted SNV values relative to a CEN180 consensus for each chromosome (calculated by Piotr)
# 2. number of CEN180 sequences in a tandem repeat array of >= 2 near-contiguous CEN180 sequences
#    (<= 10 bp apart), such that Quantile 4 contains only singletons (?)
# 3. coverage for various data sets (e.g., mean CENH3 ChIP-seq log2(ChIP/input) values)

# Usage:
# /applications/R/R-3.5.0/bin/Rscript group_CEN180_into_quantiles.R 'Chr1,Chr2,Chr3,Chr4,Chr5'

chrName <- unlist(strsplit("Chr1,Chr2,Chr3,Chr4,Chr5",
                           split = ","))

args <- commandArgs(trailingOnly = T)
chrName <- unlist(strsplit(args[1],
                           split = ","))

options(stringsAsFactors = F)
library(parallel)
library(dplyr)

# Load table of CEN180 coordinates in BED format
CEN180 <- read.table(paste0("CEN180_in_RaGOO_v2.0_",
                            paste0(chrName, collapse = "_"), ".bed"),
                     header = F)
colnames(CEN180) <- c("chr", "start0based", "end", "featureID", "wSNV", "strand")

# Get CEN180 coordinates within chrName
CEN180 <- CEN180[which(CEN180$chr %in% chrName),]

# Determine if each CEN180 sequence is part of a tandem repeat array
# of >= 2 near-contiguous CEN180 sequences (<= 10 bp apart)
# Assign an array ID number to each CEN180 within a tandem repeat array
# and count the number of tandem repeats within the array
CEN180$tandem_repeat <- as.logical("")
CEN180$array <- as.integer("")
CEN180$array_size <- as.integer("")
CEN180_TR <- NULL
for(i in seq_along(chrName)) {
  print(chrName[i])
  CEN180_chr <- CEN180[CEN180$chr == chrName[i],]
  # First CEN180 in a chromosome
  if ( CEN180_chr[2,]$start0based - CEN180_chr[1,]$end < 9 &
       CEN180_chr[2,]$strand == CEN180_chr[1,]$strand ) {
    CEN180_chr[1,]$tandem_repeat <- TRUE
    CEN180_chr[1,]$array <- as.integer(1)
  } else {
    CEN180_chr[1,]$tandem_repeat <- FALSE
  }
  # All other CEN180 sequences in a chromosome, except the last
  for(j in 2:(dim(CEN180_chr)[1]-1)) {
    if ( CEN180_chr[j,]$start0based - CEN180_chr[j-1,]$end < 9 &
         CEN180_chr[j,]$strand == CEN180_chr[j-1,]$strand ) {
      CEN180_chr[j,]$tandem_repeat <- TRUE
      CEN180_chr[j,]$array <- as.integer(CEN180_chr[j-1,]$array)
    } else if ( CEN180_chr[j+1,]$start0based - CEN180_chr[j,]$end < 9 &
                CEN180_chr[j+1,]$strand == CEN180_chr[j,]$strand ) {
      CEN180_chr[j,]$tandem_repeat <- TRUE
      if ( length(CEN180_chr[!is.na(CEN180_chr$array),]$array) > 0 ) {
        CEN180_chr[j,]$array <- as.integer(CEN180_chr[!is.na(CEN180_chr$array),]$array[
                                  length(CEN180_chr[!is.na(CEN180_chr$array),]$array)] + 1)
      } else {
        CEN180_chr[j,]$array <- as.integer(1)
      }
    } else {
      CEN180_chr[j,]$tandem_repeat <- FALSE
    }
  }
  # Last CEN180 in a chromosome
  if ( CEN180_chr[dim(CEN180_chr)[1],]$start0based - CEN180_chr[dim(CEN180_chr)[1]-1,]$end < 9 &
       CEN180_chr[dim(CEN180_chr)[1],]$strand == CEN180_chr[dim(CEN180_chr)[1]-1,]$strand ) {
    CEN180_chr[dim(CEN180_chr)[1],]$tandem_repeat <- TRUE
    CEN180_chr[dim(CEN180_chr)[1],]$array <- as.integer(CEN180_chr[dim(CEN180_chr)[1]-1,]$array)
  } else {
    CEN180_chr[dim(CEN180_chr)[1],]$tandem_repeat <- FALSE
  }
  for(k in 1:dim(CEN180_chr)[1]) {
    if ( CEN180_chr[k,]$tandem_repeat == TRUE ) {
      CEN180_chr[k,]$array_size <- as.integer(dim(CEN180_chr[
                                     which(CEN180_chr$array == CEN180_chr[k,]$array),])[1])
    } else {
      CEN180_chr[k,]$array_size <- as.integer(1)
    }
  } 
  CEN180_TR <- rbind(CEN180_TR, CEN180_chr)
}
CEN180 <- CEN180_TR

# Define set of ordering factors to be used for grouping genes into 4 quantiles
orderingFactor <- colnames(featuresNLR_pop_list[[1]])[c(7:30, 79:103, 109)]
outDir <- paste0("quantiles_by_", orderingFactor, "/")
outDir_list <- lapply(seq_along(outDir), function(w) {
  sapply(seq_along(pop_name), function(x) {
    paste0(outDir[w], pop_name[x], "/")
  })
})
plotDir_list <- lapply(seq_along(outDir), function(w) {
  sapply(seq_along(pop_name), function(x) {
    paste0(outDir_list[[w]][x], "plots/")
  })
})
#sapply(seq_along(outDir), function(w) {
# system(paste0("[ -d ", outDir[w], " ] || mkdir ", outDir[w]))
#})
#sapply(seq_along(outDir), function(w) {
#  mclapply(seq_along(pop_name), function(x) {
#    system(paste0("[ -d ", outDir_list[[w]][x], " ] || mkdir ", outDir_list[[w]][x]))
#  }, mc.cores = length(pop_name), mc.preschedule = F)
#})
#sapply(seq_along(outDir), function(w) {
#  mclapply(seq_along(pop_name), function(x) {
#    system(paste0("[ -d ", plotDir_list[[w]][x], " ] || mkdir ", plotDir_list[[w]][x]))
#  }, mc.cores = length(pop_name), mc.preschedule = F)
#})

# For each population, divide features into quantiles based on decreasing orderingFactor
for(x in 1:length(featuresNLR_pop_list)) {
  print(pop_name[x])
  NLR_quants <- data.frame(featuresNLR_pop_list[[x]],
                           NLR_quantile = as.character(""),
                           stringsAsFactors = F)
  mclapply(seq_along(orderingFactor), function(w) {
    print(orderingFactor[w])
    # Assign 0s to NA values only for coverage data
    if(grepl("_in_", orderingFactor[w])) {
      NLR_quants[,which(colnames(NLR_quants) == orderingFactor[w])][which(is.na(NLR_quants[,which(colnames(NLR_quants) == orderingFactor[w])]))] <- 0
    }
    if(grepl("all", orderingFactor[w])) {
    # 2 quantiles for population genetics stats quantiles
      quantiles <- 2
    } else {
    # 4 quantiles for coverage and cM/Mb quantiles
     quantiles <- 4
    }
    quantilesStats <- data.frame()
    for(k in 1:quantiles) {
      # First quantile should span 1 to greater than, e.g., 0.75 proportions of features
      if(k < quantiles) {
        NLR_quants[ !is.na(NLR_quants[,which(colnames(NLR_quants) == orderingFactor[w])]) &
                    percent_rank(NLR_quants[,which(colnames(NLR_quants) == orderingFactor[w])]) <= 1-((k-1)/quantiles) &
                    percent_rank(NLR_quants[,which(colnames(NLR_quants) == orderingFactor[w])]) >  1-(k/quantiles), ]$NLR_quantile <- paste0("Quantile ", k)
      } else {
      # Final quantile should span 0 to, e.g., 0.25 proportions of features
        NLR_quants[ !is.na(NLR_quants[,which(colnames(NLR_quants) == orderingFactor[w])]) &
                    percent_rank(NLR_quants[,which(colnames(NLR_quants) == orderingFactor[w])]) <= 1-((k-1)/quantiles) &
                    percent_rank(NLR_quants[,which(colnames(NLR_quants) == orderingFactor[w])]) >= 1-(k/quantiles), ]$NLR_quantile <- paste0("Quantile ", k)
      }
      write.table(NLR_quants[NLR_quants$NLR_quantile == paste0("Quantile ", k),],
                  file = paste0(outDir_list[[w]][x],
                                "quantile", k, "_of_", quantiles,
                                "_by_", orderingFactor[w],
                                "_of_NLR_genes_in_",
                                paste0(substring(chrName, first = 10, last = 16),
                                       collapse = "_"), "_",
                                substring(chrName[1][1], first = 18), "_", pop_name[x], ".txt"),
                  quote = FALSE, sep = "\t", row.names = FALSE)
      stats <- data.frame(quantile = as.integer(k),
                          n = as.integer(dim(NLR_quants[NLR_quants$NLR_quantile == paste0("Quantile ", k),])[1]),
                          mean_width = as.integer(round(mean(NLR_quants[NLR_quants$NLR_quantile == paste0("Quantile ", k),]$width, na.rm = T))),
                          total_width = as.integer(sum(NLR_quants[NLR_quants$NLR_quantile == paste0("Quantile ", k),]$width, na.rm = T)),
                          mean_orderingFactor = as.numeric(mean(NLR_quants[NLR_quants$NLR_quantile == paste0("Quantile ", k),][,which(colnames(NLR_quants) == orderingFactor[w])], na.rm = T)))
      quantilesStats <- rbind(quantilesStats, stats)
    }
    write.table(quantilesStats,
                file = paste0(outDir_list[[w]][x],
                              "summary_", quantiles, "quantiles",
                              "_by_", orderingFactor[w],
                              "_of_NLR_genes_in_",
                              paste0(substring(chrName, first = 10, last = 16),
                                     collapse = "_"), "_",
                              substring(chrName[1][1], first = 18), "_", pop_name[x], ".txt"),
                quote = FALSE, sep = "\t", row.names = FALSE)
    write.table(NLR_quants,
                file = paste0(outDir_list[[w]][x],
                              "features_", quantiles, "quantiles",
                              "_by_", orderingFactor[w],
                              "_of_NLR_genes_in_",
                              paste0(substring(chrName, first = 10, last = 16),
                                     collapse = "_"), "_",
                              substring(chrName[1][1], first = 18), "_", pop_name[x], ".txt"),
                quote = FALSE, sep = "\t", row.names = FALSE)
  }, mc.cores = length(orderingFactor), mc.preschedule = F)
}
