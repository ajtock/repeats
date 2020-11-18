#!/applications/R/R-3.5.0/bin/Rscript

# author: Andy Tock
# contact: ajt200@cam.ac.uk
# date: 17.11.2020

# Group CEN180 sequences (identified by Piotr and Ian) in RaGOO v2.0 (assembled by Matt)
# into quantiles according to decreasing:
# 1. weighted SNV values relative to a CEN180 consensus for each chromosome (calculated by Piotr)
# 2. number of CEN180 sequences in a tandem repeat array of >= 2 near-contiguous CEN180 sequences
#    (<= 10 bp apart) on the same strand
# 3. coverage for various data sets (e.g., mean CENH3 ChIP-seq log2(ChIP/input) values)

# Usage:
# /applications/R/R-3.5.0/bin/Rscript group_CEN180_into_quantiles.R 'Chr1,Chr2,Chr3,Chr4,Chr5' 4

#chrName <- unlist(strsplit("Chr1,Chr2,Chr3,Chr4,Chr5",
#                           split = ","))
#quantiles <- 4

args <- commandArgs(trailingOnly = T)
chrName <- unlist(strsplit(args[1],
                           split = ","))
quantiles <- as.integer(args[2])

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
# of >= 2 near-contiguous CEN180 sequences (<= 10 bp apart) on the same strand
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

# Calculate mean log2(CENH3 ChIP/input) in CEN180 sequences 

# Load feature matrices for CENH3 and input, calculate log2(ChIP/control) coverage
ChIPNames <- c(
               "WT_CENH3_Rep1_ChIP_SRR4430537"
              )
ChIPNamesDir <- c(
                  "CENH3_seedlings_Maheshwari_Comai_2017_GenomeRes"
                 )
log2ChIPNamesPlot <- c(
                       "CENH3"
                      )
ChIPDirs <- sapply(seq_along(ChIPNamesDir), function(x) {
  paste0("/home/ajt200/analysis/", ChIPNamesDir[x],
         "/snakemake_ChIPseq_RaGOO_v2.0/mapped/CEN180profiles/matrices/")
})

controlNames <- c(
                  "WT_REC8_Myc_Rep1_input"
                 )
controlNamesDir <- c(
                     "REC8_pooled"
                    )
controlNamesPlot <- c(
                      "Input"
                     )
controlDirs <- sapply(seq_along(controlNamesDir), function(x) {
  paste0("/home/ajt200/analysis/", controlNamesDir[x],
         "/snakemake_ChIPseq_RaGOO_v2.0/mapped/CEN180profiles/matrices/")
})

regionBodyLength <- 180
upstream <- 1000
downstream <- 1000
flankName <- "1kb"
flankNamePlot <- "1 kb"
binSize <- 10
binName <- "10bp"

#WT_CENH3_Rep1_ChIP_SRR4430537_MappedOn_Athaliana_ONT_RaGOO_v2.0_lowXM_both_sort_norm_CEN180_in_Chr1_matrix_bin10bp_flank1kb.tab
#WT_CENH3_Rep1_ChIP_SRR4430537_MappedOn_Athaliana_ONT_RaGOO_v2.0_lowXM_both_sort_norm_CEN180_in_Chr1_ranLoc_matrix_bin10bp_flank1kb.tab
## ChIP
# feature
ChIP_featureMats <- mclapply(seq_along(ChIPNames), function(x) {
  lapply(seq_along(chrName), function(y) {
    as.matrix(read.table(paste0(ChIPDirs[x],
                                ChIPNames[x],
                                "_MappedOn_Athaliana_ONT_RaGOO_v2.0_lowXM_both_sort_norm_CEN180_in_",
                                chrName[y], "_matrix_bin", binName, "_flank", flankName, ".tab"),
                         header = F, skip = 3))
  })
}, mc.cores = length(ChIPNames))
# If features from all 5 chromosomes are to be analysed,
# concatenate the 5 corresponding feature coverage matrices
ChIP_featureMats <- mclapply(seq_along(ChIP_featureMats), function(x) {
  if(length(chrName) == 5) {
    do.call(rbind, ChIP_featureMats[[x]])
  } else {
    ChIP_featureMats[[x]][[1]]
  }
}, mc.cores = length(ChIP_featureMats))

## control
# feature
control_featureMats <- mclapply(seq_along(controlNames), function(x) {
  lapply(seq_along(chrName), function(y) {
    as.matrix(read.table(paste0(controlDirs[x],
                                controlNames[x],
                                "_MappedOn_Athaliana_ONT_RaGOO_v2.0_lowXM_both_sort_norm_CEN180_in_",
                                chrName[y], "_matrix_bin", binName, "_flank", flankName, ".tab"),
                         header = F, skip = 3))
  })
}, mc.cores = length(controlNames))
# If features from all 5 chromosomes are to be analysed,
# concatenate the 5 corresponding feature coverage matrices
control_featureMats <- mclapply(seq_along(control_featureMats), function(x) {
  if(length(chrName) == 5) {
    do.call(rbind, control_featureMats[[x]])
  } else {
    control_featureMats[[x]][[1]]
  }
}, mc.cores = length(control_featureMats))

# Calculate log2(ChIP/control) coverage values
log2ChIP_featureMats <- mclapply(seq_along(ChIP_featureMats), function(x) {
  log2((ChIP_featureMats[[x]]+1)/(control_featureMats[[1]]+1))
}, mc.cores = length(ChIP_featureMats))

# Calculate mean log2(ChIP/control) coverage values for each CEN180 sequence
log2ChIP_featureMats_bodies <- lapply(seq_along(log2ChIP_featureMats), function(x) {
  log2ChIP_featureMats[[x]][,((upstream/binSize)+1):((upstream+regionBodyLength)/binSize)]
})
log2ChIP_featureMats_bodiesRowMeans <- lapply(seq_along(log2ChIP_featureMats_bodies), function(x) {
  rowMeans(log2ChIP_featureMats_bodies[[x]], na.rm = T)
})

# Add mean coverage values to CEN180 dataframe
CEN180 <- data.frame(CEN180,
                     CENH3_in_bodies = log2ChIP_featureMats_bodiesRowMeans[[1]])

# Define set of ordering factors to be used for grouping genes into quantiles
orderingFactor <- colnames(CEN180)[c(5, 9, 10)]
outDir <- paste0("quantiles_by_", orderingFactor, "/",
                 paste0(chrName, collapse = "_"), "/")
plotDir <- paste0(outDir, "/plots/")
sapply(seq_along(outDir), function(w) {
 system(paste0("[ -d ", outDir[w], " ] || mkdir -p ", outDir[w]))
})
sapply(seq_along(plotDir), function(w) {
 system(paste0("[ -d ", plotDir[w], " ] || mkdir -p ", plotDir[w]))
})

# Group features into quantiles according to decreasing orderingFactor
mclapply(seq_along(orderingFactor), function(w) {
  # Note that CEN180_DF could be defined on
  # line above or below mclapply(), with the same results
  CEN180_DF <- data.frame(CEN180,
                          quantile = as.character(""))
  print(orderingFactor[w])
  # Assign 0s to NA values only for coverage data
  if(grepl("_in_", orderingFactor[w])) {
    CEN180_DF[,which(colnames(CEN180_DF) == orderingFactor[w])][
      which(is.na(CEN180_DF[,which(colnames(CEN180_DF) == orderingFactor[w])]))] <- 0
  }
  quantilesStats <- data.frame()
  for(k in 1:quantiles) {
    # First quantile should span 1 to greater than, e.g., 0.75 proportions of features
    if(k < quantiles) {
      CEN180_DF[ !is.na(CEN180_DF[,which(colnames(CEN180_DF) == orderingFactor[w])]) &
                 percent_rank(CEN180_DF[,which(colnames(CEN180_DF) == orderingFactor[w])]) <= 1-((k-1)/quantiles) &
                 percent_rank(CEN180_DF[,which(colnames(CEN180_DF) == orderingFactor[w])]) >  1-(k/quantiles), ]$quantile <- paste0("Quantile ", k)
    } else {
    # Final quantile should span 0 to, e.g., 0.25 proportions of features
      CEN180_DF[ !is.na(CEN180_DF[,which(colnames(CEN180_DF) == orderingFactor[w])]) &
                 percent_rank(CEN180_DF[,which(colnames(CEN180_DF) == orderingFactor[w])]) <= 1-((k-1)/quantiles) &
                 percent_rank(CEN180_DF[,which(colnames(CEN180_DF) == orderingFactor[w])]) >= 1-(k/quantiles), ]$quantile <- paste0("Quantile ", k)
    }
    write.table(CEN180_DF[CEN180_DF$quantile == paste0("Quantile ", k),],
                file = paste0(outDir[w],
                              "quantile", k, "_of_", quantiles,
                              "_by_", orderingFactor[w],
                              "_of_CEN180_in_RaGOO_v2.0_",
                              paste0(chrName, collapse = "_"), ".tsv"),
                quote = FALSE, sep = "\t", row.names = FALSE)
    stats <- data.frame(quantile = as.integer(k),
                        n = as.integer(dim(CEN180_DF[CEN180_DF$quantile == paste0("Quantile ", k),])[1]),
                        mean_width = as.integer(round(mean(
                          CEN180_DF[CEN180_DF$quantile == paste0("Quantile ", k),]$end -
                          CEN180_DF[CEN180_DF$quantile == paste0("Quantile ", k),]$start0based, na.rm = T))),
                        total_width = as.integer(sum(
                          CEN180_DF[CEN180_DF$quantile == paste0("Quantile ", k),]$end -
                          CEN180_DF[CEN180_DF$quantile == paste0("Quantile ", k),]$start0based, na.rm = T)),
                        mean_orderingFactor = as.numeric(mean(CEN180_DF[CEN180_DF$quantile == paste0("Quantile ", k),][,which(colnames(CEN180_DF) == orderingFactor[w])], na.rm = T)))
    quantilesStats <- rbind(quantilesStats, stats)
  }
  write.table(quantilesStats,
              file = paste0(outDir[w],
                            "summary_", quantiles, "quantiles",
                            "_by_", orderingFactor[w],
                            "_of_CEN180_in_RaGOO_v2.0_",
                            paste0(chrName, collapse = "_"), ".tsv"),
              quote = FALSE, sep = "\t", row.names = FALSE)
  write.table(CEN180_DF,
              file = paste0(outDir[w],
                            "features_", quantiles, "quantiles",
                            "_by_", orderingFactor[w],
                            "_of_CEN180_in_RaGOO_v2.0_",
                            paste0(chrName, collapse = "_"), ".tsv"),
              quote = FALSE, sep = "\t", row.names = FALSE)
}, mc.cores = length(orderingFactor), mc.preschedule = F)
