#!/applications/R/R-4.0.0/bin/Rscript

# author: Andy Tock
# contact: ajt200@cam.ac.uk
# date: 12.01.2021

# Group CEN180 sequences (identified by Piotr and Ian) in T2T_Col (assembled partly by Matt)
# into quantiles according to decreasing:
# 1. weighted SNV values relative to a genome-wide CEN180 consensus (calculated by Piotr)
# 2. number of CEN180 sequences in a tandem repeat array of >= 2 near-contiguous CEN180 sequences
#    (<= 10 bp apart) on the same strand
# 3. repeat activity (calculated by Piotr; column titled "HORlengthsSum" by Piotr; i.e., the sum of repeat units that make up the HORs of which a given repeat unit is a member)
# 4. coverage for various data sets (e.g., mean CENH3 ChIP-seq log2(ChIP/input) values)
# 5. mappability given various k-mers

# Usage:
# /applications/R/R-4.0.0/bin/Rscript group_CEN180_into_quantiles.R 'Chr1,Chr2,Chr3,Chr4,Chr5' 4 180 1000 1kb '1 kb' 10 10bp

#chrName <- unlist(strsplit("Chr1,Chr2,Chr3,Chr4,Chr5",
#                           split = ","))
#quantiles <- 4
#regionBodyLength <- 180
#upstream <- 1000
#downstream <- 1000
#flankName <- "1kb"
#flankNamePlot <- "1 kb"
#binSize <- 10
#binName <- "10bp"

args <- commandArgs(trailingOnly = T)
chrName <- unlist(strsplit(args[1],
                           split = ","))
quantiles <- as.integer(args[2])
regionBodyLength <- as.integer(args[3])
upstream <- as.integer(args[4])
downstream <- as.integer(args[4])
flankName <- args[5]
flankNamePlot <- args[6]
binSize <- as.integer(args[7])
binName <- args[8]

options(stringsAsFactors = F)
library(parallel)
library(dplyr)
library(Hmisc) # includes rcorr() function which computes significance levels for Pearson and Spearman correlations
library(reshape)
library(ggplot2)
library(ggcorrplot)

# Load table of CEN180 coordinates in BED format
CEN180 <- read.table(paste0("/home/ajt200/analysis/repeats/CEN180_in_T2T_Col/CEN180_in_T2T_Col_",
                            paste0(chrName, collapse = "_"), ".bed"),
                     header = F)
# Convert 0-based start coordinates (BED)
# into 1-based start coordinates (for output as TSV below)
CEN180[,2] <- CEN180[,2]+1
colnames(CEN180) <- c("chr", "start", "end", "featureID", "wSNV", "strand", "HORlengthsSum", "HORcount")
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
  if ( CEN180_chr[2,]$start - CEN180_chr[1,]$end <= 10 &
       CEN180_chr[2,]$strand == CEN180_chr[1,]$strand ) {
    CEN180_chr[1,]$tandem_repeat <- TRUE
    CEN180_chr[1,]$array <- as.integer(1)
  } else {
    CEN180_chr[1,]$tandem_repeat <- FALSE
  }
  # All other CEN180 sequences in a chromosome, except the last
  for(j in 2:(dim(CEN180_chr)[1]-1)) {
    if ( CEN180_chr[j,]$start - CEN180_chr[j-1,]$end <= 10 &
         CEN180_chr[j,]$strand == CEN180_chr[j-1,]$strand ) {
      CEN180_chr[j,]$tandem_repeat <- TRUE
      CEN180_chr[j,]$array <- as.integer(CEN180_chr[j-1,]$array)
    } else if ( CEN180_chr[j+1,]$start - CEN180_chr[j,]$end <= 10 &
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
  if ( CEN180_chr[dim(CEN180_chr)[1],]$start - CEN180_chr[dim(CEN180_chr)[1]-1,]$end <= 10 &
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

# Load table of coordinates for random loci in BED format
ranLoc <- read.table(paste0("/home/ajt200/analysis/repeats/CEN180_in_T2T_Col/CEN180_in_T2T_Col_",
                            paste0(chrName, collapse = "_"), "_randomLoci.bed"),
                     header = F)
# Convert 0-based start coordinates (BED)
# into 1-based start coordinates (for output as TSV below)
ranLoc[,2] <- ranLoc[,2]+1
colnames(ranLoc) <- c("chr", "start", "end", "featureID", "wSNV", "strand")
ranLoc <- data.frame(ranLoc,
                     HORlengthsSum = NA,
                     HORcount = NA,
                     tandem_repeat = NA,
                     array = NA,
                     array_size = NA)
# Get ranLoc coordinates within chrName
ranLoc <- ranLoc[which(ranLoc$chr %in% chrName),]

# Calculate mean coverage and mappability in CEN180 sequences and random loci

# Load feature matrices for CENH3 and input, calculate log2(ChIP/control) coverage
ChIPNames <- c(
               "WT_CENH3_Rep1_ChIP_SRR4430537",
               "WT_H3K9me2_Rep1_ChIP",
               "WT_H3K27me1_Rep1_ChIP",
               "WT_H3K4me1_Rep1_ChIP",
               "WT_H3K4me2_Rep1_ChIP",
               "WT_H3K4me3_ChIP14",
               "H2AW6_ChIP_SRR5298545",
               "H2AW7_ChIP_SRR5298546",
               "WT_MNase_Rep1",
               "WT_REC8_HA_Rep2_ChIP",
               "WT_ASY1_Rep1_ChIP",
               "WT_MTOPVIB_HA_Rep1_ChIP",
               "WT_MTOPVIB_HA_Rep2_ChIP",
               "WT_DMC1_V5_Rep1_ChIP",
               "WT_DMC1_V5_Rep2_ChIP",
               "WT_SPO11oligos_Rep1"
              )
ChIPNamesDir <- c(
                  "CENH3_seedlings_Maheshwari_Comai_2017_GenomeRes/snakemake_ChIPseq_T2T_Col",
                  "170101_Chris_H3K9me2_ChIP/WT/snakemake_ChIPseq_T2T_Col",
                  rep("170920_Chris_histone_ChIP/snakemake_ChIPseq_T2T_Col", 3),
                  "160601_Kyuha_H3K4me3_ChIP/WT/snakemake_ChIPseq_T2T_Col",
                  rep("HTA6_HTA7_leaf_Lorkovic_Berger_2017_CurrBiol/snakemake_ChIPseq_T2T_Col", 2),
                  "150701_Kyuha_MNase/WT/snakemake_ChIPseq_T2T_Col",
                  "REC8_pooled/snakemake_ChIPseq_T2T_Col",
                  "20190722_cal66_Athaliana_ChIPseq_ASY1/fastq_pooled/snakemake_ChIPseq_T2T_Col",
                  rep("20190819_dh580_Athaliana_ChIPseq_MTOPVIB/fastq_pooled/snakemake_ChIPseq_T2T_Col", 2),
                  rep("20190917_dh580_Athaliana_ChIPseq_DMC1/fastq_pooled/snakemake_ChIPseq_T2T_Col", 2),
                  "160518_Kyuha_SPO11oligos/WT/snakemake_SPO11oligos_T2T_Col"
                 )
log2ChIPNamesPlot <- c(
                       "CENH3",
                       "H3K9me2",
                       "H3K27me1",
                       "H3K4me1",
                       "H3K4me2",
                       "H3K4me3",
                       "H2A.W.6",
                       "H2A.W.7",
                       "MNase",
                       "REC8-HA",
                       "ASY1",
                       "MTOPVIB-HA Rep1",
                       "MTOPVIB-HA Rep2",
                       "DMC1-V5 Rep1",
                       "DMC1-V5 Rep2",
                       "SPO11-1-oligos"
                      )
ChIPNamesPlot <- log2ChIPNamesPlot
ChIPDirs <- sapply(seq_along(ChIPNamesDir), function(x) {
  paste0("/home/ajt200/analysis/", ChIPNamesDir[x],
         "/mapped/CEN180profiles/matrices/")
})

controlNames <- c(
                  "WT_REC8_Myc_Rep1_input",
                  "H2AW_input_SRR5298544",
                  "WT_gDNA_Rep1",
                  "WT_gDNA_Rep1_R1",
                  "map_K40_E2",
                  "map_K45_E2",
                  "map_K50_E2",
                  "map_K150_E4",
                  "map_K200_E4",
                  "map_K300_E4"
                 )
controlNamesDir <- c(
                     "REC8_pooled/snakemake_ChIPseq_T2T_Col",
                     "HTA6_HTA7_leaf_Lorkovic_Berger_2017_CurrBiol/snakemake_ChIPseq_T2T_Col",
                     "150701_Natasha_gDNA/WT/snakemake_ChIPseq_T2T_Col",
                     "150701_Natasha_gDNA/WT/R1/snakemake_SPO11oligos_T2T_Col",
                     rep("nanopore/T2T_Col/genmap_mappability", 6)
                    )
controlNamesPlot <- c(
                      "PE input (sonic)",
                      "PE input (MNase)",
                      "PE gDNA",
                      "SE gDNA",
                      "k=40 e=2 mappability",
                      "k=45 e=2 mappability",
                      "k=50 e=2 mappability",
                      "k=150 e=4 mappability",
                      "k=200 e=4 mappability",
                      "k=300 e=4 mappability"
                     )
controlDirs <- sapply(seq_along(controlNamesDir), function(x) {
  paste0("/home/ajt200/analysis/", controlNamesDir[x],
         "/mapped/CEN180profiles/matrices/")
})

## DNAmeth
DNAmethNames <- c(
                  "WT_nanopolishDNAmeth_95_10kb",
                  rep("Col0_BSseq_Rep1", 3),
                  rep("WT_BSseq_Rep2_2013", 3)
                 )
DNAmethNamesDir <- c(
                     "nanopore/T2T_Col/nanopolish_DNAmeth",
                     rep("BSseq_seedling_Yang_Zhu_2016_CellRes/snakemake_BSseq_T2T_Col/coverage", 3),
                     rep("BSseq_leaf_Stroud_Jacobsen_2013_Cell_2014_NSMB/snakemake_BSseq_T2T_Col/coverage", 3)
                    )
DNAmethContexts <- c(
                     "CpG",
                     "CpG",
                     "CHG",
                     "CHH",
                     "CpG",
                     "CHG",
                     "CHH"
                    )
DNAmethNamesPlot <- c(
                      "mCG (Nanopolish)",
                      "mCG (PE BS-seq)",
                      "mCHG (PE BS-seq)",
                      "mCHH (PE BS-seq)",
                      "mCG (SE BS-seq)",
                      "mCHG (SE BS-seq)",
                      "mCHH (SE BS-seq)"
                     )
DNAmethDirs <- sapply(seq_along(DNAmethNamesDir), function(x) {
  paste0("/home/ajt200/analysis/",
         DNAmethNamesDir[x],
         "/CEN180profiles/matrices/")
})

## DNAmeth
# feature
DNAmeth_featureMats <- mclapply(seq_along(DNAmethContexts), function(x) {
  lapply(seq_along(chrName), function(y) {
    as.matrix(read.table(paste0(DNAmethDirs[x],
                                DNAmethNames[x],
                                "_MappedOn_T2T_Col_", DNAmethContexts[x], "_CEN180_in_",
                                chrName[y], "_matrix_bin", binName, "_flank", flankName, ".tab"),
                         header = F, skip = 3))
  })
}, mc.cores = length(DNAmethContexts))
# If features from multiple chromosomes are to be analysed,
# concatenate the corresponding feature coverage matrices
DNAmeth_featureMats <- mclapply(seq_along(DNAmeth_featureMats), function(x) {
  if(length(chrName) > 1) {
    do.call(rbind, DNAmeth_featureMats[[x]])
  } else {
    DNAmeth_featureMats[[x]][[1]]
  }
}, mc.cores = length(DNAmeth_featureMats))

## ChIP
# feature
ChIP_featureMats <- mclapply(seq_along(ChIPNames), function(x) {
  lapply(seq_along(chrName), function(y) {
    as.matrix(read.table(paste0(ChIPDirs[x],
                                ChIPNames[x],
                                "_MappedOn_T2T_Col_lowXM_both_sort_norm_CEN180_in_",
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
    if( grepl("map_K", controlNames[x]) ) {
      as.matrix(read.table(paste0(controlDirs[x],
                                  controlNames[x],
                                  "_MappedOn_T2T_Col_CEN180_in_",
                                  chrName[y], "_matrix_bin", binName, "_flank", flankName, ".tab"),
                           header = F, skip = 3))
    } else {
      as.matrix(read.table(paste0(controlDirs[x],
                                  controlNames[x],
                                  "_MappedOn_T2T_Col_lowXM_both_sort_norm_CEN180_in_",
                                  chrName[y], "_matrix_bin", binName, "_flank", flankName, ".tab"),
                           header = F, skip = 3))
      }
  })
}, mc.cores = length(controlNames))
# If features from multiple chromosomes are to be analysed,
# concatenate the corresponding feature coverage matrices
control_featureMats <- mclapply(seq_along(control_featureMats), function(x) {
  if(length(chrName) > 1) {
    do.call(rbind, control_featureMats[[x]])
  } else {
    control_featureMats[[x]][[1]]
  }
}, mc.cores = length(control_featureMats))

# Conditionally calculate log2(ChIP/control)
# for each matrix depending on library
# feature
log2ChIP_featureMats <- mclapply(seq_along(ChIP_featureMats), function(x) {
  if ( grepl("SRR529854", ChIPNames[x]) ) {
    print(paste0(ChIPNames[x], " library; using ", controlNames[2], " for log2((MNase+1)/(gDNA+1)) calculation"))
    log2((ChIP_featureMats[[x]]+1)/(control_featureMats[[2]]+1))
  } else if ( grepl("MNase", ChIPNames[x]) ) {
    print(paste0(ChIPNames[x], " library; using ", controlNames[3], " for log2((SPO11-1-oligos+1)/(gDNA+1)) calculation"))
    log2((ChIP_featureMats[[x]]+1)/(control_featureMats[[3]]+1))
  } else if ( grepl("SPO11oligos", ChIPNames[x]) ) {
    print(paste0(ChIPNames[x], " library; using ", controlNames[4], " for log2((SPO11-1-oligos+1)/(gDNA+1)) calculation"))
    log2((ChIP_featureMats[[x]]+1)/(control_featureMats[[4]]+1))
  } else {
    print(paste0(ChIPNames[x], " library; using ", controlNames[1], " for log2((ChIP+1)/(input+1)) calculation"))
    log2((ChIP_featureMats[[x]]+1)/(control_featureMats[[1]]+1))
  }
}, mc.cores = length(ChIP_featureMats))

# Calculate mean DNAmeth-dataset values for each CEN180 sequence
DNAmeth_featureMats_bodies <- lapply(seq_along(DNAmeth_featureMats), function(x) {
  DNAmeth_featureMats[[x]][,((upstream/binSize)+1):((upstream+regionBodyLength)/binSize)]
})
DNAmeth_featureMats_bodiesRowMeans <- lapply(seq_along(DNAmeth_featureMats_bodies), function(x) {
  rowMeans(DNAmeth_featureMats_bodies[[x]], na.rm = T)
})

# Calculate mean log2(ChIP/control) coverage values for each CEN180 sequence
log2ChIP_featureMats_bodies <- lapply(seq_along(log2ChIP_featureMats), function(x) {
  log2ChIP_featureMats[[x]][,((upstream/binSize)+1):((upstream+regionBodyLength)/binSize)]
})
log2ChIP_featureMats_bodiesRowMeans <- lapply(seq_along(log2ChIP_featureMats_bodies), function(x) {
  rowMeans(log2ChIP_featureMats_bodies[[x]], na.rm = T)
})

# Calculate mean control-dataset values for each CEN180 sequence
control_featureMats_bodies <- lapply(seq_along(control_featureMats), function(x) {
  control_featureMats[[x]][,((upstream/binSize)+1):((upstream+regionBodyLength)/binSize)]
})
control_featureMats_bodiesRowMeans <- lapply(seq_along(control_featureMats_bodies), function(x) {
  rowMeans(control_featureMats_bodies[[x]], na.rm = T)
})

# Add mean coverage values to CEN180 dataframe
CEN180 <- data.frame(CEN180,
                     CENH3_in_bodies = log2ChIP_featureMats_bodiesRowMeans[[1]],
                     H3K9me2_in_bodies = log2ChIP_featureMats_bodiesRowMeans[[2]],
                     H3K27me1_in_bodies = log2ChIP_featureMats_bodiesRowMeans[[3]],
                     H3K4me1_in_bodies = log2ChIP_featureMats_bodiesRowMeans[[4]],
                     H3K4me2_in_bodies = log2ChIP_featureMats_bodiesRowMeans[[5]],
                     H3K4me3_in_bodies = log2ChIP_featureMats_bodiesRowMeans[[6]],
                     H2AW6_in_bodies = log2ChIP_featureMats_bodiesRowMeans[[7]],
                     H2AW7_in_bodies = log2ChIP_featureMats_bodiesRowMeans[[8]],
                     MNase_in_bodies = log2ChIP_featureMats_bodiesRowMeans[[9]],
                     REC8_HA_Rep2_in_bodies = log2ChIP_featureMats_bodiesRowMeans[[10]],
                     ASY1_Rep1_in_bodies = log2ChIP_featureMats_bodiesRowMeans[[11]],
                     MTOPVIB_HA_Rep1_in_bodies = log2ChIP_featureMats_bodiesRowMeans[[12]],
                     MTOPVIB_HA_Rep2_in_bodies = log2ChIP_featureMats_bodiesRowMeans[[13]],
                     DMC1_V5_Rep1_in_bodies = log2ChIP_featureMats_bodiesRowMeans[[14]],
                     DMC1_V5_Rep2_in_bodies = log2ChIP_featureMats_bodiesRowMeans[[15]],
                     SPO11oligos_in_bodies = log2ChIP_featureMats_bodiesRowMeans[[16]],
                     mCG_Nanopolish_in_bodies = DNAmeth_featureMats_bodiesRowMeans[[1]], 
                     mCG_PE_BSseq_in_bodies = DNAmeth_featureMats_bodiesRowMeans[[2]], 
                     mCHG_PE_BSseq_in_bodies = DNAmeth_featureMats_bodiesRowMeans[[3]], 
                     mCHH_PE_BSseq_in_bodies = DNAmeth_featureMats_bodiesRowMeans[[4]], 
                     mCG_SE_BSseq_in_bodies = DNAmeth_featureMats_bodiesRowMeans[[5]], 
                     mCHG_SE_BSseq_in_bodies = DNAmeth_featureMats_bodiesRowMeans[[6]], 
                     mCHH_SE_BSseq_in_bodies = DNAmeth_featureMats_bodiesRowMeans[[7]], 
                     PE_input_sonic_in_bodies = control_featureMats_bodiesRowMeans[[1]],
                     PE_input_MNase_in_bodies = control_featureMats_bodiesRowMeans[[2]],
                     PE_gDNA_in_bodies = control_featureMats_bodiesRowMeans[[3]],
                     SE_gDNA_in_bodies = control_featureMats_bodiesRowMeans[[4]],
                     map_K40_E2_in_bodies = control_featureMats_bodiesRowMeans[[5]],
                     map_K45_E2_in_bodies = control_featureMats_bodiesRowMeans[[6]],
                     map_K50_E2_in_bodies = control_featureMats_bodiesRowMeans[[7]],
                     map_K150_E4_in_bodies = control_featureMats_bodiesRowMeans[[8]],
                     map_K200_E4_in_bodies = control_featureMats_bodiesRowMeans[[9]],
                     map_K300_E4_in_bodies = control_featureMats_bodiesRowMeans[[10]])

## DNAmeth
# ranLoc
DNAmeth_ranLocMats <- mclapply(seq_along(DNAmethContexts), function(x) {
  lapply(seq_along(chrName), function(y) {
    as.matrix(read.table(paste0(DNAmethDirs[x],
                                DNAmethNames[x],
                                "_MappedOn_T2T_Col_", DNAmethContexts[x], "_CEN180_in_",
                                chrName[y], "_ranLoc_matrix_bin", binName, "_flank", flankName, ".tab"),
                         header = F, skip = 3))
  })
}, mc.cores = length(DNAmethContexts))
# If ranLocs from multiple chromosomes are to be analysed,
# concatenate the corresponding ranLoc coverage matrices
DNAmeth_ranLocMats <- mclapply(seq_along(DNAmeth_ranLocMats), function(x) {
  if(length(chrName) > 1) {
    do.call(rbind, DNAmeth_ranLocMats[[x]])
  } else {
    DNAmeth_ranLocMats[[x]][[1]]
  }
}, mc.cores = length(DNAmeth_ranLocMats))

## ChIP
# ranLoc
ChIP_ranLocMats <- mclapply(seq_along(ChIPNames), function(x) {
  lapply(seq_along(chrName), function(y) {
    as.matrix(read.table(paste0(ChIPDirs[x],
                                ChIPNames[x],
                                "_MappedOn_T2T_Col_lowXM_both_sort_norm_CEN180_in_",
                                chrName[y], "_ranLoc_matrix_bin", binName, "_flank", flankName, ".tab"),
                         header = F, skip = 3))
  })
}, mc.cores = length(ChIPNames))
# If ranLocs from all 5 chromosomes are to be analysed,
# concatenate the 5 corresponding ranLoc coverage matrices
ChIP_ranLocMats <- mclapply(seq_along(ChIP_ranLocMats), function(x) {
  if(length(chrName) == 5) {
    do.call(rbind, ChIP_ranLocMats[[x]])
  } else {
    ChIP_ranLocMats[[x]][[1]]
  }
}, mc.cores = length(ChIP_ranLocMats))

## control
# ranLoc
control_ranLocMats <- mclapply(seq_along(controlNames), function(x) {
  lapply(seq_along(chrName), function(y) {
    if( grepl("map_K", controlNames[x]) ) {
      as.matrix(read.table(paste0(controlDirs[x],
                                  controlNames[x],
                                  "_MappedOn_T2T_Col_CEN180_in_",
                                  chrName[y], "_ranLoc_matrix_bin", binName, "_flank", flankName, ".tab"),
                           header = F, skip = 3))
    } else {
      as.matrix(read.table(paste0(controlDirs[x],
                                  controlNames[x],
                                  "_MappedOn_T2T_Col_lowXM_both_sort_norm_CEN180_in_",
                                  chrName[y], "_ranLoc_matrix_bin", binName, "_flank", flankName, ".tab"),
                           header = F, skip = 3))
      }
  })
}, mc.cores = length(controlNames))
# If ranLocs from multiple chromosomes are to be analysed,
# concatenate the corresponding ranLoc coverage matrices
control_ranLocMats <- mclapply(seq_along(control_ranLocMats), function(x) {
  if(length(chrName) > 1) {
    do.call(rbind, control_ranLocMats[[x]])
  } else {
    control_ranLocMats[[x]][[1]]
  }
}, mc.cores = length(control_ranLocMats))

# Conditionally calculate log2(ChIP/control)
# for each matrix depending on library
# ranLoc
log2ChIP_ranLocMats <- mclapply(seq_along(ChIP_ranLocMats), function(x) {
  if ( grepl("SRR529854", ChIPNames[x]) ) {
    print(paste0(ChIPNames[x], " library; using ", controlNames[2], " for log2((MNase+1)/(gDNA+1)) calculation"))
    log2((ChIP_ranLocMats[[x]]+1)/(control_ranLocMats[[2]]+1))
  } else if ( grepl("MNase", ChIPNames[x]) ) {
    print(paste0(ChIPNames[x], " library; using ", controlNames[3], " for log2((SPO11-1-oligos+1)/(gDNA+1)) calculation"))
    log2((ChIP_ranLocMats[[x]]+1)/(control_ranLocMats[[3]]+1))
  } else if ( grepl("SPO11oligos", ChIPNames[x]) ) {
    print(paste0(ChIPNames[x], " library; using ", controlNames[4], " for log2((SPO11-1-oligos+1)/(gDNA+1)) calculation"))
    log2((ChIP_ranLocMats[[x]]+1)/(control_ranLocMats[[4]]+1))
  } else {
    print(paste0(ChIPNames[x], " library; using ", controlNames[1], " for log2((ChIP+1)/(input+1)) calculation"))
    log2((ChIP_ranLocMats[[x]]+1)/(control_ranLocMats[[1]]+1))
  }
}, mc.cores = length(ChIP_ranLocMats))

# Calculate mean DNAmeth-dataset values for each ranLoc
DNAmeth_ranLocMats_bodies <- lapply(seq_along(DNAmeth_ranLocMats), function(x) {
  DNAmeth_ranLocMats[[x]][,((upstream/binSize)+1):((upstream+regionBodyLength)/binSize)]
})
DNAmeth_ranLocMats_bodiesRowMeans <- lapply(seq_along(DNAmeth_ranLocMats_bodies), function(x) {
  rowMeans(DNAmeth_ranLocMats_bodies[[x]], na.rm = T)
})

# Calculate mean log2(ChIP/control) coverage values for each ranLoc
log2ChIP_ranLocMats_bodies <- lapply(seq_along(log2ChIP_ranLocMats), function(x) {
  log2ChIP_ranLocMats[[x]][,((upstream/binSize)+1):((upstream+regionBodyLength)/binSize)]
})
log2ChIP_ranLocMats_bodiesRowMeans <- lapply(seq_along(log2ChIP_ranLocMats_bodies), function(x) {
  rowMeans(log2ChIP_ranLocMats_bodies[[x]], na.rm = T)
})

# Calculate mean control-dataset values for each ranLoc
control_ranLocMats_bodies <- lapply(seq_along(control_ranLocMats), function(x) {
  control_ranLocMats[[x]][,((upstream/binSize)+1):((upstream+regionBodyLength)/binSize)]
})
control_ranLocMats_bodiesRowMeans <- lapply(seq_along(control_ranLocMats_bodies), function(x) {
  rowMeans(control_ranLocMats_bodies[[x]], na.rm = T)
})

# Add mean coverage values to ranLoc dataframe
ranLoc <- data.frame(CEN180,
                     CENH3_in_bodies = log2ChIP_ranLocMats_bodiesRowMeans[[1]],
                     H3K9me2_in_bodies = log2ChIP_ranLocMats_bodiesRowMeans[[2]],
                     H3K27me1_in_bodies = log2ChIP_ranLocMats_bodiesRowMeans[[3]],
                     H3K4me1_in_bodies = log2ChIP_ranLocMats_bodiesRowMeans[[4]],
                     H3K4me2_in_bodies = log2ChIP_ranLocMats_bodiesRowMeans[[5]],
                     H3K4me3_in_bodies = log2ChIP_ranLocMats_bodiesRowMeans[[6]],
                     H2AW6_in_bodies = log2ChIP_ranLocMats_bodiesRowMeans[[7]],
                     H2AW7_in_bodies = log2ChIP_ranLocMats_bodiesRowMeans[[8]],
                     MNase_in_bodies = log2ChIP_ranLocMats_bodiesRowMeans[[9]],
                     REC8_HA_Rep2_in_bodies = log2ChIP_ranLocMats_bodiesRowMeans[[10]],
                     ASY1_Rep1_in_bodies = log2ChIP_ranLocMats_bodiesRowMeans[[11]],
                     MTOPVIB_HA_Rep1_in_bodies = log2ChIP_ranLocMats_bodiesRowMeans[[12]],
                     MTOPVIB_HA_Rep2_in_bodies = log2ChIP_ranLocMats_bodiesRowMeans[[13]],
                     DMC1_V5_Rep1_in_bodies = log2ChIP_ranLocMats_bodiesRowMeans[[14]],
                     DMC1_V5_Rep2_in_bodies = log2ChIP_ranLocMats_bodiesRowMeans[[15]],
                     SPO11oligos_in_bodies = log2ChIP_ranLocMats_bodiesRowMeans[[16]],
                     mCG_Nanopolish_in_bodies = DNAmeth_ranLocMats_bodiesRowMeans[[1]], 
                     mCG_PE_BSseq_in_bodies = DNAmeth_ranLocMats_bodiesRowMeans[[2]], 
                     mCHG_PE_BSseq_in_bodies = DNAmeth_ranLocMats_bodiesRowMeans[[3]], 
                     mCHH_PE_BSseq_in_bodies = DNAmeth_ranLocMats_bodiesRowMeans[[4]], 
                     mCG_SE_BSseq_in_bodies = DNAmeth_ranLocMats_bodiesRowMeans[[5]], 
                     mCHG_SE_BSseq_in_bodies = DNAmeth_ranLocMats_bodiesRowMeans[[6]], 
                     mCHH_SE_BSseq_in_bodies = DNAmeth_ranLocMats_bodiesRowMeans[[7]], 
                     PE_input_sonic_in_bodies = control_ranLocMats_bodiesRowMeans[[1]],
                     PE_input_MNase_in_bodies = control_ranLocMats_bodiesRowMeans[[2]],
                     PE_gDNA_in_bodies = control_ranLocMats_bodiesRowMeans[[3]],
                     SE_gDNA_in_bodies = control_ranLocMats_bodiesRowMeans[[4]],
                     map_K40_E2_in_bodies = control_ranLocMats_bodiesRowMeans[[5]],
                     map_K45_E2_in_bodies = control_ranLocMats_bodiesRowMeans[[6]],
                     map_K50_E2_in_bodies = control_ranLocMats_bodiesRowMeans[[7]],
                     map_K150_E4_in_bodies = control_ranLocMats_bodiesRowMeans[[8]],
                     map_K200_E4_in_bodies = control_ranLocMats_bodiesRowMeans[[9]],
                     map_K300_E4_in_bodies = control_ranLocMats_bodiesRowMeans[[10]])

# Define set of ordering factors to be used for grouping genes into quantiles
orderingFactor <- colnames(CEN180)[c(5, 7, 8, 11, 12, 28, 29, 35:37, 39:length(colnames(CEN180)))]
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
  ## Assign 0s to NA values only for coverage data
  #if(grepl("_in_", orderingFactor[w])) {
  #  CEN180_DF[,which(colnames(CEN180_DF) == orderingFactor[w])][
  #    which(is.na(CEN180_DF[,which(colnames(CEN180_DF) == orderingFactor[w])]))] <- 0
  #}
  if(orderingFactor[w] == "array_size" & chrName[length(chrName)] == "Chr1") {
    quantiles <- 2
  } else if(orderingFactor[w] == "array_size" & chrName[length(chrName)] == "Chr2") {
    quantiles <- 3
  } 
  quantilesStats <- data.frame()
  for(k in 1:quantiles) {
    # First quantile should span 1 to greater than, e.g., 0.75 proportions of features
    if(k < quantiles) {
      CEN180_DF[ !is.na(CEN180_DF[,which(colnames(CEN180_DF) == orderingFactor[w])]) &
                 rank(CEN180_DF[,which(colnames(CEN180_DF) == orderingFactor[w])]) /
                 length(CEN180_DF[,which(colnames(CEN180_DF) == orderingFactor[w])]) <=
                 1-((k-1)/quantiles) &
                 rank(CEN180_DF[,which(colnames(CEN180_DF) == orderingFactor[w])]) /
                 length(CEN180_DF[,which(colnames(CEN180_DF) == orderingFactor[w])]) >
                 1-(k/quantiles), ]$quantile <- paste0("Quantile ", k)
    } else {
    # Final quantile should span 0 to, e.g., 0.25 proportions of features
      CEN180_DF[ !is.na(CEN180_DF[,which(colnames(CEN180_DF) == orderingFactor[w])]) &
                 rank(CEN180_DF[,which(colnames(CEN180_DF) == orderingFactor[w])]) /
                 length(CEN180_DF[,which(colnames(CEN180_DF) == orderingFactor[w])]) <=
                 1-((k-1)/quantiles) &
                 rank(CEN180_DF[,which(colnames(CEN180_DF) == orderingFactor[w])]) /
                 length(CEN180_DF[,which(colnames(CEN180_DF) == orderingFactor[w])]) >=
                 1-(k/quantiles), ]$quantile <- paste0("Quantile ", k)
    }
    write.table(CEN180_DF[CEN180_DF$quantile == paste0("Quantile ", k),],
                file = paste0(outDir[w],
                              "quantile", k, "_of_", quantiles,
                              "_by_", orderingFactor[w],
                              "_of_CEN180_in_T2T_Col_",
                              paste0(chrName, collapse = "_"), ".tsv"),
                quote = FALSE, sep = "\t", row.names = FALSE)
    stats <- data.frame(quantile = as.integer(k),
                        n = as.integer(dim(CEN180_DF[CEN180_DF$quantile == paste0("Quantile ", k),])[1]),
                        mean_width = as.integer(round(mean(
                          (CEN180_DF[CEN180_DF$quantile == paste0("Quantile ", k),]$end -
                           CEN180_DF[CEN180_DF$quantile == paste0("Quantile ", k),]$start) + 1, na.rm = T))),
                        total_width = as.integer(sum(
                          (CEN180_DF[CEN180_DF$quantile == paste0("Quantile ", k),]$end -
                           CEN180_DF[CEN180_DF$quantile == paste0("Quantile ", k),]$start) + 1, na.rm = T)),
                        mean_orderingFactor = as.numeric(mean(CEN180_DF[CEN180_DF$quantile == paste0("Quantile ", k),][,which(colnames(CEN180_DF) == orderingFactor[w])], na.rm = T)))
    quantilesStats <- rbind(quantilesStats, stats)
  }
  write.table(quantilesStats,
              file = paste0(outDir[w],
                            "summary_", quantiles, "quantiles",
                            "_by_", orderingFactor[w],
                            "_of_CEN180_in_T2T_Col_",
                            paste0(chrName, collapse = "_"), ".tsv"),
              quote = FALSE, sep = "\t", row.names = FALSE)
  write.table(CEN180_DF,
              file = paste0(outDir[w],
                            "features_", quantiles, "quantiles",
                            "_by_", orderingFactor[w],
                            "_of_CEN180_in_T2T_Col_",
                            paste0(chrName, collapse = "_"), ".tsv"),
              quote = FALSE, sep = "\t", row.names = FALSE)

  # Divide ranLoc into quantiles based on feature quantile indices
  ranLoc_DF <- data.frame(ranLoc,
                          random = as.character(""))
  # Get row indices for each feature quantile
  quantileIndices <- lapply(1:quantiles, function(k) {
    which(CEN180_DF$quantile == paste0("Quantile ", k))
  })
  for(k in 1:quantiles) {
    ranLoc_DF[quantileIndices[[k]],]$random <- paste0("Random ", k)
  }
  write.table(ranLoc_DF,
              file = paste0(outDir[w],
                            "features_", quantiles, "quantiles",
                            "_by_", orderingFactor[w],
                            "_of_CEN180_in_T2T_Col_",
                            paste0(chrName, collapse = "_"), "_ranLoc.tsv"),
              quote = FALSE, sep = "\t", row.names = FALSE)
}, mc.cores = length(orderingFactor), mc.preschedule = F)


# Make correlation matrix (colour-gradient heatmap)
# Combine profiles into one data.frame in which each profile is a column
profilesVal <- lapply(seq_along(profilesGR), function(x) {
  profilesGR[[x]]$value
})

profilesDF <- as.data.frame(do.call(cbind, profilesVal),
                            stringsAsFactors = F)
colnames(profilesDF) <- profileNames

# Create correlation matrix
corMat <- round(cor(profilesDF,
                    method = "spearman",
                    use = "pairwise.complete.obs"),
                digits = 2)

# Set duplicates to NA
for(x in 1:dim(corMat)[1]) {
  corMat[x, x] <- NA
  if(x > 1) {
    corMat[x, 1:x-1] <- NA
  }
}
corMat <- corMat[,-1]

# Convert into reshape::melt formatted data.frame
# and remove duplicate pairs
corDat <- melt(corMat)
corDat <- corDat[-which(is.na(corDat[,3])),]

# Order the data.frame for plotting
profileNamesList <- as.list(profileNames)
names(profileNamesList) <- profileNames
levels(corDat$X1) <- rev(profileNamesList)
levels(corDat$X2) <- profileNamesList[-1]

# Get P-values for correlation matrix
corMatSig <- rcorr(as.matrix(profilesDF),
                   type = "spearman")$P
# Set duplicates to NA
for(x in 1:dim(corMatSig)[1]) {
  corMatSig[x, x] <- NA
  if(x > 1) {
    corMatSig[x, 1:x-1] <- NA
  }
}
corMatSig <- corMatSig[,-1]

# Convert into reshape::melt formatted data.frame
# and remove duplicate pairs
corDatSig <- melt(corMatSig)
corDatSig <- corDatSig[-which(is.na(corDatSig[,3])),]

# Standardise P-values to a sample size of 100 (q-values) as proposed by
# Good (1982) Standardized tail-area probabilities. Journal of Computation and Simulation 16: 65-66
# and summarised by Woolley (2003):
# https://stats.stackexchange.com/questions/22233/how-to-choose-significance-level-for-a-large-data-set
# http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.518.5341&rep=rep1&type=pdf
# Woolley (2003): "Clearly, the meaningfulness of the p-value diminishes as the sample size increases";
# Anne Z. (2012, Pearson eCollege, Denver): "In the real world, there are unlikely to be semi-partial correlations
# that are exactly zero, which is the null hypothesis in testing significance of a regression coefficient."
# Formally, the standardised p-value is defined as:
# q = min(0.5, p * sqrt( (n/100) ))
# Woolley (2003): "The value of 0.5 is somewhat arbitrary, though its purpose is to avoid q-values of greater than 1."
n <- dim(profilesDF)[1]
corDatSig$value <- sapply(corDatSig$value, function(x) {
  round(min(0.5, x * sqrt( (n/100) )),
        digits = 2)
})

# Order the data.frame for plotting
levels(corDatSig$X1) <- rev(profileNamesList)
levels(corDatSig$X2) <- profileNamesList[-1]

# Plot
ggObj <- ggplot(data = corDat,
                mapping = aes(X2, X1, fill = value)) +
  geom_tile() +
#  geom_text(mapping = aes(X2, X1, label = value), size = 5) +
  geom_text(data = corDatSig,
            mapping = aes(X2, X1, label = value), size = 8) +
  scale_fill_gradient2(name = bquote("Spearman's" ~ italic(r[s])),
                       low = "blue", mid = "white", high = "red",
                       midpoint = 0, breaks = seq(-1, 1, by = 0.4), limits = c(-1, 1)) +
  scale_x_discrete(expand = c(0, 0), position = "top") +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = "", y = "") +
  guides(fill = guide_colourbar(barwidth = 40, barheight = 6,
                                title.position = "top", title.hjust = 0.5)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0, size = 39, colour = "black"),
        axis.text.y = element_text(angle = 0, vjust = 1, hjust = 1, size = 39, colour = "black"),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_text(size = 40),
        legend.text = element_text(size = 40),
        legend.justification = c(1, 0),
        legend.position = c(0.65, 0.05),
        legend.direction = "horizontal",
        legend.background = element_rect(fill = "transparent"),
        plot.margin = unit(c(5.5, 70.5, 5.5, 5.5), "pt"),
        plot.title = element_text(hjust = 0.5, size = 30, colour = "black")) +
  ggtitle(bquote(.(winSize/1e6) * "-Mb Spearman's" ~ italic(r[s]) ~ "for" ~
          .(paste0(genomeName, collapse = "-, ")) * "-genome" ~
          .(region) ~ "regions (" * .(smoothing) * ")"))
ggsave(paste0(plotDir,
              "Spearman_correlation_matrix_", winName,
              "_log2ChIPcontrol_cMMb_MNase_DNAmeth_genes_TEsuperfams_in_",
              paste0(genomeName, collapse = "_"), "_genome_", region, "_", smoothing, "_qVals.pdf"),
       plot = ggObj, height = 20, width = 20)

