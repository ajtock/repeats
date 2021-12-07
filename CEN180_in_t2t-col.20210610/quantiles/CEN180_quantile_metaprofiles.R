#!/applications/R/R-4.0.0/bin/Rscript

# author: Andy Tock
# contact: ajt200@cam.ac.uk
# date: 07.12.2021

# Calculate and plot metaprofiles of ChIP-seq, MNase-seq, etc.
# (CEN180 windowed means and 95% confidence intervals, CIs)
# for each group of CEN180s, defined either by
# decreasing orderingFactor or randomly

# Usage:
# /applications/R/R-4.0.0/bin/Rscript CEN180_quantile_metaprofiles.R 'Chr1,Chr2,Chr3,Chr4,Chr5' wSNV 4 both 180 2000 2kb 10 '0.02,0.96'

#chrName <- unlist(strsplit("Chr1,Chr2,Chr3,Chr4,Chr5",
#                           split = ","))
#orderingFactor <- "CENH3_in_bodies"
#quantiles <- 4
#align <- "both"
#bodyLength <- 180
#upstream <- 2000
#downstream <- 2000
#flankName <- "2kb"
#binSize <- 10
## top left
#legendPos <- as.numeric(unlist(strsplit("0.02,0.96",
#                                        split = ",")))
## top centre
#legendPos <- as.numeric(unlist(strsplit("0.38,0.96",
#                                        split = ",")))
## top right
#legendPos <- as.numeric(unlist(strsplit("0.75,0.96",
#                                        split = ",")))
## bottom left
#legendPos <- as.numeric(unlist(strsplit("0.02,0.30",
#                                        split = ",")))

args <- commandArgs(trailingOnly = T)
chrName <- unlist(strsplit(args[1],
                           split = ","))
orderingFactor <- args[2]
quantiles <- as.numeric(args[3])
align <- args[4]
bodyLength <- as.numeric(args[5])
upstream <- as.numeric(args[6])
downstream <- as.numeric(args[6])
flankName <- args[7]
binSize <- as.numeric(args[8])
legendPos <- as.numeric(unlist(strsplit(args[9],
                                        split = ",")))

options(stringsAsFactors = F)
library(parallel)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(grid)
library(gridExtra)
library(extrafont)
library(viridis)

outDir <- paste0("quantiles_by_", orderingFactor, "/",
                 paste0(chrName, collapse = "_"), "/")
plotDir <- paste0(outDir, "plots/")
system(paste0("[ -d ", outDir, " ] || mkdir -p ", outDir))
system(paste0("[ -d ", plotDir, " ] || mkdir -p ", plotDir))

# Define plot titles
if(grepl("_in_", orderingFactor)) {
  featureNamePlot <- paste0(sub("_in_\\w+", "", orderingFactor), " quantiles")
} else if(grepl("SNV", orderingFactor)) {
  featureNamePlot <- "Variant quantiles"
} else if(orderingFactor == "array_size") {
  featureNamePlot <- "Array-size quantiles"
} else if(orderingFactor == "HORlengthsSum") {
  featureNamePlot <- "Repetitiveness quantiles"
} else if(orderingFactor == "HORcount") {
  featureNamePlot <- "HORcount quantiles"
} else if(orderingFactor == "HORavgSize") {
  featureNamePlot <- "HOR size quantiles"
}
ranFeatNamePlot <- "Random quantiles"
ranLocNamePlot <- "RanLoc quantiles"

# Define quantile colours
#quantileColours <- c("red", "purple", "blue", "navy")
quantileColours <- rev(viridis(quantiles))
quantileColours[1] <- "orange"
#quantileColours <- rev(plasma(quantiles))
#quantileColours[1] <- "gold"

# Define feature start and end labels for plotting
featureStartLab <- "Start"
featureEndLab <- "End"

# Load table of features grouped into quantiles
featuresDF <- read.table(paste0(outDir,
                                "features_", quantiles, "quantiles",
                                "_by_", orderingFactor,
                                "_of_CEN180_in_t2t-col.20210610_",
                                paste0(chrName, collapse = "_"), ".tsv"),
                         header = T, sep = "\t", row.names = NULL)

# Load features to confirm feature (row) ordering in "featuresDF" is the same
# as in "features" (which was used for generating the coverage matrices)
features <- lapply(seq_along(chrName), function(y) {
  read.table(paste0("/home/ajt200/analysis/nanopore/t2t-col.20210610/annotation/CEN180/",
                    "CEN180_in_t2t-col.20210610_",
                    chrName[y], ".bed"),
             header = F)
})
# If features from multiple subgenomes and/or compartments are to be analysed,
# concatenate the corresponding feature data.frames
if(length(chrName) > 1) {
  features <- do.call(rbind, features)
} else {
  features <- features[[1]]
}
colnames(features) <- c("chr", "start", "end", "featureID", "wSNV", "strand", "HORlengthsSum", "HORcount")
stopifnot(identical(as.character(featuresDF$featureID),
                    as.character(features$featureID)))
rm(features); gc()

# Get row indices for each feature quantile
quantileIndices <- lapply(1:quantiles, function(k) {
  which(featuresDF$quantile == paste0("Quantile ", k))
})

## Random feature quantiles
# Define function to randomly select n rows from
# a data.frame
selectRandomFeatures <- function(features, n) {
  return(features[sample(x = dim(features)[1],
                         size = n,
                         replace = FALSE),])
}

# Divide features into random sets of equal number,
# with the same number of CEN180s per chromosome as
# above-defined orderingFactor-defined feature quantiles
randomPCIndices <- lapply(1:quantiles, function(k) {
  randomPCIndicesk <- NULL
  for(i in 1:length(chrName)) {
    # Define seed so that random selections are reproducible
    set.seed(93750174)
    randomPCfeatureskChr <- selectRandomFeatures(features = featuresDF[featuresDF$chr == chrName[i],],
                                                 n = dim(featuresDF[featuresDF$quantile == paste0("Quantile ", k) &
                                                                    featuresDF$chr == chrName[i],])[1])
    randomPCIndicesk <- c(randomPCIndicesk, as.integer(rownames(randomPCfeatureskChr)))
  }
  randomPCIndicesk
})
# Confirm per-chromosome feature numbers are the same for quantiles and random groupings
lapply(seq_along(1:quantiles), function(k) {
  sapply(seq_along(chrName), function(x) {
    if(!identical(dim(featuresDF[randomPCIndices[[k]],][featuresDF[randomPCIndices[[k]],]$chr == chrName[x],]),
                  dim(featuresDF[quantileIndices[[k]],][featuresDF[quantileIndices[[k]],]$chr == chrName[x],])))    {
      stop("Quantile features and random features do not consist of the same number of features per chromosome")
    }
  })
})


# Load feature matrices for each chromatin dataset, calculate log2(ChIP/control),
# and sort by decreasing log2mat1RegionRowMeans
# Load feature matrices for CENH3 and input, calculate log2(ChIP/control) coverage
ChIPNames <- c(
               "WT_CENH3_Rep1_ChIP_SRR4430537",
               "Col_CENH3_PE150_Rep1_ChIP",
               "met1_CENH3_PE150_Rep1_ChIP",
               "Col_CENH3_PE250_Rep1_ChIP",
               "met1_CENH3_PE250_Rep1_ChIP",
               "cmt23_CENH3_PE250_Rep1_ChIP",
               "ddc_CENH3_PE250_Rep1_ChIP",
               "kss_CENH3_PE250_Rep1_ChIP"
#               "WT_H3K9me2_Rep1_ChIP",
#               "WT_H3K27me1_Rep1_ChIP",
#               "WT_H3K4me1_Rep1_ChIP",
#               "WT_H3K4me2_Rep1_ChIP",
#               "WT_H3K4me3_ChIP14",
#               "H2AW6_ChIP_SRR5298545",
#               "H2AW7_ChIP_SRR5298546",
#               "WT_MNase_Rep1",
#               "WT_REC8_HA_Rep2_ChIP",
#               "WT_ASY1_Rep1_ChIP",
#               "WT_MTOPVIB_HA_Rep1_ChIP",
#               "WT_MTOPVIB_HA_Rep2_ChIP",
#               "WT_DMC1_V5_Rep1_ChIP",
#               "WT_DMC1_V5_Rep2_ChIP",
#               "WT_SPO11oligos_Rep1"
              )
ChIPNamesDir <- c(
                  "CENH3_seedlings_Maheshwari_Comai_2017_GenomeRes/snakemake_ChIPseq_t2t-col.20210610",
                  rep("CENH3_PE150_mn359_20210818/snakemake_ChIPseq_t2t-col.20210610", 2),
                  rep("CENH3_PE250_mn359_20210908/snakemake_ChIPseq_t2t-col.20210610", 5)
#                  "170101_Chris_H3K9me2_ChIP/WT/snakemake_ChIPseq_t2t-col.20210610",
#                  rep("170920_Chris_histone_ChIP/snakemake_ChIPseq_t2t-col.20210610", 3),
#                  "160601_Kyuha_H3K4me3_ChIP/WT/snakemake_ChIPseq_t2t-col.20210610",
#                  rep("HTA6_HTA7_leaf_Lorkovic_Berger_2017_CurrBiol/snakemake_ChIPseq_t2t-col.20210610", 2),
#                  "150701_Kyuha_MNase/WT/snakemake_ChIPseq_t2t-col.20210610",
#                  "REC8_pooled/snakemake_ChIPseq_t2t-col.20210610",
#                  "20190722_cal66_Athaliana_ChIPseq_ASY1/fastq_pooled/snakemake_ChIPseq_t2t-col.20210610",
#                  rep("20190819_dh580_Athaliana_ChIPseq_MTOPVIB/fastq_pooled/snakemake_ChIPseq_t2t-col.20210610", 2),
#                  rep("20190917_dh580_Athaliana_ChIPseq_DMC1/fastq_pooled/snakemake_ChIPseq_t2t-col.20210610", 2),
#                  "160518_Kyuha_SPO11oligos/WT/snakemake_SPO11oligos_t2t-col.20210610"
                 )
log2ChIPNamesPlot <- c(
                       "Col CENH3 (2×100 bp)",
                       "Col CENH3 (2×150 bp)",
                       "met1 CENH3 (2×150 bp)",
                       "Col CENH3 (2×250 bp)",
                       "met1 CENH3 (2×250 bp)",
                       "cmt23 CENH3 (2×250 bp)",
                       "ddc CENH3 (2×250 bp)",
                       "kss CENH3 (2×250 bp)"
#                       "H3K9me2",
#                       "H3K27me1",
#                       "H3K4me1",
#                       "H3K4me2",
#                       "H3K4me3",
#                       "H2A.W.6",
#                       "H2A.W.7",
#                       "MNase",
#                       "REC8",
#                       "ASY1",
#                       "MTOPVIB Rep1",
#                       "MTOPVIB Rep2",
#                       "DMC1 Rep1",
#                       "DMC1 Rep2",
#                       "SPO11-1"
                      )
ChIPNamesPlot <- log2ChIPNamesPlot
ChIPDirs <- sapply(seq_along(ChIPNamesDir), function(x) {
  paste0("/home/ajt200/analysis/", ChIPNamesDir[x],
         "/mapped/CEN180profiles/matrices/")
})
log2ChIPColours <- c(
                     rep("black", length(log2ChIPNamesPlot))
                    )
ChIPColours <- log2ChIPColours

controlNames <- c(
                  "WT_CENH3_Rep1_input_SRR4430555",
                  "Col_CENH3_PE150_Rep1_input",
                  "met1_CENH3_PE150_Rep1_input",
                  "Col_CENH3_PE250_Rep1_input",
                  "cmt23_CENH3_PE250_Rep1_input",
                  "ddc_CENH3_PE250_Rep1_input",
                  "kss_CENH3_PE250_Rep1_input"
#                  "WT_REC8_Myc_Rep1_input",
#                  "H2AW_input_SRR5298544",
#                  "WT_gDNA_Rep1",
#                  "WT_gDNA_Rep1_R1",
#                  "map_K40_E2",
#                  "map_K45_E2",
#                  "map_K50_E2",
#                  "map_K150_E4",
#                  "map_K200_E4",
#                  "map_K300_E4"
                 )
controlNamesDir <- c(
                     "CENH3_seedlings_Maheshwari_Comai_2017_GenomeRes/snakemake_ChIPseq_t2t-col.20210610",
                     rep("CENH3_PE150_mn359_20210818/snakemake_ChIPseq_t2t-col.20210610", 2),
                     rep("CENH3_PE250_mn359_20210908/snakemake_ChIPseq_t2t-col.20210610", 4)
#                     "REC8_pooled/snakemake_ChIPseq_t2t-col.20210610",
#                     "HTA6_HTA7_leaf_Lorkovic_Berger_2017_CurrBiol/snakemake_ChIPseq_t2t-col.20210610",
#                     "150701_Natasha_gDNA/WT/snakemake_ChIPseq_t2t-col.20210610",
#                     "150701_Natasha_gDNA/WT/R1/snakemake_SPO11oligos_t2t-col.20210610",
#                     rep("nanopore/t2t-col.20210610/genmap_mappability", 6)
                    )
controlNamesPlot <- c(
                      "Col input (2×100 bp)",
                      "Col input (2×150 bp)",
                      "met1 input (2×150 bp)",
                      "Col input (2×250 bp)",
                      "cmt23 input (2×250 bp)",
                      "ddc input (2×250 bp)",
                      "kss input (2×250 bp)"
#                      "Input (sonic.)",
#                      "Input (MNase)",
#                      "PE gDNA",
#                      "SE gDNA",
#                      "k=40 e=2",
#                      "k=45 e=2",
#                      "k=50 e=2",
#                      "k=150 e=4",
#                      "k=200 e=4",
#                      "k=300 e=4"
                     )
controlDirs <- sapply(seq_along(controlNamesDir), function(x) {
  paste0("/home/ajt200/analysis/", controlNamesDir[x],
         "/mapped/CEN180profiles/matrices/")
})
controlColours <- c(
                    rep("black", length(controlNamesPlot))
                   )

## ChIP
# feature
ChIP_featureMats <- mclapply(seq_along(ChIPNames), function(x) {
  lapply(seq_along(chrName), function(y) {
    as.matrix(read.table(paste0(ChIPDirs[x],
                                ChIPNames[x],
                                "_MappedOn_t2t-col.20210610_lowXM_", align, "_sort_norm_CEN180_in_",
                                chrName[y], "_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                         header = F, skip = 3,
                         colClasses = c(rep("NULL", 1000/binSize), rep(NA, ((1000*2)+(180))/binSize), rep("NULL", 1000/binSize))))
  })
}, mc.cores = length(ChIPNames))
# If features from multiple chromosomes are to be analysed,
# concatenate the corresponding feature coverage matrices
ChIP_featureMats <- mclapply(seq_along(ChIP_featureMats), function(x) {
  if(length(chrName) > 1) {
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
                                  "_MappedOn_t2t-col.20210610_CEN180_in_",
                                  chrName[y], "_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                           header = F, skip = 3,
                           colClasses = c(rep("NULL", 1000/binSize), rep(NA, ((1000*2)+(180))/binSize), rep("NULL", 1000/binSize))))
    } else {
      as.matrix(read.table(paste0(controlDirs[x],
                                  controlNames[x],
                                  "_MappedOn_t2t-col.20210610_lowXM_", align, "_sort_norm_CEN180_in_",
                                  chrName[y], "_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                           header = F, skip = 3,
                           colClasses = c(rep("NULL", 1000/binSize), rep(NA, ((1000*2)+(180))/binSize), rep("NULL", 1000/binSize))))
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

## ChIP
# ranLoc
ChIP_ranLocMats <- mclapply(seq_along(ChIPNames), function(x) {
  lapply(seq_along(chrName), function(y) {
    as.matrix(read.table(paste0(ChIPDirs[x],
                                ChIPNames[x],
                                "_MappedOn_t2t-col.20210610_lowXM_", align, "_sort_norm_CEN180_in_",
                                chrName[y], "_ranLoc_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                         header = F, skip = 3,
                         colClasses = c(rep("NULL", 1000/binSize), rep(NA, ((1000*2)+(180))/binSize), rep("NULL", 1000/binSize))))
  })
}, mc.cores = length(ChIPNames))
# If ranLocs from multiple chromosomes are to be analysed,
# concatenate the corresponding ranLoc coverage matrices
ChIP_ranLocMats <- mclapply(seq_along(ChIP_ranLocMats), function(x) {
  if(length(chrName) > 1) {
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
                                  "_MappedOn_t2t-col.20210610_CEN180_in_",
                                  chrName[y], "_ranLoc_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                           header = F, skip = 3,
                           colClasses = c(rep("NULL", 1000/binSize), rep(NA, ((1000*2)+(180))/binSize), rep("NULL", 1000/binSize))))
    } else {
      as.matrix(read.table(paste0(controlDirs[x],
                                  controlNames[x],
                                  "_MappedOn_t2t-col.20210610_lowXM_", align, "_sort_norm_CEN180_in_",
                                  chrName[y], "_ranLoc_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                           header = F, skip = 3,
                           colClasses = c(rep("NULL", 1000/binSize), rep(NA, ((1000*2)+(180))/binSize), rep("NULL", 1000/binSize))))
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
# feature
log2ChIP_featureMats <- mclapply(seq_along(ChIP_featureMats), function(x) {
  if ( ChIPNames[x] %in% c("WT_CENH3_Rep1_ChIP_SRR4430537") ) {
    print(paste0(ChIPNames[x], " library; using ", controlNames[1], " for log2((ChIP+1)/(input+1)) calculation"))
    log2((ChIP_featureMats[[x]]+1)/(control_featureMats[[1]]+1))
  } else if ( ChIPNames[x] %in% c("Col_CENH3_PE150_Rep1_ChIP") ) {
    print(paste0(ChIPNames[x], " library; using ", controlNames[2], " for log2((ChIP+1)/(input+1)) calculation"))
    log2((ChIP_featureMats[[x]]+1)/(control_featureMats[[2]]+1))
  } else if ( ChIPNames[x] %in% c("met1_CENH3_PE150_Rep1_ChIP") ) {
    print(paste0(ChIPNames[x], " library; using ", controlNames[3], " for log2((ChIP+1)/(input+1)) calculation"))
    log2((ChIP_featureMats[[x]]+1)/(control_featureMats[[3]]+1))
  } else if ( ChIPNames[x] %in% c("Col_CENH3_PE250_Rep1_ChIP", "met1_CENH3_PE250_Rep1_ChIP") ) {
    print(paste0(ChIPNames[x], " library; using ", controlNames[4], " for log2((ChIP+1)/(input+1)) calculation"))
    log2((ChIP_featureMats[[x]]+1)/(control_featureMats[[4]]+1))
  } else if ( ChIPNames[x] %in% c("cmt23_CENH3_PE250_Rep1_ChIP") ) {
    print(paste0(ChIPNames[x], " library; using ", controlNames[5], " for log2((ChIP+1)/(input+1)) calculation"))
    log2((ChIP_featureMats[[x]]+1)/(control_featureMats[[5]]+1))
  } else if ( ChIPNames[x] %in% c("ddc_CENH3_PE250_Rep1_ChIP") ) {
    print(paste0(ChIPNames[x], " library; using ", controlNames[6], " for log2((ChIP+1)/(input+1)) calculation"))
    log2((ChIP_featureMats[[x]]+1)/(control_featureMats[[6]]+1))
  } else if ( ChIPNames[x] %in% c("kss_CENH3_PE250_Rep1_ChIP") ) {
    print(paste0(ChIPNames[x], " library; using ", controlNames[7], " for log2((ChIP+1)/(input+1)) calculation"))
    log2((ChIP_featureMats[[x]]+1)/(control_featureMats[[7]]+1))
  }
#  } else if ( grepl("MNaseX", ChIPNames[x]) ) {
#    print(paste0(ChIPNames[x], " library; using ", controlNames[3], " for log2((MNase+1)/(gDNA+1)) calculation"))
#    log2((ChIP_featureMats[[x]]+1)/(control_featureMats[[3]]+1))
#  } else if ( grepl("SPO11oligos", ChIPNames[x]) ) {
#    print(paste0(ChIPNames[x], " library; using ", controlNames[4], " for log2((SPO11-1-oligos+1)/(gDNA+1)) calculation"))
#    log2((ChIP_featureMats[[x]]+1)/(control_featureMats[[4]]+1))
#  } else {
#    print(paste0(ChIPNames[x], " library; using ", controlNames[1], " for log2((ChIP+1)/(input+1)) calculation"))
#    log2((ChIP_featureMats[[x]]+1)/(control_featureMats[[1]]+1))
#  }
}, mc.cores = length(ChIP_featureMats))

# Conditionally calculate log2(ChIP/control)
# for each matrix depending on library
# ranLoc
log2ChIP_ranLocMats <- mclapply(seq_along(ChIP_ranLocMats), function(x) {
  if ( ChIPNames[x] %in% c("WT_CENH3_Rep1_ChIP_SRR4430537") ) {
    print(paste0(ChIPNames[x], " library; using ", controlNames[1], " for log2((ChIP+1)/(input+1)) calculation"))
    log2((ChIP_ranLocMats[[x]]+1)/(control_ranLocMats[[1]]+1))
  } else if ( ChIPNames[x] %in% c("Col_CENH3_PE150_Rep1_ChIP") ) {
    print(paste0(ChIPNames[x], " library; using ", controlNames[2], " for log2((ChIP+1)/(input+1)) calculation"))
    log2((ChIP_ranLocMats[[x]]+1)/(control_ranLocMats[[2]]+1))
  } else if ( ChIPNames[x] %in% c("met1_CENH3_PE150_Rep1_ChIP") ) {
    print(paste0(ChIPNames[x], " library; using ", controlNames[3], " for log2((ChIP+1)/(input+1)) calculation"))
    log2((ChIP_ranLocMats[[x]]+1)/(control_ranLocMats[[3]]+1))
  } else if ( ChIPNames[x] %in% c("Col_CENH3_PE250_Rep1_ChIP", "met1_CENH3_PE250_Rep1_ChIP") ) {
    print(paste0(ChIPNames[x], " library; using ", controlNames[4], " for log2((ChIP+1)/(input+1)) calculation"))
    log2((ChIP_ranLocMats[[x]]+1)/(control_ranLocMats[[4]]+1))
  } else if ( ChIPNames[x] %in% c("cmt23_CENH3_PE250_Rep1_ChIP") ) {
    print(paste0(ChIPNames[x], " library; using ", controlNames[5], " for log2((ChIP+1)/(input+1)) calculation"))
    log2((ChIP_ranLocMats[[x]]+1)/(control_ranLocMats[[5]]+1))
  } else if ( ChIPNames[x] %in% c("ddc_CENH3_PE250_Rep1_ChIP") ) {
    print(paste0(ChIPNames[x], " library; using ", controlNames[6], " for log2((ChIP+1)/(input+1)) calculation"))
    log2((ChIP_ranLocMats[[x]]+1)/(control_ranLocMats[[6]]+1))
  } else if ( ChIPNames[x] %in% c("kss_CENH3_PE250_Rep1_ChIP") ) {
    print(paste0(ChIPNames[x], " library; using ", controlNames[7], " for log2((ChIP+1)/(input+1)) calculation"))
    log2((ChIP_ranLocMats[[x]]+1)/(control_ranLocMats[[7]]+1))
  }
#  } else if ( grepl("MNaseX", ChIPNames[x]) ) {
#    print(paste0(ChIPNames[x], " library; using ", controlNames[3], " for log2((MNase+1)/(gDNA+1)) calculation"))
#    log2((ChIP_ranLocMats[[x]]+1)/(control_ranLocMats[[3]]+1))
#  } else if ( grepl("SPO11oligos", ChIPNames[x]) ) {
#    print(paste0(ChIPNames[x], " library; using ", controlNames[4], " for log2((SPO11-1-oligos+1)/(gDNA+1)) calculation"))
#    log2((ChIP_ranLocMats[[x]]+1)/(control_ranLocMats[[4]]+1))
#  } else {
#    print(paste0(ChIPNames[x], " library; using ", controlNames[1], " for log2((ChIP+1)/(input+1)) calculation"))
#    log2((ChIP_ranLocMats[[x]]+1)/(control_ranLocMats[[1]]+1))
#  }
}, mc.cores = length(ChIP_ranLocMats))

# log2ChIP
# Add column names
for(x in seq_along(log2ChIP_featureMats)) {
  colnames(log2ChIP_featureMats[[x]]) <- c(paste0("u", 1:((upstream-1000)/binSize)),
                                           paste0("t", (((upstream-1000)/binSize)+1):(((upstream-1000)+bodyLength)/binSize)),
                                           paste0("d", ((((upstream-1000)+bodyLength)/binSize)+1):((((upstream-1000)+bodyLength)/binSize)+((downstream-1000)/binSize))))
  colnames(log2ChIP_ranLocMats[[x]]) <- c(paste0("u", 1:((upstream-1000)/binSize)),
                                          paste0("t", (((upstream-1000)/binSize)+1):(((upstream-1000)+bodyLength)/binSize)),
                                          paste0("d", ((((upstream-1000)+bodyLength)/binSize)+1):((((upstream-1000)+bodyLength)/binSize)+((downstream-1000)/binSize))))
}

# Subdivide coverage matrices into above-defined quantiles and random groupings
log2ChIP_mats_quantiles <- mclapply(seq_along(log2ChIP_featureMats), function(x) {
  list(
       # feature quantiles
       lapply(1:quantiles, function(k) {
         log2ChIP_featureMats[[x]][quantileIndices[[k]],]
       }),
       # feature random groupings
       lapply(1:quantiles, function(k) {
         log2ChIP_featureMats[[x]][randomPCIndices[[k]],]
       }),
       # random loci groupings
       lapply(1:quantiles, function(k) {
         log2ChIP_ranLocMats[[x]][quantileIndices[[k]],]
       })
      ) 
}, mc.cores = length(log2ChIP_featureMats))

# Transpose matrix and convert into dataframe
# in which first column is window name
wideDFfeature_list_log2ChIP <- mclapply(seq_along(log2ChIP_mats_quantiles), function(x) {
  lapply(seq_along(log2ChIP_mats_quantiles[[x]]), function(y) {
    lapply(seq_along(log2ChIP_mats_quantiles[[x]][[y]]), function(k) {
      data.frame(window = colnames(log2ChIP_mats_quantiles[[x]][[y]][[k]]),
                 t(log2ChIP_mats_quantiles[[x]][[y]][[k]]))
    })
  })
}, mc.cores = length(log2ChIP_mats_quantiles))

# Convert into tidy data.frame (long format)
tidyDFfeature_list_log2ChIP  <- mclapply(seq_along(wideDFfeature_list_log2ChIP), function(x) {
  lapply(seq_along(log2ChIP_mats_quantiles[[x]]), function(y) {
    lapply(seq_along(log2ChIP_mats_quantiles[[x]][[y]]), function(k) {
      gather(data  = wideDFfeature_list_log2ChIP[[x]][[y]][[k]],
             key   = feature,
             value = coverage,
             -window)
    })
  }) 
}, mc.cores = length(wideDFfeature_list_log2ChIP))

# Order levels of factor "window" so that sequential levels
# correspond to sequential windows
for(x in seq_along(tidyDFfeature_list_log2ChIP)) {
  for(y in seq_along(log2ChIP_mats_quantiles[[x]])) {
    for(k in seq_along(log2ChIP_mats_quantiles[[x]][[y]])) {
      tidyDFfeature_list_log2ChIP[[x]][[y]][[k]]$window <- factor(tidyDFfeature_list_log2ChIP[[x]][[y]][[k]]$window,
                                                                  levels = as.character(wideDFfeature_list_log2ChIP[[x]][[y]][[k]]$window))
    }
  }
}

# Create summary data.frame in which each row corresponds to a window (Column 1),
# Column2 is the number of coverage values (features) per window,
# Column3 is the mean of coverage values per window,
# Column4 is the standard deviation of coverage values per window,
# Column5 is the standard error of the mean of coverage values per window,
# Column6 is the lower bound of the 95% confidence interval, and
# Column7 is the upper bound of the 95% confidence interval
summaryDFfeature_list_log2ChIP  <- mclapply(seq_along(tidyDFfeature_list_log2ChIP), function(x) {
  lapply(seq_along(log2ChIP_mats_quantiles[[x]]), function(y) {
    lapply(seq_along(log2ChIP_mats_quantiles[[x]][[y]]), function(k) {
      data.frame(window = as.character(wideDFfeature_list_log2ChIP[[x]][[y]][[k]]$window),
                 n      = tapply(X     = tidyDFfeature_list_log2ChIP[[x]][[y]][[k]]$coverage,
                                 INDEX = tidyDFfeature_list_log2ChIP[[x]][[y]][[k]]$window,
                                 FUN   = length),
                 mean   = tapply(X     = tidyDFfeature_list_log2ChIP[[x]][[y]][[k]]$coverage,
                                 INDEX = tidyDFfeature_list_log2ChIP[[x]][[y]][[k]]$window,
                                 FUN   = mean,
                                 na.rm = TRUE),
                 sd     = tapply(X     = tidyDFfeature_list_log2ChIP[[x]][[y]][[k]]$coverage,
                                 INDEX = tidyDFfeature_list_log2ChIP[[x]][[y]][[k]]$window,
                                 FUN   = sd,
                                 na.rm = TRUE))
    })
  })
}, mc.cores = length(tidyDFfeature_list_log2ChIP))

for(x in seq_along(summaryDFfeature_list_log2ChIP)) {
  for(y in seq_along(log2ChIP_mats_quantiles[[x]])) {
    for(k in seq_along(log2ChIP_mats_quantiles[[x]][[y]])) {
      summaryDFfeature_list_log2ChIP[[x]][[y]][[k]]$window <- factor(summaryDFfeature_list_log2ChIP[[x]][[y]][[k]]$window,
                                                                     levels = as.character(wideDFfeature_list_log2ChIP[[x]][[y]][[k]]$window))
      summaryDFfeature_list_log2ChIP[[x]][[y]][[k]]$winNo <- factor(1:dim(summaryDFfeature_list_log2ChIP[[x]][[y]][[k]])[1])
      summaryDFfeature_list_log2ChIP[[x]][[y]][[k]]$sem <- summaryDFfeature_list_log2ChIP[[x]][[y]][[k]]$sd/sqrt(summaryDFfeature_list_log2ChIP[[x]][[y]][[k]]$n-1)
      summaryDFfeature_list_log2ChIP[[x]][[y]][[k]]$CI_lower <- summaryDFfeature_list_log2ChIP[[x]][[y]][[k]]$mean -
        qt(0.975, df = summaryDFfeature_list_log2ChIP[[x]][[y]][[k]]$n-1)*summaryDFfeature_list_log2ChIP[[x]][[y]][[k]]$sem
      summaryDFfeature_list_log2ChIP[[x]][[y]][[k]]$CI_upper <- summaryDFfeature_list_log2ChIP[[x]][[y]][[k]]$mean +
        qt(0.975, df = summaryDFfeature_list_log2ChIP[[x]][[y]][[k]]$n-1)*summaryDFfeature_list_log2ChIP[[x]][[y]][[k]]$sem
    }
  }
}

quantileNames <- paste0(rep("Quantile ", quantiles), 1:quantiles)
randomPCNames <- paste0(rep("Random ", quantiles), 1:quantiles)
for(x in seq_along(summaryDFfeature_list_log2ChIP)) {
  # feature quantiles
  names(summaryDFfeature_list_log2ChIP[[x]][[1]]) <- quantileNames
  # feature random groupings
  names(summaryDFfeature_list_log2ChIP[[x]][[2]]) <- randomPCNames
  # random loci groupings
  names(summaryDFfeature_list_log2ChIP[[x]][[3]]) <- randomPCNames
}

# Convert list of lists of lists of feature quantiles summaryDFfeature_list_log2ChIP into
# a list of lists of single data.frames containing all feature quantiles for plotting
summaryDFfeature_log2ChIP  <- mclapply(seq_along(summaryDFfeature_list_log2ChIP), function(x) {
  lapply(seq_along(log2ChIP_mats_quantiles[[x]]), function(y) {
    bind_rows(summaryDFfeature_list_log2ChIP[[x]][[y]], .id = "quantile")
  })
}, mc.cores = length(summaryDFfeature_list_log2ChIP))
for(x in seq_along(summaryDFfeature_log2ChIP)) {
  # feature quantiles
  summaryDFfeature_log2ChIP[[x]][[1]]$quantile <- factor(summaryDFfeature_log2ChIP[[x]][[1]]$quantile,
                                                         levels = names(summaryDFfeature_list_log2ChIP[[x]][[1]]))
  # feature random groupings
  summaryDFfeature_log2ChIP[[x]][[2]]$quantile <- factor(summaryDFfeature_log2ChIP[[x]][[2]]$quantile,
                                                         levels = names(summaryDFfeature_list_log2ChIP[[x]][[2]]))
  # random loci groupings
  summaryDFfeature_log2ChIP[[x]][[3]]$quantile <- factor(summaryDFfeature_log2ChIP[[x]][[3]]$quantile,
                                                         levels = names(summaryDFfeature_list_log2ChIP[[x]][[3]]))
}

# Define y-axis limits
ymin_list_log2ChIP <- min(unlist(lapply(seq_along(summaryDFfeature_log2ChIP), function(x) {
  min(c(summaryDFfeature_log2ChIP[[x]][[1]]$CI_lower,
        summaryDFfeature_log2ChIP[[x]][[2]]$CI_lower,
        summaryDFfeature_log2ChIP[[x]][[3]]$CI_lower))
})))
ymax_list_log2ChIP <- max(unlist(lapply(seq_along(summaryDFfeature_log2ChIP), function(x) {
  max(c(summaryDFfeature_log2ChIP[[x]][[1]]$CI_upper,
        summaryDFfeature_log2ChIP[[x]][[2]]$CI_upper,
        summaryDFfeature_log2ChIP[[x]][[3]]$CI_upper))
})))

# Define legend labels
legendLabs_feature <- lapply(seq_along(quantileNames), function(x) {
  grobTree(textGrob(bquote(.(quantileNames[x])),
                    x = legendPos[1], y = legendPos[2]-((x-1)*0.06), just = "left",
                    gp = gpar(col = quantileColours[x], fontsize = 18)))
})
legendLabs_ranFeat <- lapply(seq_along(randomPCNames), function(x) {
  grobTree(textGrob(bquote(.(randomPCNames[x])),
                    x = legendPos[1], y = legendPos[2]-((x-1)*0.06), just = "left",
                    gp = gpar(col = quantileColours[x], fontsize = 18)))
})
legendLabs_ranLoc <- lapply(seq_along(randomPCNames), function(x) {
  grobTree(textGrob(bquote(.(randomPCNames[x])),
                    x = legendPos[1], y = legendPos[2]-((x-1)*0.06), just = "left",
                    gp = gpar(col = quantileColours[x], fontsize = 18)))
})

# Plot average profiles with 95% CI ribbon
## feature
ggObj1_combined_log2ChIP <- mclapply(seq_along(log2ChIPNamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_log2ChIP[[x]][[1]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = quantile)
        ) +
  geom_line(data = summaryDFfeature,
            mapping = aes(colour = quantile),
            size = 1) +
  scale_colour_manual(values = quantileColours) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = quantile),
              alpha = 0.4) +
  scale_fill_manual(values = quantileColours) +
  scale_y_continuous(limits = c(ymin_list_log2ChIP, ymax_list_log2ChIP),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              ((upstream-1000)/binSize)+1,
                              (dim(summaryDFfeature_log2ChIP[[x]][[1]])[1]/quantiles)-((downstream-1000)/binSize),
                              dim(summaryDFfeature_log2ChIP[[x]][[1]])[1]/quantiles),
                   labels = c(paste0("-", "1kb"),
                              featureStartLab,
                              featureEndLab,
                              paste0("+", "1kb"))) +
  geom_vline(xintercept = c(((upstream-1000)/binSize)+1,
                            (dim(summaryDFfeature_log2ChIP[[x]][[1]])[1]/quantiles)-((downstream-1000)/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = log2ChIPNamesPlot[x]) +
  annotation_custom(legendLabs_feature[[1]]) +
  annotation_custom(legendLabs_feature[[2]]) +
  annotation_custom(legendLabs_feature[[3]]) +
  annotation_custom(legendLabs_feature[[4]]) +
  theme_bw() +
  theme(
        axis.ticks = element_line(size = 1.0, colour = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.x = element_text(size = 22, colour = "black"),
        axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
        axis.title = element_text(size = 30, colour = log2ChIPColours[x]),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_rect(size = 3.5, colour = "black"),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
        plot.title = element_text(hjust = 1.0, size = 30)) +
  ggtitle(bquote(.(featureNamePlot) ~ "(" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFfeature$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
}, mc.cores = length(log2ChIPNamesPlot))

## ranFeat
ggObj2_combined_log2ChIP <- mclapply(seq_along(log2ChIPNamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_log2ChIP[[x]][[2]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = quantile)
        ) +
  geom_line(data = summaryDFfeature,
            mapping = aes(colour = quantile),
            size = 1) +
  scale_colour_manual(values = quantileColours) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = quantile),
              alpha = 0.4) +
  scale_fill_manual(values = quantileColours) +
  scale_y_continuous(limits = c(ymin_list_log2ChIP, ymax_list_log2ChIP),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              ((upstream-1000)/binSize)+1,
                              (dim(summaryDFfeature_log2ChIP[[x]][[2]])[1]/quantiles)-((downstream-1000)/binSize),
                              dim(summaryDFfeature_log2ChIP[[x]][[2]])[1]/quantiles),
                   labels = c(paste0("-", "1kb"),
                              featureStartLab,
                              featureEndLab,
                              paste0("+", "1kb"))) +
  geom_vline(xintercept = c(((upstream-1000)/binSize)+1,
                            (dim(summaryDFfeature_log2ChIP[[x]][[2]])[1]/quantiles)-((downstream-1000)/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = log2ChIPNamesPlot[x]) +
  annotation_custom(legendLabs_ranFeat[[1]]) +
  annotation_custom(legendLabs_ranFeat[[2]]) +
  annotation_custom(legendLabs_ranFeat[[3]]) +
  annotation_custom(legendLabs_ranFeat[[4]]) +
  theme_bw() +
  theme(
        axis.ticks = element_line(size = 1.0, colour = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.x = element_text(size = 22, colour = "black"),
        axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
        axis.title = element_text(size = 30, colour = log2ChIPColours[x]),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_rect(size = 3.5, colour = "black"),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
        plot.title = element_text(hjust = 1.0, size = 30)) +
  ggtitle(bquote(.(ranFeatNamePlot) ~ "(" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFfeature$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
}, mc.cores = length(log2ChIPNamesPlot))

## ranLoc
ggObj3_combined_log2ChIP <- mclapply(seq_along(log2ChIPNamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_log2ChIP[[x]][[3]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = quantile)
        ) +
  geom_line(data = summaryDFfeature,
            mapping = aes(colour = quantile),
            size = 1) +
  scale_colour_manual(values = quantileColours) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = quantile),
              alpha = 0.4) +
  scale_fill_manual(values = quantileColours) +
  scale_y_continuous(limits = c(ymin_list_log2ChIP, ymax_list_log2ChIP),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              ((upstream-1000)/binSize)+1,
                              (dim(summaryDFfeature_log2ChIP[[x]][[3]])[1]/quantiles)-((downstream-1000)/binSize),
                              dim(summaryDFfeature_log2ChIP[[x]][[3]])[1]/quantiles),
                   labels = c(paste0("-", "1kb"),
                              featureStartLab,
                              featureEndLab,
                              paste0("+", "1kb"))) +
  geom_vline(xintercept = c(((upstream-1000)/binSize)+1,
                            (dim(summaryDFfeature_log2ChIP[[x]][[3]])[1]/quantiles)-((downstream-1000)/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = log2ChIPNamesPlot[x]) +
  annotation_custom(legendLabs_ranLoc[[1]]) +
  annotation_custom(legendLabs_ranLoc[[2]]) +
  annotation_custom(legendLabs_ranLoc[[3]]) +
  annotation_custom(legendLabs_ranLoc[[4]]) +
  theme_bw() +
  theme(
        axis.ticks = element_line(size = 1.0, colour = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.x = element_text(size = 22, colour = "black"),
        axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
        axis.title = element_text(size = 30, colour = log2ChIPColours[x]),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_rect(size = 3.5, colour = "black"),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
        plot.title = element_text(hjust = 1.0, size = 30)) +
  ggtitle(bquote(.(ranLocNamePlot) ~ "(" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFfeature$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
}, mc.cores = length(log2ChIPNamesPlot))

ggObjGA_combined <- grid.arrange(grobs = c(
                                           ggObj1_combined_log2ChIP,
                                           ggObj2_combined_log2ChIP,
                                           ggObj3_combined_log2ChIP
                                          ),
                                 layout_matrix = cbind(
                                                       1:length(c(log2ChIPNamesPlot)),
                                                       (length(c(log2ChIPNamesPlot))+1):(length(c(log2ChIPNamesPlot))*2),
                                                       ((length(c(log2ChIPNamesPlot))*2)+1):(length(c(log2ChIPNamesPlot))*3)
                                                      ))
ggsave(paste0(plotDir,
              "log2ChIPcontrol_avgProfiles_around_", quantiles, "quantiles",
               "_by_", orderingFactor,
               "_of_CEN180_in_t2t-col.20210610_",
               paste0(chrName, collapse = "_"), "_", align, ".pdf"),
       plot = ggObjGA_combined,
       height = 6.5*length(c(log2ChIPNamesPlot)), width = 21, limitsize = FALSE)

#### Free up memory by removing no longer required objects
rm(
   log2ChIP_featureMats, log2ChIP_ranLocMats,
   log2ChIP_mats_quantiles,
   wideDFfeature_list_log2ChIP,
   tidyDFfeature_list_log2ChIP,
   summaryDFfeature_list_log2ChIP,
   summaryDFfeature_log2ChIP
  ) 
gc()
#####


# ChIP
# Add column names
for(x in seq_along(ChIP_featureMats)) {
  colnames(ChIP_featureMats[[x]]) <- c(paste0("u", 1:((upstream-1000)/binSize)),
                                       paste0("t", (((upstream-1000)/binSize)+1):(((upstream-1000)+bodyLength)/binSize)),
                                       paste0("d", ((((upstream-1000)+bodyLength)/binSize)+1):((((upstream-1000)+bodyLength)/binSize)+((downstream-1000)/binSize))))
  colnames(ChIP_ranLocMats[[x]]) <- c(paste0("u", 1:((upstream-1000)/binSize)),
                                      paste0("t", (((upstream-1000)/binSize)+1):(((upstream-1000)+bodyLength)/binSize)),
                                      paste0("d", ((((upstream-1000)+bodyLength)/binSize)+1):((((upstream-1000)+bodyLength)/binSize)+((downstream-1000)/binSize))))
}

# Subdivide coverage matrices into above-defined quantiles and random groupings
ChIP_mats_quantiles <- mclapply(seq_along(ChIP_featureMats), function(x) {
  list(
       # feature quantiles
       lapply(1:quantiles, function(k) {
         ChIP_featureMats[[x]][quantileIndices[[k]],]
       }),
       # feature random groupings
       lapply(1:quantiles, function(k) {
         ChIP_featureMats[[x]][randomPCIndices[[k]],]
       }),
       # random loci groupings
       lapply(1:quantiles, function(k) {
         ChIP_ranLocMats[[x]][quantileIndices[[k]],]
       })
      ) 
}, mc.cores = length(ChIP_featureMats))

# Transpose matrix and convert into dataframe
# in which first column is window name
wideDFfeature_list_ChIP <- mclapply(seq_along(ChIP_mats_quantiles), function(x) {
  lapply(seq_along(ChIP_mats_quantiles[[x]]), function(y) {
    lapply(seq_along(ChIP_mats_quantiles[[x]][[y]]), function(k) {
      data.frame(window = colnames(ChIP_mats_quantiles[[x]][[y]][[k]]),
                 t(ChIP_mats_quantiles[[x]][[y]][[k]]))
    })
  })
}, mc.cores = length(ChIP_mats_quantiles))

# Convert into tidy data.frame (long format)
tidyDFfeature_list_ChIP  <- mclapply(seq_along(wideDFfeature_list_ChIP), function(x) {
  lapply(seq_along(ChIP_mats_quantiles[[x]]), function(y) {
    lapply(seq_along(ChIP_mats_quantiles[[x]][[y]]), function(k) {
      gather(data  = wideDFfeature_list_ChIP[[x]][[y]][[k]],
             key   = feature,
             value = coverage,
             -window)
    })
  }) 
}, mc.cores = length(wideDFfeature_list_ChIP))

# Order levels of factor "window" so that sequential levels
# correspond to sequential windows
for(x in seq_along(tidyDFfeature_list_ChIP)) {
  for(y in seq_along(ChIP_mats_quantiles[[x]])) {
    for(k in seq_along(ChIP_mats_quantiles[[x]][[y]])) {
      tidyDFfeature_list_ChIP[[x]][[y]][[k]]$window <- factor(tidyDFfeature_list_ChIP[[x]][[y]][[k]]$window,
                                                                  levels = as.character(wideDFfeature_list_ChIP[[x]][[y]][[k]]$window))
    }
  }
}

# Create summary data.frame in which each row corresponds to a window (Column 1),
# Column2 is the number of coverage values (features) per window,
# Column3 is the mean of coverage values per window,
# Column4 is the standard deviation of coverage values per window,
# Column5 is the standard error of the mean of coverage values per window,
# Column6 is the lower bound of the 95% confidence interval, and
# Column7 is the upper bound of the 95% confidence interval
summaryDFfeature_list_ChIP  <- mclapply(seq_along(tidyDFfeature_list_ChIP), function(x) {
  lapply(seq_along(ChIP_mats_quantiles[[x]]), function(y) {
    lapply(seq_along(ChIP_mats_quantiles[[x]][[y]]), function(k) {
      data.frame(window = as.character(wideDFfeature_list_ChIP[[x]][[y]][[k]]$window),
                 n      = tapply(X     = tidyDFfeature_list_ChIP[[x]][[y]][[k]]$coverage,
                                 INDEX = tidyDFfeature_list_ChIP[[x]][[y]][[k]]$window,
                                 FUN   = length),
                 mean   = tapply(X     = tidyDFfeature_list_ChIP[[x]][[y]][[k]]$coverage,
                                 INDEX = tidyDFfeature_list_ChIP[[x]][[y]][[k]]$window,
                                 FUN   = mean,
                                 na.rm = TRUE),
                 sd     = tapply(X     = tidyDFfeature_list_ChIP[[x]][[y]][[k]]$coverage,
                                 INDEX = tidyDFfeature_list_ChIP[[x]][[y]][[k]]$window,
                                 FUN   = sd,
                                 na.rm = TRUE))
    })
  })
}, mc.cores = length(tidyDFfeature_list_ChIP))

for(x in seq_along(summaryDFfeature_list_ChIP)) {
  for(y in seq_along(ChIP_mats_quantiles[[x]])) {
    for(k in seq_along(ChIP_mats_quantiles[[x]][[y]])) {
      summaryDFfeature_list_ChIP[[x]][[y]][[k]]$window <- factor(summaryDFfeature_list_ChIP[[x]][[y]][[k]]$window,
                                                                     levels = as.character(wideDFfeature_list_ChIP[[x]][[y]][[k]]$window))
      summaryDFfeature_list_ChIP[[x]][[y]][[k]]$winNo <- factor(1:dim(summaryDFfeature_list_ChIP[[x]][[y]][[k]])[1])
      summaryDFfeature_list_ChIP[[x]][[y]][[k]]$sem <- summaryDFfeature_list_ChIP[[x]][[y]][[k]]$sd/sqrt(summaryDFfeature_list_ChIP[[x]][[y]][[k]]$n-1)
      summaryDFfeature_list_ChIP[[x]][[y]][[k]]$CI_lower <- summaryDFfeature_list_ChIP[[x]][[y]][[k]]$mean -
        qt(0.975, df = summaryDFfeature_list_ChIP[[x]][[y]][[k]]$n-1)*summaryDFfeature_list_ChIP[[x]][[y]][[k]]$sem
      summaryDFfeature_list_ChIP[[x]][[y]][[k]]$CI_upper <- summaryDFfeature_list_ChIP[[x]][[y]][[k]]$mean +
        qt(0.975, df = summaryDFfeature_list_ChIP[[x]][[y]][[k]]$n-1)*summaryDFfeature_list_ChIP[[x]][[y]][[k]]$sem
    }
  }
}

quantileNames <- paste0(rep("Quantile ", quantiles), 1:quantiles)
randomPCNames <- paste0(rep("Random ", quantiles), 1:quantiles)
for(x in seq_along(summaryDFfeature_list_ChIP)) {
  # feature quantiles
  names(summaryDFfeature_list_ChIP[[x]][[1]]) <- quantileNames
  # feature random groupings
  names(summaryDFfeature_list_ChIP[[x]][[2]]) <- randomPCNames
  # random loci groupings
  names(summaryDFfeature_list_ChIP[[x]][[3]]) <- randomPCNames
}

# Convert list of lists of lists of feature quantiles summaryDFfeature_list_ChIP into
# a list of lists of single data.frames containing all feature quantiles for plotting
summaryDFfeature_ChIP  <- mclapply(seq_along(summaryDFfeature_list_ChIP), function(x) {
  lapply(seq_along(ChIP_mats_quantiles[[x]]), function(y) {
    bind_rows(summaryDFfeature_list_ChIP[[x]][[y]], .id = "quantile")
  })
}, mc.cores = length(summaryDFfeature_list_ChIP))
for(x in seq_along(summaryDFfeature_ChIP)) {
  # feature quantiles
  summaryDFfeature_ChIP[[x]][[1]]$quantile <- factor(summaryDFfeature_ChIP[[x]][[1]]$quantile,
                                                         levels = names(summaryDFfeature_list_ChIP[[x]][[1]]))
  # feature random groupings
  summaryDFfeature_ChIP[[x]][[2]]$quantile <- factor(summaryDFfeature_ChIP[[x]][[2]]$quantile,
                                                         levels = names(summaryDFfeature_list_ChIP[[x]][[2]]))
  # random loci groupings
  summaryDFfeature_ChIP[[x]][[3]]$quantile <- factor(summaryDFfeature_ChIP[[x]][[3]]$quantile,
                                                         levels = names(summaryDFfeature_list_ChIP[[x]][[3]]))
}

# Define y-axis limits
ymin_list_ChIP <- min(unlist(lapply(seq_along(summaryDFfeature_ChIP), function(x) {
  min(c(summaryDFfeature_ChIP[[x]][[1]]$CI_lower,
        summaryDFfeature_ChIP[[x]][[2]]$CI_lower,
        summaryDFfeature_ChIP[[x]][[3]]$CI_lower))
})))
ymax_list_ChIP <- max(unlist(lapply(seq_along(summaryDFfeature_ChIP), function(x) {
  max(c(summaryDFfeature_ChIP[[x]][[1]]$CI_upper,
        summaryDFfeature_ChIP[[x]][[2]]$CI_upper,
        summaryDFfeature_ChIP[[x]][[3]]$CI_upper))
})))

# Define legend labels
legendLabs_feature <- lapply(seq_along(quantileNames), function(x) {
  grobTree(textGrob(bquote(.(quantileNames[x])),
                    x = legendPos[1], y = legendPos[2]-((x-1)*0.06), just = "left",
                    gp = gpar(col = quantileColours[x], fontsize = 18)))
})
legendLabs_ranFeat <- lapply(seq_along(randomPCNames), function(x) {
  grobTree(textGrob(bquote(.(randomPCNames[x])),
                    x = legendPos[1], y = legendPos[2]-((x-1)*0.06), just = "left",
                    gp = gpar(col = quantileColours[x], fontsize = 18)))
})
legendLabs_ranLoc <- lapply(seq_along(randomPCNames), function(x) {
  grobTree(textGrob(bquote(.(randomPCNames[x])),
                    x = legendPos[1], y = legendPos[2]-((x-1)*0.06), just = "left",
                    gp = gpar(col = quantileColours[x], fontsize = 18)))
})

# Plot average profiles with 95% CI ribbon
## feature
ggObj1_combined_ChIP <- mclapply(seq_along(ChIPNamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_ChIP[[x]][[1]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = quantile)
        ) +
  geom_line(data = summaryDFfeature,
            mapping = aes(colour = quantile),
            size = 1) +
  scale_colour_manual(values = quantileColours) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = quantile),
              alpha = 0.4) +
  scale_fill_manual(values = quantileColours) +
  scale_y_continuous(limits = c(ymin_list_ChIP, ymax_list_ChIP),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              ((upstream-1000)/binSize)+1,
                              (dim(summaryDFfeature_ChIP[[x]][[1]])[1]/quantiles)-((downstream-1000)/binSize),
                              dim(summaryDFfeature_ChIP[[x]][[1]])[1]/quantiles),
                   labels = c(paste0("-", "1kb"),
                              featureStartLab,
                              featureEndLab,
                              paste0("+", "1kb"))) +
  geom_vline(xintercept = c(((upstream-1000)/binSize)+1,
                            (dim(summaryDFfeature_ChIP[[x]][[1]])[1]/quantiles)-((downstream-1000)/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = ChIPNamesPlot[x]) +
  annotation_custom(legendLabs_feature[[1]]) +
  annotation_custom(legendLabs_feature[[2]]) +
  annotation_custom(legendLabs_feature[[3]]) +
  annotation_custom(legendLabs_feature[[4]]) +
  theme_bw() +
  theme(
        axis.ticks = element_line(size = 1.0, colour = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.x = element_text(size = 22, colour = "black"),
        axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
        axis.title = element_text(size = 30, colour = ChIPColours[x]),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_rect(size = 3.5, colour = "black"),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
        plot.title = element_text(hjust = 1.0, size = 30)) +
  ggtitle(bquote(.(featureNamePlot) ~ "(" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFfeature$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
}, mc.cores = length(ChIPNamesPlot))

## ranFeat
ggObj2_combined_ChIP <- mclapply(seq_along(ChIPNamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_ChIP[[x]][[2]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = quantile)
        ) +
  geom_line(data = summaryDFfeature,
            mapping = aes(colour = quantile),
            size = 1) +
  scale_colour_manual(values = quantileColours) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = quantile),
              alpha = 0.4) +
  scale_fill_manual(values = quantileColours) +
  scale_y_continuous(limits = c(ymin_list_ChIP, ymax_list_ChIP),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              ((upstream-1000)/binSize)+1,
                              (dim(summaryDFfeature_ChIP[[x]][[2]])[1]/quantiles)-((downstream-1000)/binSize),
                              dim(summaryDFfeature_ChIP[[x]][[2]])[1]/quantiles),
                   labels = c(paste0("-", "1kb"),
                              featureStartLab,
                              featureEndLab,
                              paste0("+", "1kb"))) +
  geom_vline(xintercept = c(((upstream-1000)/binSize)+1,
                            (dim(summaryDFfeature_ChIP[[x]][[2]])[1]/quantiles)-((downstream-1000)/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = ChIPNamesPlot[x]) +
  annotation_custom(legendLabs_ranFeat[[1]]) +
  annotation_custom(legendLabs_ranFeat[[2]]) +
  annotation_custom(legendLabs_ranFeat[[3]]) +
  annotation_custom(legendLabs_ranFeat[[4]]) +
  theme_bw() +
  theme(
        axis.ticks = element_line(size = 1.0, colour = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.x = element_text(size = 22, colour = "black"),
        axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
        axis.title = element_text(size = 30, colour = ChIPColours[x]),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_rect(size = 3.5, colour = "black"),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
        plot.title = element_text(hjust = 1.0, size = 30)) +
  ggtitle(bquote(.(ranFeatNamePlot) ~ "(" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFfeature$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
}, mc.cores = length(ChIPNamesPlot))

## ranLoc
ggObj3_combined_ChIP <- mclapply(seq_along(ChIPNamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_ChIP[[x]][[3]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = quantile)
        ) +
  geom_line(data = summaryDFfeature,
            mapping = aes(colour = quantile),
            size = 1) +
  scale_colour_manual(values = quantileColours) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = quantile),
              alpha = 0.4) +
  scale_fill_manual(values = quantileColours) +
  scale_y_continuous(limits = c(ymin_list_ChIP, ymax_list_ChIP),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              ((upstream-1000)/binSize)+1,
                              (dim(summaryDFfeature_ChIP[[x]][[3]])[1]/quantiles)-((downstream-1000)/binSize),
                              dim(summaryDFfeature_ChIP[[x]][[3]])[1]/quantiles),
                   labels = c(paste0("-", "1kb"),
                              featureStartLab,
                              featureEndLab,
                              paste0("+", "1kb"))) +
  geom_vline(xintercept = c(((upstream-1000)/binSize)+1,
                            (dim(summaryDFfeature_ChIP[[x]][[3]])[1]/quantiles)-((downstream-1000)/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = ChIPNamesPlot[x]) +
  annotation_custom(legendLabs_ranLoc[[1]]) +
  annotation_custom(legendLabs_ranLoc[[2]]) +
  annotation_custom(legendLabs_ranLoc[[3]]) +
  annotation_custom(legendLabs_ranLoc[[4]]) +
  theme_bw() +
  theme(
        axis.ticks = element_line(size = 1.0, colour = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.x = element_text(size = 22, colour = "black"),
        axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
        axis.title = element_text(size = 30, colour = ChIPColours[x]),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_rect(size = 3.5, colour = "black"),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
        plot.title = element_text(hjust = 1.0, size = 30)) +
  ggtitle(bquote(.(ranLocNamePlot) ~ "(" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFfeature$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
}, mc.cores = length(ChIPNamesPlot))

ggObjGA_combined <- grid.arrange(grobs = c(
                                           ggObj1_combined_ChIP,
                                           ggObj2_combined_ChIP,
                                           ggObj3_combined_ChIP
                                          ),
                                 layout_matrix = cbind(
                                                       1:length(c(ChIPNamesPlot)),
                                                       (length(c(ChIPNamesPlot))+1):(length(c(ChIPNamesPlot))*2),
                                                       ((length(c(ChIPNamesPlot))*2)+1):(length(c(ChIPNamesPlot))*3)
                                                      ))
ggsave(paste0(plotDir,
              "ChIP_avgProfiles_around_", quantiles, "quantiles",
               "_by_", orderingFactor,
               "_of_CEN180_in_t2t-col.20210610_",
               paste0(chrName, collapse = "_"), "_", align, ".pdf"),
       plot = ggObjGA_combined,
       height = 6.5*length(c(ChIPNamesPlot)), width = 21, limitsize = FALSE)

#### Free up memory by removing no longer required objects
rm(
   ChIP_featureMats, ChIP_ranLocMats,
   ChIP_mats_quantiles,
   wideDFfeature_list_ChIP,
   tidyDFfeature_list_ChIP,
   summaryDFfeature_list_ChIP,
   summaryDFfeature_ChIP
  ) 
gc()
#####


# control
# Add column names
for(x in seq_along(control_featureMats)) {
  colnames(control_featureMats[[x]]) <- c(paste0("u", 1:((upstream-1000)/binSize)),
                                          paste0("t", (((upstream-1000)/binSize)+1):(((upstream-1000)+bodyLength)/binSize)),
                                          paste0("d", ((((upstream-1000)+bodyLength)/binSize)+1):((((upstream-1000)+bodyLength)/binSize)+((downstream-1000)/binSize))))
  colnames(control_ranLocMats[[x]]) <- c(paste0("u", 1:((upstream-1000)/binSize)),
                                         paste0("t", (((upstream-1000)/binSize)+1):(((upstream-1000)+bodyLength)/binSize)),
                                         paste0("d", ((((upstream-1000)+bodyLength)/binSize)+1):((((upstream-1000)+bodyLength)/binSize)+((downstream-1000)/binSize))))
}

# Subdivide coverage matrices into above-defined quantiles and random groupings
control_mats_quantiles <- mclapply(seq_along(control_featureMats), function(x) {
  list(
       # feature quantiles
       lapply(1:quantiles, function(k) {
         control_featureMats[[x]][quantileIndices[[k]],]
       }),
       # feature random groupings
       lapply(1:quantiles, function(k) {
         control_featureMats[[x]][randomPCIndices[[k]],]
       }),
       # random loci groupings
       lapply(1:quantiles, function(k) {
         control_ranLocMats[[x]][quantileIndices[[k]],]
       })
      ) 
}, mc.cores = length(control_featureMats))

# Transpose matrix and convert into dataframe
# in which first column is window name
wideDFfeature_list_control <- mclapply(seq_along(control_mats_quantiles), function(x) {
  lapply(seq_along(control_mats_quantiles[[x]]), function(y) {
    lapply(seq_along(control_mats_quantiles[[x]][[y]]), function(k) {
      data.frame(window = colnames(control_mats_quantiles[[x]][[y]][[k]]),
                 t(control_mats_quantiles[[x]][[y]][[k]]))
    })
  })
}, mc.cores = length(control_mats_quantiles))

# Convert into tidy data.frame (long format)
tidyDFfeature_list_control  <- mclapply(seq_along(wideDFfeature_list_control), function(x) {
  lapply(seq_along(control_mats_quantiles[[x]]), function(y) {
    lapply(seq_along(control_mats_quantiles[[x]][[y]]), function(k) {
      gather(data  = wideDFfeature_list_control[[x]][[y]][[k]],
             key   = feature,
             value = coverage,
             -window)
    })
  }) 
}, mc.cores = length(wideDFfeature_list_control))

# Order levels of factor "window" so that sequential levels
# correspond to sequential windows
for(x in seq_along(tidyDFfeature_list_control)) {
  for(y in seq_along(control_mats_quantiles[[x]])) {
    for(k in seq_along(control_mats_quantiles[[x]][[y]])) {
      tidyDFfeature_list_control[[x]][[y]][[k]]$window <- factor(tidyDFfeature_list_control[[x]][[y]][[k]]$window,
                                                                  levels = as.character(wideDFfeature_list_control[[x]][[y]][[k]]$window))
    }
  }
}

# Create summary data.frame in which each row corresponds to a window (Column 1),
# Column2 is the number of coverage values (features) per window,
# Column3 is the mean of coverage values per window,
# Column4 is the standard deviation of coverage values per window,
# Column5 is the standard error of the mean of coverage values per window,
# Column6 is the lower bound of the 95% confidence interval, and
# Column7 is the upper bound of the 95% confidence interval
summaryDFfeature_list_control  <- mclapply(seq_along(tidyDFfeature_list_control), function(x) {
  lapply(seq_along(control_mats_quantiles[[x]]), function(y) {
    lapply(seq_along(control_mats_quantiles[[x]][[y]]), function(k) {
      data.frame(window = as.character(wideDFfeature_list_control[[x]][[y]][[k]]$window),
                 n      = tapply(X     = tidyDFfeature_list_control[[x]][[y]][[k]]$coverage,
                                 INDEX = tidyDFfeature_list_control[[x]][[y]][[k]]$window,
                                 FUN   = length),
                 mean   = tapply(X     = tidyDFfeature_list_control[[x]][[y]][[k]]$coverage,
                                 INDEX = tidyDFfeature_list_control[[x]][[y]][[k]]$window,
                                 FUN   = mean,
                                 na.rm = TRUE),
                 sd     = tapply(X     = tidyDFfeature_list_control[[x]][[y]][[k]]$coverage,
                                 INDEX = tidyDFfeature_list_control[[x]][[y]][[k]]$window,
                                 FUN   = sd,
                                 na.rm = TRUE))
    })
  })
}, mc.cores = length(tidyDFfeature_list_control))

for(x in seq_along(summaryDFfeature_list_control)) {
  for(y in seq_along(control_mats_quantiles[[x]])) {
    for(k in seq_along(control_mats_quantiles[[x]][[y]])) {
      summaryDFfeature_list_control[[x]][[y]][[k]]$window <- factor(summaryDFfeature_list_control[[x]][[y]][[k]]$window,
                                                                     levels = as.character(wideDFfeature_list_control[[x]][[y]][[k]]$window))
      summaryDFfeature_list_control[[x]][[y]][[k]]$winNo <- factor(1:dim(summaryDFfeature_list_control[[x]][[y]][[k]])[1])
      summaryDFfeature_list_control[[x]][[y]][[k]]$sem <- summaryDFfeature_list_control[[x]][[y]][[k]]$sd/sqrt(summaryDFfeature_list_control[[x]][[y]][[k]]$n-1)
      summaryDFfeature_list_control[[x]][[y]][[k]]$CI_lower <- summaryDFfeature_list_control[[x]][[y]][[k]]$mean -
        qt(0.975, df = summaryDFfeature_list_control[[x]][[y]][[k]]$n-1)*summaryDFfeature_list_control[[x]][[y]][[k]]$sem
      summaryDFfeature_list_control[[x]][[y]][[k]]$CI_upper <- summaryDFfeature_list_control[[x]][[y]][[k]]$mean +
        qt(0.975, df = summaryDFfeature_list_control[[x]][[y]][[k]]$n-1)*summaryDFfeature_list_control[[x]][[y]][[k]]$sem
    }
  }
}

quantileNames <- paste0(rep("Quantile ", quantiles), 1:quantiles)
randomPCNames <- paste0(rep("Random ", quantiles), 1:quantiles)
for(x in seq_along(summaryDFfeature_list_control)) {
  # feature quantiles
  names(summaryDFfeature_list_control[[x]][[1]]) <- quantileNames
  # feature random groupings
  names(summaryDFfeature_list_control[[x]][[2]]) <- randomPCNames
  # random loci groupings
  names(summaryDFfeature_list_control[[x]][[3]]) <- randomPCNames
}

# Convert list of lists of lists of feature quantiles summaryDFfeature_list_control into
# a list of lists of single data.frames containing all feature quantiles for plotting
summaryDFfeature_control  <- mclapply(seq_along(summaryDFfeature_list_control), function(x) {
  lapply(seq_along(control_mats_quantiles[[x]]), function(y) {
    bind_rows(summaryDFfeature_list_control[[x]][[y]], .id = "quantile")
  })
}, mc.cores = length(summaryDFfeature_list_control))
for(x in seq_along(summaryDFfeature_control)) {
  # feature quantiles
  summaryDFfeature_control[[x]][[1]]$quantile <- factor(summaryDFfeature_control[[x]][[1]]$quantile,
                                                         levels = names(summaryDFfeature_list_control[[x]][[1]]))
  # feature random groupings
  summaryDFfeature_control[[x]][[2]]$quantile <- factor(summaryDFfeature_control[[x]][[2]]$quantile,
                                                         levels = names(summaryDFfeature_list_control[[x]][[2]]))
  # random loci groupings
  summaryDFfeature_control[[x]][[3]]$quantile <- factor(summaryDFfeature_control[[x]][[3]]$quantile,
                                                         levels = names(summaryDFfeature_list_control[[x]][[3]]))
}

# Define y-axis limits
ymin_list_control <- min(unlist(lapply(seq_along(summaryDFfeature_control), function(x) {
  min(c(summaryDFfeature_control[[x]][[1]]$CI_lower,
        summaryDFfeature_control[[x]][[2]]$CI_lower,
        summaryDFfeature_control[[x]][[3]]$CI_lower))
})))
ymax_list_control <- max(unlist(lapply(seq_along(summaryDFfeature_control), function(x) {
  max(c(summaryDFfeature_control[[x]][[1]]$CI_upper,
        summaryDFfeature_control[[x]][[2]]$CI_upper,
        summaryDFfeature_control[[x]][[3]]$CI_upper))
})))

# Define legend labels
legendLabs_feature <- lapply(seq_along(quantileNames), function(x) {
  grobTree(textGrob(bquote(.(quantileNames[x])),
                    x = legendPos[1], y = legendPos[2]-((x-1)*0.06), just = "left",
                    gp = gpar(col = quantileColours[x], fontsize = 18)))
})
legendLabs_ranFeat <- lapply(seq_along(randomPCNames), function(x) {
  grobTree(textGrob(bquote(.(randomPCNames[x])),
                    x = legendPos[1], y = legendPos[2]-((x-1)*0.06), just = "left",
                    gp = gpar(col = quantileColours[x], fontsize = 18)))
})
legendLabs_ranLoc <- lapply(seq_along(randomPCNames), function(x) {
  grobTree(textGrob(bquote(.(randomPCNames[x])),
                    x = legendPos[1], y = legendPos[2]-((x-1)*0.06), just = "left",
                    gp = gpar(col = quantileColours[x], fontsize = 18)))
})

# Plot average profiles with 95% CI ribbon
## feature
ggObj1_combined_control <- mclapply(seq_along(controlNamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_control[[x]][[1]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = quantile)
        ) +
  geom_line(data = summaryDFfeature,
            mapping = aes(colour = quantile),
            size = 1) +
  scale_colour_manual(values = quantileColours) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = quantile),
              alpha = 0.4) +
  scale_fill_manual(values = quantileColours) +
  scale_y_continuous(limits = c(ymin_list_control, ymax_list_control),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              ((upstream-1000)/binSize)+1,
                              (dim(summaryDFfeature_control[[x]][[1]])[1]/quantiles)-((downstream-1000)/binSize),
                              dim(summaryDFfeature_control[[x]][[1]])[1]/quantiles),
                   labels = c(paste0("-", "1kb"),
                              featureStartLab,
                              featureEndLab,
                              paste0("+", "1kb"))) +
  geom_vline(xintercept = c(((upstream-1000)/binSize)+1,
                            (dim(summaryDFfeature_control[[x]][[1]])[1]/quantiles)-((downstream-1000)/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = controlNamesPlot[x]) +
  annotation_custom(legendLabs_feature[[1]]) +
  annotation_custom(legendLabs_feature[[2]]) +
  annotation_custom(legendLabs_feature[[3]]) +
  annotation_custom(legendLabs_feature[[4]]) +
  theme_bw() +
  theme(
        axis.ticks = element_line(size = 1.0, colour = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.x = element_text(size = 22, colour = "black"),
        axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
        axis.title = element_text(size = 30, colour = controlColours[x]),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_rect(size = 3.5, colour = "black"),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
        plot.title = element_text(hjust = 1.0, size = 30)) +
  ggtitle(bquote(.(featureNamePlot) ~ "(" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFfeature$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
}, mc.cores = length(controlNamesPlot))

## ranFeat
ggObj2_combined_control <- mclapply(seq_along(controlNamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_control[[x]][[2]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = quantile)
        ) +
  geom_line(data = summaryDFfeature,
            mapping = aes(colour = quantile),
            size = 1) +
  scale_colour_manual(values = quantileColours) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = quantile),
              alpha = 0.4) +
  scale_fill_manual(values = quantileColours) +
  scale_y_continuous(limits = c(ymin_list_control, ymax_list_control),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              ((upstream-1000)/binSize)+1,
                              (dim(summaryDFfeature_control[[x]][[2]])[1]/quantiles)-((downstream-1000)/binSize),
                              dim(summaryDFfeature_control[[x]][[2]])[1]/quantiles),
                   labels = c(paste0("-", "1kb"),
                              featureStartLab,
                              featureEndLab,
                              paste0("+", "1kb"))) +
  geom_vline(xintercept = c(((upstream-1000)/binSize)+1,
                            (dim(summaryDFfeature_control[[x]][[2]])[1]/quantiles)-((downstream-1000)/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = controlNamesPlot[x]) +
  annotation_custom(legendLabs_ranFeat[[1]]) +
  annotation_custom(legendLabs_ranFeat[[2]]) +
  annotation_custom(legendLabs_ranFeat[[3]]) +
  annotation_custom(legendLabs_ranFeat[[4]]) +
  theme_bw() +
  theme(
        axis.ticks = element_line(size = 1.0, colour = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.x = element_text(size = 22, colour = "black"),
        axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
        axis.title = element_text(size = 30, colour = controlColours[x]),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_rect(size = 3.5, colour = "black"),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
        plot.title = element_text(hjust = 1.0, size = 30)) +
  ggtitle(bquote(.(ranFeatNamePlot) ~ "(" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFfeature$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
}, mc.cores = length(controlNamesPlot))

## ranLoc
ggObj3_combined_control <- mclapply(seq_along(controlNamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_control[[x]][[3]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = quantile)
        ) +
  geom_line(data = summaryDFfeature,
            mapping = aes(colour = quantile),
            size = 1) +
  scale_colour_manual(values = quantileColours) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = quantile),
              alpha = 0.4) +
  scale_fill_manual(values = quantileColours) +
  scale_y_continuous(limits = c(ymin_list_control, ymax_list_control),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              ((upstream-1000)/binSize)+1,
                              (dim(summaryDFfeature_control[[x]][[3]])[1]/quantiles)-((downstream-1000)/binSize),
                              dim(summaryDFfeature_control[[x]][[3]])[1]/quantiles),
                   labels = c(paste0("-", "1kb"),
                              featureStartLab,
                              featureEndLab,
                              paste0("+", "1kb"))) +
  geom_vline(xintercept = c(((upstream-1000)/binSize)+1,
                            (dim(summaryDFfeature_control[[x]][[3]])[1]/quantiles)-((downstream-1000)/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = controlNamesPlot[x]) +
  annotation_custom(legendLabs_ranLoc[[1]]) +
  annotation_custom(legendLabs_ranLoc[[2]]) +
  annotation_custom(legendLabs_ranLoc[[3]]) +
  annotation_custom(legendLabs_ranLoc[[4]]) +
  theme_bw() +
  theme(
        axis.ticks = element_line(size = 1.0, colour = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.x = element_text(size = 22, colour = "black"),
        axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
        axis.title = element_text(size = 30, colour = controlColours[x]),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_rect(size = 3.5, colour = "black"),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
        plot.title = element_text(hjust = 1.0, size = 30)) +
  ggtitle(bquote(.(ranLocNamePlot) ~ "(" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFfeature$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
}, mc.cores = length(controlNamesPlot))

ggObjGA_combined <- grid.arrange(grobs = c(
                                           ggObj1_combined_control,
                                           ggObj2_combined_control,
                                           ggObj3_combined_control
                                          ),
                                 layout_matrix = cbind(
                                                       1:length(c(controlNamesPlot)),
                                                       (length(c(controlNamesPlot))+1):(length(c(controlNamesPlot))*2),
                                                       ((length(c(controlNamesPlot))*2)+1):(length(c(controlNamesPlot))*3)
                                                      ))
ggsave(paste0(plotDir,
              "control_avgProfiles_around_", quantiles, "quantiles",
               "_by_", orderingFactor,
               "_of_CEN180_in_t2t-col.20210610_",
               paste0(chrName, collapse = "_"), "_", align, ".pdf"),
       plot = ggObjGA_combined,
       height = 6.5*length(c(controlNamesPlot)), width = 21, limitsize = FALSE)

#### Free up memory by removing no longer required objects
rm(
   control_featureMats, control_ranLocMats,
   control_mats_quantiles,
   wideDFfeature_list_control,
   tidyDFfeature_list_control,
   summaryDFfeature_list_control,
   summaryDFfeature_control
  ) 
gc()
#####

