#!/applications/R/R-3.5.0/bin/Rscript

# author: Andy Tock
# contact: ajt200@cam.ac.uk
# date: 18.11.2020

# Calculate and plot metaprofiles of ChIP-seq, MNase-seq, etc.
# (CEN180 windowed means and 95% confidence intervals, CIs)
# for each group of CEN180s, defined either by
# decreasing orderingFactor or randomly

# Usage:
# /applications/R/R-3.5.0/bin/Rscript CEN180_quantile_metaprofiles.R 'Chr1,Chr2,Chr3,Chr4,Chr5' wSNV 4 both 180 1000 1kb 10 '0.02,0.30'

#chrName <- unlist(strsplit("Chr1,Chr2,Chr3,Chr4,Chr5",
#                           split = ","))
#orderingFactor <- "wSNV"
#quantiles <- 4
#align <- "both"
#bodyLength <- 180
#upstream <- 1000
#downstream <- 1000
#flankName <- "1kb"
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

outDir <- paste0("quantiles_by_", orderingFactor, "/",
                 paste0(chrName, collapse = "_"), "/")
plotDir <- paste0(outDir, "plots/")
system(paste0("[ -d ", outDir, " ] || mkdir -p ", outDir))
system(paste0("[ -d ", plotDir, " ] || mkdir -p ", plotDir))

# Define plot titles
if(grepl("_in_", orderingFactor)) {
  featureNamePlot <- paste0(sub("_in_\\w+", "", orderingFactor), " quantiles")
} else if(grepl("SNV", orderingFactor)) {
  featureNamePlot <- "SNV quantiles"
} else if(orderingFactor == "array_size") {
  featureNamePlot <- "Array-size quantiles"
}
ranFeatNamePlot <- "Random quantiles"
ranLocNamePlot <- "RanLoc quantiles"

# Define quantile colours
quantileColours <- c("red", "purple", "blue", "navy")

# Define feature start and end labels for plotting
featureStartLab <- "Start"
featureEndLab <- "End"

# Load table of features grouped into quantiles
featuresDF <- read.table(paste0(outDir,
                                "features_", quantiles, "quantiles",
                                "_by_", orderingFactor,
                                "_of_CEN180_in_RaGOO_v2.0_",
                                paste0(chrName, collapse = "_"), ".tsv"),
                         header = T, sep = "\t", row.names = NULL)

# Load features to confirm feature (row) ordering in "featuresDF" is the same
# as in "features" (which was used for generating the coverage matrices)
features <- lapply(seq_along(chrName), function(y) {
  read.table(paste0("CEN180_in_RaGOO_v2.0_",
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
colnames(features) <- c("chr", "start", "end", "featureID", "wSNV", "strand")
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

# Define seed so that random selections are reproducible
set.seed(93750174)

# Divide features into random sets of equal number,
# with the same number of CEN180s per chromosome as
# above-defined orderingFactor-defined feature quantiles
randomPCIndices <- lapply(1:quantiles, function(k) {
  randomPCIndicesk <- NULL
  for(i in 1:length(chrName)) {
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
ChIPNames <- c(
               "WT_CENH3_Rep1_ChIP_SRR4430537",
               "LoCENH3_Rep1_ChIP_SRR4430543",
               "ZmCENH3_Rep1_ChIP_SRR4430549",
               "WT_MNase_Rep1",
               "WT_REC8_HA_Rep2_ChIP",
               "WT_ASY1_Rep1_ChIP",
               "WT_MTOPVIB_HA_Rep1_ChIP",
               "WT_MTOPVIB_HA_Rep2_ChIP",
               "WT_DMC1_V5_Rep1_ChIP",
               "WT_SPO11oligos_Rep1",
#               "WT_SPO11oligos_Rep2",
#               "WT_SPO11oligos_Rep3",
               "WT_H3K9me2_Rep1_ChIP",
               "WT_H3K4me1_Rep1_ChIP",
               "WT_H3K4me2_Rep1_ChIP",
               "WT_H3K27me1_Rep1_ChIP",
               "WT_H3K4me3_ChIP14",
               "WT_H3K27me3_ChIP_SRR1509478",
#               "H2A_ChIP_set1",
               "H2AW6_ChIP_set1",
               "H2AX_ChIP_set1",
               "H2AZ_ChIP_set1",
#               "H3_ChIP_set1",
               "H2AW7_ChIP_set2"
#               "H3K27me3_ChIP_set2",
#               "H3K4me1_ChIP_set2",
#               "H3K4me3_ChIP_set2",
#               "H3_ChIP_set2",
#               "H4K20me1_ChIP_set2",
#               "H1_ChIP_set3",
#               "H3_ChIP_set3",
#               "H3K27me1_ChIP_set5",
#               "H3K36me3_ChIP_set5",
#               "H3K9me1_ChIP_set5",
#               "H3K9me2_ChIP_set5",
#               "H3_ChIP_set5",
#               "H3_ChIP_set7",
#               "HTB1_ChIP_set7",
#               "HTB2_ChIP_set7",
#               "HTB3_ChIP_set7",
#               "HTB4_ChIP_set7"
              )
ChIPNamesDir <- c(
                  rep("CENH3_seedlings_Maheshwari_Comai_2017_GenomeRes/snakemake_ChIPseq_RaGOO_v2.0", 3),
                  "150701_Kyuha_MNase/WT/snakemake_ChIPseq_RaGOO_v2.0",
                  "REC8_pooled/snakemake_ChIPseq_RaGOO_v2.0",
                  "20190722_cal66_Athaliana_ChIPseq_ASY1/fastq_pooled/snakemake_ChIPseq_RaGOO_v2.0",
                  rep("20190819_dh580_Athaliana_ChIPseq_MTOPVIB/fastq_pooled/snakemake_ChIPseq_RaGOO_v2.0", 2),
                  "20190917_dh580_Athaliana_ChIPseq_DMC1/fastq_pooled/snakemake_ChIPseq_RaGOO_v2.0",
                  rep("160518_Kyuha_SPO11oligos/WT/snakemake_SPO11oligos_RaGOO_v2.0", 1),
                  "170101_Chris_H3K9me2_ChIP/WT/snakemake_ChIPseq_RaGOO_v2.0",
                  rep("170920_Chris_histone_ChIP/snakemake_ChIPseq_RaGOO_v2.0", 3),
                  "160601_Kyuha_H3K4me3_ChIP/WT/snakemake_ChIPseq_RaGOO_v2.0",
                  "H3K27me3_bud_UWMadison_2015/snakemake_ChIPseq_RaGOO_v2.0",
                  rep("histones_10dps_seedling_Berger_2020/snakemake_ChIPseq_RaGOO_v2.0", 4)
                 )
log2ChIPNamesPlot <- c(
                       "CENH3",
                       "LoCENH3",
                       "ZmCENH3",
                       "MNase",
                       "REC8-HA",
                       "ASY1",
                       "MTOPVIB-HA Rep1",
                       "MTOPVIB-HA Rep2",
                       "DMC1-V5",
                       "SPO11-1-oligos Rep1",
#                       "SPO11-1-oligos Rep2",
#                       "SPO11-1-oligos Rep3",
                       "H3K9me2 (CL)",
                       "H3K4me1 (CL)",
                       "H3K4me2 (CL)",
                       "H3K27me1 (CL)",
                       "H3K4me3 (KC)",
                       "H3K27me3 (UWM)",
#                       "H2A set1",
                       "H2A.W6 set1",
                       "H2A.X set1",
                       "H2A.Z set1",
#                       "H3 set1",
                       "H2A.W7 set2"
#                       "H3K27me3 set2",
#                       "H3K4me1 set2",
#                       "H3K4me3 set2",
#                       "H3 set2",
#                       "H4K20me1 set2",
#                       "H1 set3",
#                       "H3 set3",
#                       "H3K27me1 set5",
#                       "H3K36me3 set5",
#                       "H3K9me1 set5",
#                       "H3K9me2 set5",
#                       "H3 set5",
#                       "H3 set7",
#                       "HTB1 set7",
#                       "HTB2 set7",
#                       "HTB3 set7",
#                       "HTB4 set7"
                      )
log2ChIPColours <- c(
                     rep("black", length(log2ChIPNamesPlot))
                    )
ChIPDirs <- sapply(seq_along(ChIPNames), function(x) {
  paste0("/home/ajt200/analysis/",
         ChIPNamesDir[x],
         "/mapped/CEN180profiles/matrices/")
})

controlNames <- c(
                  "WT_REC8_Myc_Rep1_input",
                  "WT_gDNA_Rep1",
                  "WT_gDNA_Rep1_R1",
                  "input_set1",
                  "input_set2"
#                  "input_set3",
#                  "input_set5",
#                  "input_set7"
                 )
controlNamesDir <- c(
                     "REC8_pooled/snakemake_ChIPseq_RaGOO_v2.0",
                     "150701_Natasha_gDNA/WT/snakemake_ChIPseq_RaGOO_v2.0",
                     "150701_Natasha_gDNA/WT/R1/snakemake_SPO11oligos_RaGOO_v2.0",
                     rep("histones_10dps_seedling_Berger_2020/snakemake_ChIPseq_RaGOO_v2.0", 2)
                    )
controlNamesPlot <- c(
                      "Input",
                      "gDNA",
                      "gDNA",
                      rep("Input", 2)
                     )
controlColours <- c(
                    rep("black", length(controlNamesPlot))
                   )
controlDirs <- sapply(seq_along(controlNames), function(x) {
  paste0("/home/ajt200/analysis/",
         controlNamesDir[x],
         "/mapped/CEN180profiles/matrices/")
})

## ChIP
# feature
ChIP_featureMats <- mclapply(seq_along(ChIPNames), function(x) {
  lapply(seq_along(chrName), function(y) {
    as.matrix(read.table(paste0(ChIPDirs[x],
                                ChIPNames[x],
                                "_MappedOn_Athaliana_ONT_RaGOO_v2.0_lowXM_", align, "_sort_norm_CEN180_in_",
                                chrName[y], "_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                         header = F, skip = 3))
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
    as.matrix(read.table(paste0(controlDirs[x],
                                controlNames[x],
                                "_MappedOn_Athaliana_ONT_RaGOO_v2.0_lowXM_", align, "_sort_norm_CEN180_in_",
                                chrName[y], "_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                         header = F, skip = 3))
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
                                "_MappedOn_Athaliana_ONT_RaGOO_v2.0_lowXM_", align, "_sort_norm_CEN180_in_",
                                chrName[y], "_ranLoc_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                         header = F, skip = 3))
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
    as.matrix(read.table(paste0(controlDirs[x],
                                controlNames[x],
                                "_MappedOn_Athaliana_ONT_RaGOO_v2.0_lowXM_", align, "_sort_norm_CEN180_in_",
                                chrName[y], "_ranLoc_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                         header = F, skip = 3))
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
  if ( grepl("MNase", ChIPNames[x]) ) {
    print(paste0(ChIPNames[x], " library; using ", controlNames[2], " for log2((MNase+1)/(gDNA+1)) calculation"))
    log2((ChIP_featureMats[[x]]+1)/(control_featureMats[[2]]+1))
  } else if ( grepl("SPO11oligos", ChIPNames[x]) ) {
    print(paste0(ChIPNames[x], " library; using ", controlNames[3], " for log2((SPO11-1-oligos+1)/(gDNA+1)) calculation"))
    log2((ChIP_featureMats[[x]]+1)/(control_featureMats[[3]]+1))
  } else if ( grepl("set1", ChIPNames[x]) ) {
    print(paste0(ChIPNames[x], " library; using ", controlNames[4], " for log2((ChIP+1)/(input+1)) calculation"))
    log2((ChIP_featureMats[[x]]+1)/(control_featureMats[[4]]+1))
  } else if ( grepl("set2", ChIPNames[x]) ) {
    print(paste0(ChIPNames[x], " library; using ", controlNames[5], " for log2((ChIP+1)/(input+1)) calculation"))
    log2((ChIP_featureMats[[x]]+1)/(control_featureMats[[5]]+1))
  } else if ( grepl("set3", ChIPNames[x]) ) {
    print(paste0(ChIPNames[x], " library; using ", controlNames[6], " for log2((ChIP+1)/(input+1)) calculation"))
    log2((ChIP_featureMats[[x]]+1)/(control_featureMats[[6]]+1))
  } else if ( grepl("set5", ChIPNames[x]) ) {
    print(paste0(ChIPNames[x], " library; using ", controlNames[7], " for log2((ChIP+1)/(input+1)) calculation"))
    log2((ChIP_featureMats[[x]]+1)/(control_featureMats[[7]]+1))
  } else if ( grepl("set7", ChIPNames[x]) ) {
    print(paste0(ChIPNames[x], " library; using ", controlNames[8], " for log2((ChIP+1)/(input+1)) calculation"))
    log2((ChIP_featureMats[[x]]+1)/(control_featureMats[[8]]+1))
  } else {
    print(paste0(ChIPNames[x], " library; using ", controlNames[1], " for log2((ChIP+1)/(input+1)) calculation"))
    log2((ChIP_featureMats[[x]]+1)/(control_featureMats[[1]]+1))
  }
}, mc.cores = length(ChIP_featureMats))

# Conditionally calculate log2(ChIP/control)
# for each matrix depending on library
# ranLoc
log2ChIP_ranLocMats <- mclapply(seq_along(ChIP_ranLocMats), function(x) {
  if ( grepl("MNase", ChIPNames[x]) ) {
    print(paste0(ChIPNames[x], " library; using ", controlNames[2], " for log2((MNase+1)/(gDNA+1)) calculation"))
    log2((ChIP_ranLocMats[[x]]+1)/(control_ranLocMats[[2]]+1))
  } else if ( grepl("SPO11oligos", ChIPNames[x]) ) {
    print(paste0(ChIPNames[x], " library; using ", controlNames[3], " for log2((SPO11-1-oligos+1)/(gDNA+1)) calculation"))
    log2((ChIP_ranLocMats[[x]]+1)/(control_ranLocMats[[3]]+1))
  } else if ( grepl("set1", ChIPNames[x]) ) {
    print(paste0(ChIPNames[x], " library; using ", controlNames[4], " for log2((ChIP+1)/(input+1)) calculation"))
    log2((ChIP_ranLocMats[[x]]+1)/(control_ranLocMats[[4]]+1))
  } else if ( grepl("set2", ChIPNames[x]) ) {
    print(paste0(ChIPNames[x], " library; using ", controlNames[5], " for log2((ChIP+1)/(input+1)) calculation"))
    log2((ChIP_ranLocMats[[x]]+1)/(control_ranLocMats[[5]]+1))
  } else if ( grepl("set3", ChIPNames[x]) ) {
    print(paste0(ChIPNames[x], " library; using ", controlNames[6], " for log2((ChIP+1)/(input+1)) calculation"))
    log2((ChIP_ranLocMats[[x]]+1)/(control_ranLocMats[[6]]+1))
  } else if ( grepl("set5", ChIPNames[x]) ) {
    print(paste0(ChIPNames[x], " library; using ", controlNames[7], " for log2((ChIP+1)/(input+1)) calculation"))
    log2((ChIP_ranLocMats[[x]]+1)/(control_ranLocMats[[7]]+1))
  } else if ( grepl("set7", ChIPNames[x]) ) {
    print(paste0(ChIPNames[x], " library; using ", controlNames[8], " for log2((ChIP+1)/(input+1)) calculation"))
    log2((ChIP_ranLocMats[[x]]+1)/(control_ranLocMats[[8]]+1))
  } else {
    print(paste0(ChIPNames[x], " library; using ", controlNames[1], " for log2((ChIP+1)/(input+1)) calculation"))
    log2((ChIP_ranLocMats[[x]]+1)/(control_ranLocMats[[1]]+1))
  }
}, mc.cores = length(ChIP_ranLocMats))


# Add column names
for(x in seq_along(log2ChIP_featureMats)) {
  colnames(log2ChIP_featureMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                           paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                           paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))
  colnames(log2ChIP_ranLocMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                          paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                          paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))
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
ymin_list_log2ChIP <- lapply(seq_along(summaryDFfeature_log2ChIP), function(x) {
  min(c(summaryDFfeature_log2ChIP[[x]][[1]]$CI_lower,
        summaryDFfeature_log2ChIP[[x]][[2]]$CI_lower,
        summaryDFfeature_log2ChIP[[x]][[3]]$CI_lower))
})
ymax_list_log2ChIP <- lapply(seq_along(summaryDFfeature_log2ChIP), function(x) {
  max(c(summaryDFfeature_log2ChIP[[x]][[1]]$CI_upper,
        summaryDFfeature_log2ChIP[[x]][[2]]$CI_upper,
        summaryDFfeature_log2ChIP[[x]][[3]]$CI_upper))
})

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
  scale_y_continuous(limits = c(ymin_list_log2ChIP[[x]], ymax_list_log2ChIP[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_log2ChIP[[x]][[1]])[1]/quantiles)-(downstream/binSize),
                              dim(summaryDFfeature_log2ChIP[[x]][[1]])[1]/quantiles),
                   labels = c(paste0("-", flankName),
                              featureStartLab,
                              featureEndLab,
                              paste0("+", flankName))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_log2ChIP[[x]][[1]])[1]/quantiles)-(downstream/binSize)),
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
  scale_y_continuous(limits = c(ymin_list_log2ChIP[[x]], ymax_list_log2ChIP[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_log2ChIP[[x]][[2]])[1]/quantiles)-(downstream/binSize),
                              dim(summaryDFfeature_log2ChIP[[x]][[2]])[1]/quantiles),
                   labels = c(paste0("-", flankName),
                              featureStartLab,
                              featureEndLab,
                              paste0("+", flankName))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_log2ChIP[[x]][[2]])[1]/quantiles)-(downstream/binSize)),
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
  scale_y_continuous(limits = c(ymin_list_log2ChIP[[x]], ymax_list_log2ChIP[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_log2ChIP[[x]][[3]])[1]/quantiles)-(downstream/binSize),
                              dim(summaryDFfeature_log2ChIP[[x]][[3]])[1]/quantiles),
                   labels = c(paste0("-", flankName),
                              "Start",
                              "End",
                              paste0("+", flankName))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_log2ChIP[[x]][[3]])[1]/quantiles)-(downstream/binSize)),
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
               "_of_CEN180_in_RaGOO_v2.0_",
               paste0(chrName, collapse = "_"), "_", align, ".pdf"),
       plot = ggObjGA_combined,
       height = 6.5*length(c(log2ChIPNamesPlot)), width = 21, limitsize = FALSE)

#### Free up memory by removing no longer required objects
rm(
   ChIP_featureMats, ChIP_ranLocMats,
   control_featureMats, control_ranLocMats,
   log2ChIP_featureMats, log2ChIP_ranLocMats,
   log2ChIP_mats_quantiles,
   wideDFfeature_list_log2ChIP,
   tidyDFfeature_list_log2ChIP,
   summaryDFfeature_list_log2ChIP,
   summaryDFfeature_log2ChIP
  ) 
gc()
#####
