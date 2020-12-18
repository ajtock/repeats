#!/applications/R/R-4.0.0/bin/Rscript

# author: Andy Tock
# contact: ajt200@cam.ac.uk
# date: 18.11.2020

# Calculate and plot metaprofiles of ChIP-seq, MNase-seq, etc.
# (CEN180 windowed means and 95% confidence intervals, CIs)
# for all CEN180 sequences and randomly positioned loci

# Usage:
# /applications/R/R-4.0.0/bin/Rscript CEN180_CENAthila_CENgap_metaprofiles.R 'Chr1,Chr2,Chr3,Chr4,Chr5' both 180 11000 1000 1kb 10 1000 '10bp' '1kb' '0.02,0.96'

#chrName <- unlist(strsplit("Chr1,Chr2,Chr3,Chr4,Chr5",
#                           split = ","))
#align <- "both"
#bodyLength <- 180
#Athila_bodyLength <- 11000
#upstream <- 1000
#downstream <- 1000
#flankName <- "1kb"
#binSize <- 10
#Athila_binSize <- 1000
#binName <- "10bp"
#Athila_binName <- "1kb"
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
align <- args[2]
bodyLength <- as.numeric(args[3])
Athila_bodyLength <- as.numeric(args[4])
upstream <- as.numeric(args[5])
downstream <- as.numeric(args[5])
flankName <- args[6]
binSize <- as.numeric(args[7])
Athila_binSize <- as.numeric(args[8])
binName <- args[9]
Athila_binName <- args[10]
legendPos <- as.numeric(unlist(strsplit(args[11],
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
#extrafont::loadfonts()

outDir <- paste0(paste0(chrName, collapse = "_"), "/")
plotDir <- paste0(outDir, "plots/")
system(paste0("[ -d ", outDir, " ] || mkdir -p ", outDir))
system(paste0("[ -d ", plotDir, " ] || mkdir -p ", plotDir))

if(length(chrName) == 5) {
  featureNamePlot <- "All CEN180"
  gapNamePlot <- "All CENgap"
  AthilaNamePlot <- "All CENAthila"
} else {
  featureNamePlot <- paste0(paste0(chrName, collapse = ","), " CEN180")
  gapNamePlot <- paste0(paste0(chrName, collapse = ","), " CENgap")
  AthilaNamePlot <- paste0(paste0(chrName, collapse = ","), " CENAthila")
}

# Define feature start and end labels for plotting
featureStartLab <- "Start"
featureEndLab <- "End"

# Load feature matrices for each chromatin dataset, calculate log2(ChIP/control),
# and sort by decreasing log2mat1RegionRowMeans
ChIPNames <- c(
               "WT_CENH3_Rep1_ChIP_SRR4430537",
               "LoCENH3_Rep1_ChIP_SRR4430543",
               "ZmCENH3_Rep1_ChIP_SRR4430549",
#               "WT_MNase_Rep1",
               "WT_REC8_HA_Rep2_ChIP",
               "WT_ASY1_Rep1_ChIP",
#               "WT_MTOPVIB_HA_Rep1_ChIP",
#               "WT_MTOPVIB_HA_Rep2_ChIP",
#               "WT_DMC1_V5_Rep1_ChIP",
               "WT_SPO11oligos_Rep1",
               "WT_SPO11oligos_Rep2",
               "WT_SPO11oligos_Rep3",
#               "met1_SPO11oligos_Rep1",
#               "met1_SPO11oligos_Rep2",
#               "met1_SPO11oligos_Rep3",
#               "kss_SPO11oligos_Rep1",
#               "kss_SPO11oligos_Rep2",
               "WT_H3K9me2_Rep1_ChIP",
               "WT_H3K4me1_Rep1_ChIP",
               "WT_H3K4me2_Rep1_ChIP",
               "WT_H3K27me1_Rep1_ChIP"
#               "WT_H3K4me3_ChIP14",
#               "WT_H3K27me3_ChIP_SRR1509478",
#               "H2A_ChIP_set1",
#               "H2AW6_ChIP_set1",
#               "H2AX_ChIP_set1",
#               "H2AZ_ChIP_set1",
#               "H3_ChIP_set1",
#               "H2AW7_ChIP_set2",
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
#               "H3K9me2_ChIP_set5"
#               "H3_ChIP_set5",
#               "H3_ChIP_set7",
#               "HTB1_ChIP_set7",
#               "HTB2_ChIP_set7",
#               "HTB3_ChIP_set7",
#               "HTB4_ChIP_set7"
              )
ChIPNamesDir <- c(
                  rep("CENH3_seedlings_Maheshwari_Comai_2017_GenomeRes/snakemake_ChIPseq_RaGOO_v2.0", 3),
#                  "150701_Kyuha_MNase/WT/snakemake_ChIPseq_RaGOO_v2.0",
                  "REC8_pooled/snakemake_ChIPseq_RaGOO_v2.0",
                  "20190722_cal66_Athaliana_ChIPseq_ASY1/fastq_pooled/snakemake_ChIPseq_RaGOO_v2.0",
#                  rep("20190819_dh580_Athaliana_ChIPseq_MTOPVIB/fastq_pooled/snakemake_ChIPseq_RaGOO_v2.0", 2),
#                  "20190917_dh580_Athaliana_ChIPseq_DMC1/fastq_pooled/snakemake_ChIPseq_RaGOO_v2.0",
                  rep("160518_Kyuha_SPO11oligos/WT/snakemake_SPO11oligos_RaGOO_v2.0", 3),
#                  rep("160518_Kyuha_SPO11oligos/met1/snakemake_SPO11oligos_RaGOO_v2.0", 3),
#                  rep("160518_Kyuha_SPO11oligos/kss/snakemake_SPO11oligos_RaGOO_v2.0", 2),
                  "170101_Chris_H3K9me2_ChIP/WT/snakemake_ChIPseq_RaGOO_v2.0",
                  rep("170920_Chris_histone_ChIP/snakemake_ChIPseq_RaGOO_v2.0", 3)
#                  "160601_Kyuha_H3K4me3_ChIP/WT/snakemake_ChIPseq_RaGOO_v2.0",
#                  "H3K27me3_bud_UWMadison_2015/snakemake_ChIPseq_RaGOO_v2.0",
#                  rep("histones_10dps_seedling_Berger_2020/snakemake_ChIPseq_RaGOO_v2.0", 17)
                 )
log2ChIPNamesPlot <- c(
                       "CENH3",
                       "LoCENH3",
                       "ZmCENH3",
#                       "MNase",
                       "REC8-HA",
                       "ASY1",
#                       "MTOPVIB-HA Rep1",
#                       "MTOPVIB-HA Rep2",
#                       "DMC1-V5",
                       "SPO11-1-oligos Rep1",
                       "SPO11-1-oligos Rep2",
                       "SPO11-1-oligos Rep3",
#                       "met1 SPO11-1-oligos Rep1",
#                       "met1 SPO11-1-oligos Rep2",
#                       "met1 SPO11-1-oligos Rep3",
#                       "kss SPO11-1-oligos Rep1",
#                       "kss SPO11-1-oligos Rep2",
                       "H3K9me2 (CL)",
                       "H3K4me1 (CL)",
                       "H3K4me2 (CL)",
                       "H3K27me1 (CL)"
#                       "H3K4me3 (KC)",
#                       "H3K27me3 (UWM)",
#                       "H2A set1",
#                       "H2A.W6 set1",
#                       "H2A.X set1",
#                       "H2A.Z set1",
#                       "H3 set1",
#                       "H2A.W7 set2",
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
#                       "H3K9me2 set5"
#                       "H3 set5",
#                       "H3 set7",
#                       "HTB1 set7",
#                       "HTB2 set7",
#                       "HTB3 set7",
#                       "HTB4 set7"
                      )
ChIPNamesPlot <- log2ChIPNamesPlot
log2ChIPColours <- c(
                     rep("red", length(log2ChIPNamesPlot))
                    )
ChIPColours <- log2ChIPColours
ChIPDirs <- sapply(seq_along(ChIPNames), function(x) {
  paste0("/home/ajt200/analysis/",
         ChIPNamesDir[x],
         "/mapped/CEN180profiles/matrices/")
})

controlNames <- c(
                  "WT_REC8_Myc_Rep1_input",
#                  "WT_gDNA_Rep1",
                  "WT_gDNA_Rep1_R1",
#                  "input_set1",
#                  "input_set2",
#                  "input_set3",
#                  "input_set5",
#                  "input_set7",
                  "map_K40_E2",
                  "map_K45_E2",
                  "map_K50_E2",
                  "map_K150_E4",
                  "map_K200_E4",
                  "map_K300_E4"
                 )
controlNamesDir <- c(
                     "REC8_pooled/snakemake_ChIPseq_RaGOO_v2.0",
#                     "150701_Natasha_gDNA/WT/snakemake_ChIPseq_RaGOO_v2.0",
                     "150701_Natasha_gDNA/WT/R1/snakemake_SPO11oligos_RaGOO_v2.0",
#                     rep("histones_10dps_seedling_Berger_2020/snakemake_ChIPseq_RaGOO_v2.0", 4),
                     rep("nanopore/RaGOO_v2.0/genmap_mappability", 6)
                    )
controlNamesPlot <- c(
                      "PE input",
#                      "Paired-end gDNA",
                      "SE gDNA",
#                      "Single-end input set1",
#                      "Single-end input set2",
#                      "Single-end input set3",
#                      "Single-end input set5",
#                      "Single-end input set7",
                      "k=40 e=2 mappability",
                      "k=45 e=2 mappability",
                      "k=50 e=2 mappability",
                      "k=150 e=4 mappability",
                      "k=200 e=4 mappability",
                      "k=300 e=4 mappability"
                     )
controlColours <- c(
                    rep("red", length(controlNamesPlot))
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
                                chrName[y], "_matrix_bin", binName, "_flank", flankName, ".tab"),
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
    if( grepl("map_K", controlNames[x]) ) {
      as.matrix(read.table(paste0(controlDirs[x],
                                  controlNames[x],
                                  "_MappedOn_Athaliana_ONT_RaGOO_v2.0_CEN180_in_",
                                  chrName[y], "_matrix_bin", binName, "_flank", flankName, ".tab"),
                           header = F, skip = 3))
    } else {
      as.matrix(read.table(paste0(controlDirs[x],
                                  controlNames[x],
                                  "_MappedOn_Athaliana_ONT_RaGOO_v2.0_lowXM_", align, "_sort_norm_CEN180_in_",
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

## ChIP
# gap
ChIP_gapMats <- mclapply(seq_along(ChIPNames), function(x) {
  lapply(seq_along(chrName), function(y) {
    as.matrix(read.table(paste0(ChIPDirs[x],
                                ChIPNames[x],
                                "_MappedOn_Athaliana_ONT_RaGOO_v2.0_lowXM_", align, "_sort_norm_CENgap_in_",
                                chrName[y], "_matrix_bin", binName, "_flank", flankName, ".tab"),
                         header = F, skip = 3))
  })
}, mc.cores = length(ChIPNames))
# If gaps from multiple chromosomes are to be analysed,
# concatenate the corresponding gap coverage matrices
ChIP_gapMats <- mclapply(seq_along(ChIP_gapMats), function(x) {
  if(length(chrName) > 1) {
    do.call(rbind, ChIP_gapMats[[x]])
  } else {
    ChIP_gapMats[[x]][[1]]
  }
}, mc.cores = length(ChIP_gapMats))

## control
# gap
control_gapMats <- mclapply(seq_along(controlNames), function(x) {
  lapply(seq_along(chrName), function(y) {
    if( grepl("map_K", controlNames[x]) ) {
      as.matrix(read.table(paste0(controlDirs[x],
                                  controlNames[x],
                                  "_MappedOn_Athaliana_ONT_RaGOO_v2.0_CENgap_in_",
                                  chrName[y], "_matrix_bin", binName, "_flank", flankName, ".tab"),
                           header = F, skip = 3))
    } else {
      as.matrix(read.table(paste0(controlDirs[x],
                                  controlNames[x],
                                  "_MappedOn_Athaliana_ONT_RaGOO_v2.0_lowXM_", align, "_sort_norm_CENgap_in_",
                                  chrName[y], "_matrix_bin", binName, "_flank", flankName, ".tab"),
                           header = F, skip = 3))
      }
  })
}, mc.cores = length(controlNames))
# If gaps from multiple chromosomes are to be analysed,
# concatenate the corresponding gap coverage matrices
control_gapMats <- mclapply(seq_along(control_gapMats), function(x) {
  if(length(chrName) > 1) {
    do.call(rbind, control_gapMats[[x]])
  } else {
    control_gapMats[[x]][[1]]
  }
}, mc.cores = length(control_gapMats))

## ChIP
# Athila
ChIP_AthilaMats <- mclapply(seq_along(ChIPNames), function(x) {
  lapply(seq_along(chrName), function(y) {
    as.matrix(read.table(paste0(ChIPDirs[x],
                                ChIPNames[x],
                                "_MappedOn_Athaliana_ONT_RaGOO_v2.0_lowXM_", align, "_sort_norm_CENAthila_in_",
                                chrName[y], "_matrix_bin", binName, "_flank", flankName, ".tab"),
                         header = F, skip = 3))
  })
}, mc.cores = length(ChIPNames))
# If Athilas from multiple chromosomes are to be analysed,
# concatenate the corresponding Athila coverage matrices
ChIP_AthilaMats <- mclapply(seq_along(ChIP_AthilaMats), function(x) {
  if(length(chrName) > 1) {
    do.call(rbind, ChIP_AthilaMats[[x]])
  } else {
    ChIP_AthilaMats[[x]][[1]]
  }
}, mc.cores = length(ChIP_AthilaMats))

## control
# Athila
control_AthilaMats <- mclapply(seq_along(controlNames), function(x) {
  lapply(seq_along(chrName), function(y) {
    if( grepl("map_K", controlNames[x]) ) {
      as.matrix(read.table(paste0(controlDirs[x],
                                  controlNames[x],
                                  "_MappedOn_Athaliana_ONT_RaGOO_v2.0_CENAthila_in_",
                                  chrName[y], "_matrix_bin", binName, "_flank", flankName, ".tab"),
                           header = F, skip = 3))
    } else {
      as.matrix(read.table(paste0(controlDirs[x],
                                  controlNames[x],
                                  "_MappedOn_Athaliana_ONT_RaGOO_v2.0_lowXM_", align, "_sort_norm_CENAthila_in_",
                                  chrName[y], "_matrix_bin", binName, "_flank", flankName, ".tab"),
                           header = F, skip = 3))
      }
  })
}, mc.cores = length(controlNames))
# If Athilas from multiple chromosomes are to be analysed,
# concatenate the corresponding Athila coverage matrices
control_AthilaMats <- mclapply(seq_along(control_AthilaMats), function(x) {
  if(length(chrName) > 1) {
    do.call(rbind, control_AthilaMats[[x]])
  } else {
    control_AthilaMats[[x]][[1]]
  }
}, mc.cores = length(control_AthilaMats))

# Conditionally calculate log2(ChIP/control)
# for each matrix depending on library
# feature
log2ChIP_featureMats <- mclapply(seq_along(ChIP_featureMats), function(x) {
  if ( grepl("MNase", ChIPNames[x]) ) {
    print(paste0(ChIPNames[x], " library; using ", controlNames[1], " for log2((MNase+1)/(gDNA+1)) calculation"))
    log2((ChIP_featureMats[[x]]+1)/(control_featureMats[[1]]+1))
  } else if ( grepl("SPO11oligos", ChIPNames[x]) ) {
    print(paste0(ChIPNames[x], " library; using ", controlNames[2], " for log2((SPO11-1-oligos+1)/(gDNA+1)) calculation"))
    log2((ChIP_featureMats[[x]]+1)/(control_featureMats[[2]]+1))
#  } else if ( grepl("set1", ChIPNames[x]) ) {
#    print(paste0(ChIPNames[x], " library; using ", controlNames[4], " for log2((ChIP+1)/(input+1)) calculation"))
#    log2((ChIP_featureMats[[x]]+1)/(control_featureMats[[4]]+1))
#  } else if ( grepl("set2", ChIPNames[x]) ) {
#    print(paste0(ChIPNames[x], " library; using ", controlNames[5], " for log2((ChIP+1)/(input+1)) calculation"))
#    log2((ChIP_featureMats[[x]]+1)/(control_featureMats[[5]]+1))
#  } else if ( grepl("set3", ChIPNames[x]) ) {
#    print(paste0(ChIPNames[x], " library; using ", controlNames[6], " for log2((ChIP+1)/(input+1)) calculation"))
#    log2((ChIP_featureMats[[x]]+1)/(control_featureMats[[6]]+1))
#  } else if ( grepl("set5", ChIPNames[x]) ) {
#    print(paste0(ChIPNames[x], " library; using ", controlNames[7], " for log2((ChIP+1)/(input+1)) calculation"))
#    log2((ChIP_featureMats[[x]]+1)/(control_featureMats[[7]]+1))
#  } else if ( grepl("set7", ChIPNames[x]) ) {
#    print(paste0(ChIPNames[x], " library; using ", controlNames[8], " for log2((ChIP+1)/(input+1)) calculation"))
#    log2((ChIP_featureMats[[x]]+1)/(control_featureMats[[8]]+1))
  } else {
    print(paste0(ChIPNames[x], " library; using ", controlNames[1], " for log2((ChIP+1)/(input+1)) calculation"))
    log2((ChIP_featureMats[[x]]+1)/(control_featureMats[[1]]+1))
  }
}, mc.cores = length(ChIP_featureMats))

# Conditionally calculate log2(ChIP/control)
# for each matrix depending on library
# gap
log2ChIP_gapMats <- mclapply(seq_along(ChIP_gapMats), function(x) {
  if ( grepl("MNase", ChIPNames[x]) ) {
    print(paste0(ChIPNames[x], " library; using ", controlNames[1], " for log2((MNase+1)/(gDNA+1)) calculation"))
    log2((ChIP_gapMats[[x]]+1)/(control_gapMats[[1]]+1))
  } else if ( grepl("SPO11oligos", ChIPNames[x]) ) {
    print(paste0(ChIPNames[x], " library; using ", controlNames[2], " for log2((SPO11-1-oligos+1)/(gDNA+1)) calculation"))
    log2((ChIP_gapMats[[x]]+1)/(control_gapMats[[2]]+1))
#  } else if ( grepl("set1", ChIPNames[x]) ) {
#    print(paste0(ChIPNames[x], " library; using ", controlNames[4], " for log2((ChIP+1)/(input+1)) calculation"))
#    log2((ChIP_gapMats[[x]]+1)/(control_gapMats[[4]]+1))
#  } else if ( grepl("set2", ChIPNames[x]) ) {
#    print(paste0(ChIPNames[x], " library; using ", controlNames[5], " for log2((ChIP+1)/(input+1)) calculation"))
#    log2((ChIP_gapMats[[x]]+1)/(control_gapMats[[5]]+1))
#  } else if ( grepl("set3", ChIPNames[x]) ) {
#    print(paste0(ChIPNames[x], " library; using ", controlNames[6], " for log2((ChIP+1)/(input+1)) calculation"))
#    log2((ChIP_gapMats[[x]]+1)/(control_gapMats[[6]]+1))
#  } else if ( grepl("set5", ChIPNames[x]) ) {
#    print(paste0(ChIPNames[x], " library; using ", controlNames[7], " for log2((ChIP+1)/(input+1)) calculation"))
#    log2((ChIP_gapMats[[x]]+1)/(control_gapMats[[7]]+1))
#  } else if ( grepl("set7", ChIPNames[x]) ) {
#    print(paste0(ChIPNames[x], " library; using ", controlNames[8], " for log2((ChIP+1)/(input+1)) calculation"))
#    log2((ChIP_gapMats[[x]]+1)/(control_gapMats[[8]]+1))
  } else {
    print(paste0(ChIPNames[x], " library; using ", controlNames[1], " for log2((ChIP+1)/(input+1)) calculation"))
    log2((ChIP_gapMats[[x]]+1)/(control_gapMats[[1]]+1))
  }
}, mc.cores = length(ChIP_gapMats))

# Conditionally calculate log2(ChIP/control)
# for each matrix depending on library
# Athila
log2ChIP_AthilaMats <- mclapply(seq_along(ChIP_AthilaMats), function(x) {
  if ( grepl("MNase", ChIPNames[x]) ) {
    print(paste0(ChIPNames[x], " library; using ", controlNames[1], " for log2((MNase+1)/(gDNA+1)) calculation"))
    log2((ChIP_AthilaMats[[x]]+1)/(control_AthilaMats[[1]]+1))
  } else if ( grepl("SPO11oligos", ChIPNames[x]) ) {
    print(paste0(ChIPNames[x], " library; using ", controlNames[2], " for log2((SPO11-1-oligos+1)/(gDNA+1)) calculation"))
    log2((ChIP_AthilaMats[[x]]+1)/(control_AthilaMats[[2]]+1))
#  } else if ( grepl("set1", ChIPNames[x]) ) {
#    print(paste0(ChIPNames[x], " library; using ", controlNames[4], " for log2((ChIP+1)/(input+1)) calculation"))
#    log2((ChIP_AthilaMats[[x]]+1)/(control_AthilaMats[[4]]+1))
#  } else if ( grepl("set2", ChIPNames[x]) ) {
#    print(paste0(ChIPNames[x], " library; using ", controlNames[5], " for log2((ChIP+1)/(input+1)) calculation"))
#    log2((ChIP_AthilaMats[[x]]+1)/(control_AthilaMats[[5]]+1))
#  } else if ( grepl("set3", ChIPNames[x]) ) {
#    print(paste0(ChIPNames[x], " library; using ", controlNames[6], " for log2((ChIP+1)/(input+1)) calculation"))
#    log2((ChIP_AthilaMats[[x]]+1)/(control_AthilaMats[[6]]+1))
#  } else if ( grepl("set5", ChIPNames[x]) ) {
#    print(paste0(ChIPNames[x], " library; using ", controlNames[7], " for log2((ChIP+1)/(input+1)) calculation"))
#    log2((ChIP_AthilaMats[[x]]+1)/(control_AthilaMats[[7]]+1))
#  } else if ( grepl("set7", ChIPNames[x]) ) {
#    print(paste0(ChIPNames[x], " library; using ", controlNames[8], " for log2((ChIP+1)/(input+1)) calculation"))
#    log2((ChIP_AthilaMats[[x]]+1)/(control_AthilaMats[[8]]+1))
  } else {
    print(paste0(ChIPNames[x], " library; using ", controlNames[1], " for log2((ChIP+1)/(input+1)) calculation"))
    log2((ChIP_AthilaMats[[x]]+1)/(control_AthilaMats[[1]]+1))
  }
}, mc.cores = length(ChIP_AthilaMats))

# log2ChIP
# Add column names
for(x in seq_along(log2ChIP_featureMats)) {
  colnames(log2ChIP_featureMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                           paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                           paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))
  colnames(log2ChIP_gapMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                       paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                       paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))
  colnames(log2ChIP_AthilaMats[[x]]) <- c(paste0("u", 1:(upstream/Athila_binSize)),
                                          paste0("t", ((upstream/Athila_binSize)+1):((upstream+Athila_bodyLength)/Athila_binSize)),
                                          paste0("d", (((upstream+Athila_bodyLength)/Athila_binSize)+1):(((upstream+Athila_bodyLength)/Athila_binSize)+(downstream/Athila_binSize))))
}

# Create list of lists in which each element in the enclosing list corresponds to a library
# and the two elements in the nested list correspond to coverage matrices for features and random loci 
log2ChIP_mats <- mclapply(seq_along(log2ChIP_featureMats), function(x) {
  list(
       # features
       log2ChIP_featureMats[[x]],
       # gaps
       log2ChIP_gapMats[[x]],
       # Athilas
       log2ChIP_AthilaMats[[x]]
      ) 
}, mc.cores = length(log2ChIP_featureMats))

# Transpose matrix and convert into dataframe
# in which first column is window name
wideDFfeature_list_log2ChIP <- mclapply(seq_along(log2ChIP_mats), function(x) {
  lapply(seq_along(log2ChIP_mats[[x]]), function(y) {
    data.frame(window = colnames(log2ChIP_mats[[x]][[y]]),
               t(log2ChIP_mats[[x]][[y]]))
  })
}, mc.cores = length(log2ChIP_mats))

# Convert into tidy data.frame (long format)
tidyDFfeature_list_log2ChIP  <- mclapply(seq_along(wideDFfeature_list_log2ChIP), function(x) {
  lapply(seq_along(log2ChIP_mats[[x]]), function(y) {
    gather(data  = wideDFfeature_list_log2ChIP[[x]][[y]],
           key   = feature,
           value = coverage,
           -window)
  }) 
}, mc.cores = length(wideDFfeature_list_log2ChIP))

# Order levels of factor "window" so that sequential levels
# correspond to sequential windows
for(x in seq_along(tidyDFfeature_list_log2ChIP)) {
  for(y in seq_along(log2ChIP_mats[[x]])) {
    tidyDFfeature_list_log2ChIP[[x]][[y]]$window <- factor(tidyDFfeature_list_log2ChIP[[x]][[y]]$window,
                                                           levels = as.character(wideDFfeature_list_log2ChIP[[x]][[y]]$window))
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
  lapply(seq_along(log2ChIP_mats[[x]]), function(y) {
    data.frame(window = as.character(wideDFfeature_list_log2ChIP[[x]][[y]]$window),
               n      = tapply(X     = tidyDFfeature_list_log2ChIP[[x]][[y]]$coverage,
                               INDEX = tidyDFfeature_list_log2ChIP[[x]][[y]]$window,
                               FUN   = length),
               mean   = tapply(X     = tidyDFfeature_list_log2ChIP[[x]][[y]]$coverage,
                               INDEX = tidyDFfeature_list_log2ChIP[[x]][[y]]$window,
                               FUN   = mean,
                               na.rm = TRUE),
               sd     = tapply(X     = tidyDFfeature_list_log2ChIP[[x]][[y]]$coverage,
                               INDEX = tidyDFfeature_list_log2ChIP[[x]][[y]]$window,
                               FUN   = sd,
                               na.rm = TRUE))
  })
}, mc.cores = length(tidyDFfeature_list_log2ChIP))

for(x in seq_along(summaryDFfeature_list_log2ChIP)) {
  for(y in seq_along(log2ChIP_mats[[x]])) {
    summaryDFfeature_list_log2ChIP[[x]][[y]]$window <- factor(summaryDFfeature_list_log2ChIP[[x]][[y]]$window,
                                                              levels = as.character(wideDFfeature_list_log2ChIP[[x]][[y]]$window))
    summaryDFfeature_list_log2ChIP[[x]][[y]]$winNo <- factor(1:dim(summaryDFfeature_list_log2ChIP[[x]][[y]])[1])
    summaryDFfeature_list_log2ChIP[[x]][[y]]$sem <- summaryDFfeature_list_log2ChIP[[x]][[y]]$sd/sqrt(summaryDFfeature_list_log2ChIP[[x]][[y]]$n-1)
    summaryDFfeature_list_log2ChIP[[x]][[y]]$CI_lower <- summaryDFfeature_list_log2ChIP[[x]][[y]]$mean -
      qt(0.975, df = summaryDFfeature_list_log2ChIP[[x]][[y]]$n-1)*summaryDFfeature_list_log2ChIP[[x]][[y]]$sem
    summaryDFfeature_list_log2ChIP[[x]][[y]]$CI_upper <- summaryDFfeature_list_log2ChIP[[x]][[y]]$mean +
      qt(0.975, df = summaryDFfeature_list_log2ChIP[[x]][[y]]$n-1)*summaryDFfeature_list_log2ChIP[[x]][[y]]$sem
  }
}

summaryDFfeature_log2ChIP <- summaryDFfeature_list_log2ChIP

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

# Plot average profiles with 95% CI ribbon
## feature
ggObj1_combined_log2ChIP <- mclapply(seq_along(log2ChIPNamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_log2ChIP[[x]][[1]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = 1)
        ) +
  geom_line(data = summaryDFfeature,
            colour = log2ChIPColours[x],
            size = 1) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper),
              fill = log2ChIPColours[x],
              alpha = 0.4) +
  scale_y_continuous(limits = c(ymin_list_log2ChIP[[x]], ymax_list_log2ChIP[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_log2ChIP[[x]][[1]])[1])-(downstream/binSize),
                              dim(summaryDFfeature_log2ChIP[[x]][[1]])[1]),
                   labels = c(paste0("-", flankName),
                              featureStartLab,
                              featureEndLab,
                              paste0("+", flankName))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_log2ChIP[[x]][[1]])[1])-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = log2ChIPNamesPlot[x]) +
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
        plot.title = element_text(hjust = 0.5, size = 30)) +
  ggtitle(bquote(.(featureNamePlot) ~ "(" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFfeature$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
}, mc.cores = length(log2ChIPNamesPlot))

## gap
ggObj2_combined_log2ChIP <- mclapply(seq_along(log2ChIPNamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_log2ChIP[[x]][[2]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = 1)
        ) +
  geom_line(data = summaryDFfeature,
            colour = log2ChIPColours[x],
            size = 1) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper),
              fill = log2ChIPColours[x],
              alpha = 0.4) +
  scale_y_continuous(limits = c(ymin_list_log2ChIP[[x]], ymax_list_log2ChIP[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_log2ChIP[[x]][[2]])[1])-(downstream/binSize),
                              dim(summaryDFfeature_log2ChIP[[x]][[2]])[1]),
                   labels = c(paste0("-", flankName),
                              featureStartLab,
                              featureEndLab,
                              paste0("+", flankName))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_log2ChIP[[x]][[2]])[1])-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = log2ChIPNamesPlot[x]) +
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
        plot.title = element_text(hjust = 0.5, size = 30)) +
  ggtitle(bquote(.(gapNamePlot) ~ "(" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFfeature$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
}, mc.cores = length(log2ChIPNamesPlot))

## Athila
ggObj3_combined_log2ChIP <- mclapply(seq_along(log2ChIPNamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_log2ChIP[[x]][[3]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = 1)
        ) +
  geom_line(data = summaryDFfeature,
            colour = log2ChIPColours[x],
            size = 1) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper),
              fill = log2ChIPColours[x],
              alpha = 0.4) +
  scale_y_continuous(limits = c(ymin_list_log2ChIP[[x]], ymax_list_log2ChIP[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/Athila_binSize)+1,
                              (dim(summaryDFfeature_log2ChIP[[x]][[3]])[1])-(downstream/Athila_binSize),
                              dim(summaryDFfeature_log2ChIP[[x]][[3]])[1]),
                   labels = c(paste0("-", flankName),
                              featureStartLab,
                              featureEndLab,
                              paste0("+", flankName))) +
  geom_vline(xintercept = c((upstream/Athila_binSize)+1,
                            (dim(summaryDFfeature_log2ChIP[[x]][[3]])[1])-(downstream/Athila_binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = log2ChIPNamesPlot[x]) +
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
        plot.title = element_text(hjust = 0.5, size = 30)) +
  ggtitle(bquote(.(AthilaNamePlot) ~ "(" * italic("n") ~ "=" ~
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
              "log2ChIPcontrol_avgProfiles_around",
               "_CEN180_CENgap_CENAthila_in_RaGOO_v2.0_",
               paste0(chrName, collapse = "_"), "_", align, ".pdf"),
       plot = ggObjGA_combined,
       height = 6.5*length(c(log2ChIPNamesPlot)), width = 30, limitsize = FALSE)

#### Free up memory by removing no longer required objects
rm(
   log2ChIP_featureMats, log2ChIP_AthilaMats,
   log2ChIP_mats,
   wideDFfeature_list_log2ChIP,
   tidyDFfeature_list_log2ChIP,
   summaryDFfeature_list_log2ChIP
#   summaryDFfeature_log2ChIP
  ) 
gc()
#####

# control
# Add column names
for(x in seq_along(control_featureMats)) {
  colnames(control_featureMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                           paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                           paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))
  colnames(control_gapMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                       paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                       paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))
  colnames(control_AthilaMats[[x]]) <- c(paste0("u", 1:(upstream/Athila_binSize)),
                                          paste0("t", ((upstream/Athila_binSize)+1):((upstream+Athila_bodyLength)/Athila_binSize)),
                                          paste0("d", (((upstream+Athila_bodyLength)/Athila_binSize)+1):(((upstream+Athila_bodyLength)/Athila_binSize)+(downstream/Athila_binSize))))
}

# Create list of lists in which each element in the enclosing list corresponds to a library
# and the two elements in the nested list correspond to coverage matrices for features and random loci 
control_mats <- mclapply(seq_along(control_featureMats), function(x) {
  list(
       # features
       control_featureMats[[x]],
       # gaps
       control_gapMats[[x]],
       # ranLocs
       control_AthilaMats[[x]]
      ) 
}, mc.cores = length(control_featureMats))

# Transpose matrix and convert into dataframe
# in which first column is window name
wideDFfeature_list_control <- mclapply(seq_along(control_mats), function(x) {
  lapply(seq_along(control_mats[[x]]), function(y) {
    data.frame(window = colnames(control_mats[[x]][[y]]),
               t(control_mats[[x]][[y]]))
  })
}, mc.cores = length(control_mats))

# Convert into tidy data.frame (long format)
tidyDFfeature_list_control  <- mclapply(seq_along(wideDFfeature_list_control), function(x) {
  lapply(seq_along(control_mats[[x]]), function(y) {
    gather(data  = wideDFfeature_list_control[[x]][[y]],
           key   = feature,
           value = coverage,
           -window)
  }) 
}, mc.cores = length(wideDFfeature_list_control))

# Order levels of factor "window" so that sequential levels
# correspond to sequential windows
for(x in seq_along(tidyDFfeature_list_control)) {
  for(y in seq_along(control_mats[[x]])) {
    tidyDFfeature_list_control[[x]][[y]]$window <- factor(tidyDFfeature_list_control[[x]][[y]]$window,
                                                           levels = as.character(wideDFfeature_list_control[[x]][[y]]$window))
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
  lapply(seq_along(control_mats[[x]]), function(y) {
    data.frame(window = as.character(wideDFfeature_list_control[[x]][[y]]$window),
               n      = tapply(X     = tidyDFfeature_list_control[[x]][[y]]$coverage,
                               INDEX = tidyDFfeature_list_control[[x]][[y]]$window,
                               FUN   = length),
               mean   = tapply(X     = tidyDFfeature_list_control[[x]][[y]]$coverage,
                               INDEX = tidyDFfeature_list_control[[x]][[y]]$window,
                               FUN   = mean,
                               na.rm = TRUE),
               sd     = tapply(X     = tidyDFfeature_list_control[[x]][[y]]$coverage,
                               INDEX = tidyDFfeature_list_control[[x]][[y]]$window,
                               FUN   = sd,
                               na.rm = TRUE))
  })
}, mc.cores = length(tidyDFfeature_list_control))

for(x in seq_along(summaryDFfeature_list_control)) {
  for(y in seq_along(control_mats[[x]])) {
    summaryDFfeature_list_control[[x]][[y]]$window <- factor(summaryDFfeature_list_control[[x]][[y]]$window,
                                                              levels = as.character(wideDFfeature_list_control[[x]][[y]]$window))
    summaryDFfeature_list_control[[x]][[y]]$winNo <- factor(1:dim(summaryDFfeature_list_control[[x]][[y]])[1])
    summaryDFfeature_list_control[[x]][[y]]$sem <- summaryDFfeature_list_control[[x]][[y]]$sd/sqrt(summaryDFfeature_list_control[[x]][[y]]$n-1)
    summaryDFfeature_list_control[[x]][[y]]$CI_lower <- summaryDFfeature_list_control[[x]][[y]]$mean -
      qt(0.975, df = summaryDFfeature_list_control[[x]][[y]]$n-1)*summaryDFfeature_list_control[[x]][[y]]$sem
    summaryDFfeature_list_control[[x]][[y]]$CI_upper <- summaryDFfeature_list_control[[x]][[y]]$mean +
      qt(0.975, df = summaryDFfeature_list_control[[x]][[y]]$n-1)*summaryDFfeature_list_control[[x]][[y]]$sem
  }
}

summaryDFfeature_control <- summaryDFfeature_list_control

# Define y-axis limits
ymin_list_control <- lapply(seq_along(summaryDFfeature_control), function(x) {
  min(c(summaryDFfeature_control[[x]][[1]]$CI_lower,
        summaryDFfeature_control[[x]][[2]]$CI_lower,
        summaryDFfeature_control[[x]][[3]]$CI_lower))
})
ymax_list_control <- lapply(seq_along(summaryDFfeature_control), function(x) {
  max(c(summaryDFfeature_control[[x]][[1]]$CI_upper,
        summaryDFfeature_control[[x]][[2]]$CI_upper,
        summaryDFfeature_control[[x]][[3]]$CI_upper))
})

# Plot average profiles with 95% CI ribbon
## feature
ggObj1_combined_control <- mclapply(seq_along(controlNamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_control[[x]][[1]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = 1)
        ) +
  geom_line(data = summaryDFfeature,
            colour = controlColours[x],
            size = 1) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper),
              fill = controlColours[x],
              alpha = 0.4) +
  scale_y_continuous(limits = c(ymin_list_control[[x]], ymax_list_control[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_control[[x]][[1]])[1])-(downstream/binSize),
                              dim(summaryDFfeature_control[[x]][[1]])[1]),
                   labels = c(paste0("-", flankName),
                              featureStartLab,
                              featureEndLab,
                              paste0("+", flankName))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_control[[x]][[1]])[1])-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = controlNamesPlot[x]) +
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
        plot.title = element_text(hjust = 0.5, size = 30)) +
  ggtitle(bquote(.(featureNamePlot) ~ "(" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFfeature$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
}, mc.cores = length(controlNamesPlot))

## gap
ggObj2_combined_control <- mclapply(seq_along(controlNamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_control[[x]][[2]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = 1)
        ) +
  geom_line(data = summaryDFfeature,
            colour = controlColours[x],
            size = 1) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper),
              fill = controlColours[x],
              alpha = 0.4) +
  scale_y_continuous(limits = c(ymin_list_control[[x]], ymax_list_control[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_control[[x]][[2]])[1])-(downstream/binSize),
                              dim(summaryDFfeature_control[[x]][[2]])[1]),
                   labels = c(paste0("-", flankName),
                              featureStartLab,
                              featureEndLab,
                              paste0("+", flankName))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_control[[x]][[2]])[1])-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = controlNamesPlot[x]) +
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
        plot.title = element_text(hjust = 0.5, size = 30)) +
  ggtitle(bquote(.(gapNamePlot) ~ "(" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFfeature$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
}, mc.cores = length(controlNamesPlot))

## Athila
ggObj3_combined_control <- mclapply(seq_along(controlNamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_control[[x]][[3]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = 1)
        ) +
  geom_line(data = summaryDFfeature,
            colour = controlColours[x],
            size = 1) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper),
              fill = controlColours[x],
              alpha = 0.4) +
  scale_y_continuous(limits = c(ymin_list_control[[x]], ymax_list_control[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/Athila_binSize)+1,
                              (dim(summaryDFfeature_control[[x]][[3]])[1])-(downstream/Athila_binSize),
                              dim(summaryDFfeature_control[[x]][[3]])[1]),
                   labels = c(paste0("-", flankName),
                              featureStartLab,
                              featureEndLab,
                              paste0("+", flankName))) +
  geom_vline(xintercept = c((upstream/Athila_binSize)+1,
                            (dim(summaryDFfeature_control[[x]][[3]])[1])-(downstream/Athila_binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = controlNamesPlot[x]) +
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
        plot.title = element_text(hjust = 0.5, size = 30)) +
  ggtitle(bquote(.(AthilaNamePlot) ~ "(" * italic("n") ~ "=" ~
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
              "control_avgProfiles_around",
               "_CEN180_CENgap_CENAthila_in_RaGOO_v2.0_",
               paste0(chrName, collapse = "_"), "_", align, ".pdf"),
       plot = ggObjGA_combined,
       height = 6.5*length(c(controlNamesPlot)), width = 30, limitsize = FALSE)

#### Free up memory by removing no longer required objects
rm(
   control_featureMats, control_AthilaMats,
   control_mats,
   wideDFfeature_list_control,
   tidyDFfeature_list_control,
   summaryDFfeature_list_control
#   summaryDFfeature_control
  ) 
gc()
#####


# ChIP
# Add column names
for(x in seq_along(ChIP_featureMats)) {
  colnames(ChIP_featureMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                           paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                           paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))
  colnames(ChIP_gapMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                       paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                       paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))
  colnames(ChIP_AthilaMats[[x]]) <- c(paste0("u", 1:(upstream/Athila_binSize)),
                                          paste0("t", ((upstream/Athila_binSize)+1):((upstream+Athila_bodyLength)/Athila_binSize)),
                                          paste0("d", (((upstream+Athila_bodyLength)/Athila_binSize)+1):(((upstream+Athila_bodyLength)/Athila_binSize)+(downstream/Athila_binSize))))
}

# Create list of lists in which each element in the enclosing list corresponds to a library
# and the two elements in the nested list correspond to coverage matrices for features and random loci 
ChIP_mats <- mclapply(seq_along(ChIP_featureMats), function(x) {
  list(
       # features
       ChIP_featureMats[[x]],
       # gaps
       ChIP_gapMats[[x]],
       # ranLocs
       ChIP_AthilaMats[[x]]
      ) 
}, mc.cores = length(ChIP_featureMats))

# Transpose matrix and convert into dataframe
# in which first column is window name
wideDFfeature_list_ChIP <- mclapply(seq_along(ChIP_mats), function(x) {
  lapply(seq_along(ChIP_mats[[x]]), function(y) {
    data.frame(window = colnames(ChIP_mats[[x]][[y]]),
               t(ChIP_mats[[x]][[y]]))
  })
}, mc.cores = length(ChIP_mats))

# Convert into tidy data.frame (long format)
tidyDFfeature_list_ChIP  <- mclapply(seq_along(wideDFfeature_list_ChIP), function(x) {
  lapply(seq_along(ChIP_mats[[x]]), function(y) {
    gather(data  = wideDFfeature_list_ChIP[[x]][[y]],
           key   = feature,
           value = coverage,
           -window)
  }) 
}, mc.cores = length(wideDFfeature_list_ChIP))

# Order levels of factor "window" so that sequential levels
# correspond to sequential windows
for(x in seq_along(tidyDFfeature_list_ChIP)) {
  for(y in seq_along(ChIP_mats[[x]])) {
    tidyDFfeature_list_ChIP[[x]][[y]]$window <- factor(tidyDFfeature_list_ChIP[[x]][[y]]$window,
                                                           levels = as.character(wideDFfeature_list_ChIP[[x]][[y]]$window))
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
  lapply(seq_along(ChIP_mats[[x]]), function(y) {
    data.frame(window = as.character(wideDFfeature_list_ChIP[[x]][[y]]$window),
               n      = tapply(X     = tidyDFfeature_list_ChIP[[x]][[y]]$coverage,
                               INDEX = tidyDFfeature_list_ChIP[[x]][[y]]$window,
                               FUN   = length),
               mean   = tapply(X     = tidyDFfeature_list_ChIP[[x]][[y]]$coverage,
                               INDEX = tidyDFfeature_list_ChIP[[x]][[y]]$window,
                               FUN   = mean,
                               na.rm = TRUE),
               sd     = tapply(X     = tidyDFfeature_list_ChIP[[x]][[y]]$coverage,
                               INDEX = tidyDFfeature_list_ChIP[[x]][[y]]$window,
                               FUN   = sd,
                               na.rm = TRUE))
  })
}, mc.cores = length(tidyDFfeature_list_ChIP))

for(x in seq_along(summaryDFfeature_list_ChIP)) {
  for(y in seq_along(ChIP_mats[[x]])) {
    summaryDFfeature_list_ChIP[[x]][[y]]$window <- factor(summaryDFfeature_list_ChIP[[x]][[y]]$window,
                                                              levels = as.character(wideDFfeature_list_ChIP[[x]][[y]]$window))
    summaryDFfeature_list_ChIP[[x]][[y]]$winNo <- factor(1:dim(summaryDFfeature_list_ChIP[[x]][[y]])[1])
    summaryDFfeature_list_ChIP[[x]][[y]]$sem <- summaryDFfeature_list_ChIP[[x]][[y]]$sd/sqrt(summaryDFfeature_list_ChIP[[x]][[y]]$n-1)
    summaryDFfeature_list_ChIP[[x]][[y]]$CI_lower <- summaryDFfeature_list_ChIP[[x]][[y]]$mean -
      qt(0.975, df = summaryDFfeature_list_ChIP[[x]][[y]]$n-1)*summaryDFfeature_list_ChIP[[x]][[y]]$sem
    summaryDFfeature_list_ChIP[[x]][[y]]$CI_upper <- summaryDFfeature_list_ChIP[[x]][[y]]$mean +
      qt(0.975, df = summaryDFfeature_list_ChIP[[x]][[y]]$n-1)*summaryDFfeature_list_ChIP[[x]][[y]]$sem
  }
}

summaryDFfeature_ChIP <- summaryDFfeature_list_ChIP

# Define y-axis limits
ymin_list_ChIP <- lapply(seq_along(summaryDFfeature_ChIP), function(x) {
  min(c(summaryDFfeature_ChIP[[x]][[1]]$CI_lower,
        summaryDFfeature_ChIP[[x]][[2]]$CI_lower,
        summaryDFfeature_ChIP[[x]][[3]]$CI_lower))
})
ymax_list_ChIP <- lapply(seq_along(summaryDFfeature_ChIP), function(x) {
  max(c(summaryDFfeature_ChIP[[x]][[1]]$CI_upper,
        summaryDFfeature_ChIP[[x]][[2]]$CI_upper,
        summaryDFfeature_ChIP[[x]][[3]]$CI_upper))
})

# Plot average profiles with 95% CI ribbon
## feature
ggObj1_combined_ChIP <- mclapply(seq_along(ChIPNamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_ChIP[[x]][[1]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = 1)
        ) +
  geom_line(data = summaryDFfeature,
            colour = ChIPColours[x],
            size = 1) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper),
              fill = ChIPColours[x],
              alpha = 0.4) +
  scale_y_continuous(limits = c(ymin_list_ChIP[[x]], ymax_list_ChIP[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_ChIP[[x]][[1]])[1])-(downstream/binSize),
                              dim(summaryDFfeature_ChIP[[x]][[1]])[1]),
                   labels = c(paste0("-", flankName),
                              featureStartLab,
                              featureEndLab,
                              paste0("+", flankName))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_ChIP[[x]][[1]])[1])-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = ChIPNamesPlot[x]) +
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
        plot.title = element_text(hjust = 0.5, size = 30)) +
  ggtitle(bquote(.(featureNamePlot) ~ "(" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFfeature$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
}, mc.cores = length(ChIPNamesPlot))

## gap
ggObj2_combined_ChIP <- mclapply(seq_along(ChIPNamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_ChIP[[x]][[2]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = 1)
        ) +
  geom_line(data = summaryDFfeature,
            colour = ChIPColours[x],
            size = 1) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper),
              fill = ChIPColours[x],
              alpha = 0.4) +
  scale_y_continuous(limits = c(ymin_list_ChIP[[x]], ymax_list_ChIP[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_ChIP[[x]][[2]])[1])-(downstream/binSize),
                              dim(summaryDFfeature_ChIP[[x]][[2]])[1]),
                   labels = c(paste0("-", flankName),
                              featureStartLab,
                              featureEndLab,
                              paste0("+", flankName))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_ChIP[[x]][[2]])[1])-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = ChIPNamesPlot[x]) +
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
        plot.title = element_text(hjust = 0.5, size = 30)) +
  ggtitle(bquote(.(gapNamePlot) ~ "(" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFfeature$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
}, mc.cores = length(ChIPNamesPlot))

## Athila
ggObj3_combined_ChIP <- mclapply(seq_along(ChIPNamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_ChIP[[x]][[3]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = 1)
        ) +
  geom_line(data = summaryDFfeature,
            colour = ChIPColours[x],
            size = 1) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper),
              fill = ChIPColours[x],
              alpha = 0.4) +
  scale_y_continuous(limits = c(ymin_list_ChIP[[x]], ymax_list_ChIP[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/Athila_binSize)+1,
                              (dim(summaryDFfeature_ChIP[[x]][[3]])[1])-(downstream/Athila_binSize),
                              dim(summaryDFfeature_ChIP[[x]][[3]])[1]),
                   labels = c(paste0("-", flankName),
                              featureStartLab,
                              featureEndLab,
                              paste0("+", flankName))) +
  geom_vline(xintercept = c((upstream/Athila_binSize)+1,
                            (dim(summaryDFfeature_ChIP[[x]][[3]])[1])-(downstream/Athila_binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = ChIPNamesPlot[x]) +
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
        plot.title = element_text(hjust = 0.5, size = 30)) +
  ggtitle(bquote(.(AthilaNamePlot) ~ "(" * italic("n") ~ "=" ~
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
              "ChIP_avgProfiles_around",
               "_CEN180_CENgap_CENAthila_in_RaGOO_v2.0_",
               paste0(chrName, collapse = "_"), "_", align, ".pdf"),
       plot = ggObjGA_combined,
       height = 6.5*length(c(ChIPNamesPlot)), width = 30, limitsize = FALSE)

#### Free up memory by removing no longer required objects
rm(
   ChIP_featureMats, ChIP_AthilaMats,
   ChIP_mats,
   wideDFfeature_list_ChIP,
   tidyDFfeature_list_ChIP,
   summaryDFfeature_list_ChIP
#   summaryDFfeature_ChIP
  ) 
gc()
#####


## DNAmeth
DNAmethNames <- c(
                  "WT_nanopolishDNAmeth_95_10kb", 
                  rep("WT_BSseq_Rep1_2014", 3)
                 )
DNAmethNamesDir <- c(
                     "nanopore/RaGOO_v2.0/nanopolish_DNAmeth",
                     rep("BSseq_leaf_Stroud_Jacobsen_2013_Cell_2014_NSMB/snakemake_BSseq_RaGOO_v2.0/coverage", 3)
                    )
DNAmethContexts <- c(
                     "CpG",
                     "CpG",
                     "CHG",
                     "CHH"
                    )
DNAmethNamesPlot <- c(
                      "mCG (Nanopolish)",
                      "mCG (SE BS-seq)",
                      "mCHG (SE BS-seq)",
                      "mCHH (SE BS-seq)"
                     )
DNAmethColours <- c(
                    rep("red", length(DNAmethNamesPlot))
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
                                "_MappedOn_Athaliana_ONT_RaGOO_v2.0_", DNAmethContexts[x], "_CEN180_in_",
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

# gap
DNAmeth_gapMats <- mclapply(seq_along(DNAmethContexts), function(x) {
  lapply(seq_along(chrName), function(y) {
    as.matrix(read.table(paste0(DNAmethDirs[x],
                                DNAmethNames[x],
                                "_MappedOn_Athaliana_ONT_RaGOO_v2.0_", DNAmethContexts[x], "_CENgap_in_",
                                chrName[y], "_matrix_bin", binName, "_flank", flankName, ".tab"),
                         header = F, skip = 3))
  })
}, mc.cores = length(DNAmethContexts))
# If gaps from multiple chromosomes are to be analysed,
# concatenate the corresponding gap coverage matrices
DNAmeth_gapMats <- mclapply(seq_along(DNAmeth_gapMats), function(x) {
  if(length(chrName) > 1) {
    do.call(rbind, DNAmeth_gapMats[[x]])
  } else {
    DNAmeth_gapMats[[x]][[1]]
  }
}, mc.cores = length(DNAmeth_gapMats))

# Athila
DNAmeth_AthilaMats <- mclapply(seq_along(DNAmethContexts), function(x) {
  lapply(seq_along(chrName), function(y) {
    as.matrix(read.table(paste0(DNAmethDirs[x],
                                DNAmethNames[x],
                                "_MappedOn_Athaliana_ONT_RaGOO_v2.0_", DNAmethContexts[x], "_CENAthila_in_",
                                chrName[y], "_matrix_bin", binName, "_flank", flankName, ".tab"),
                         header = F, skip = 3))
  })
}, mc.cores = length(DNAmethContexts))
# If Athilas from multiple chromosomes are to be analysed,
# concatenate the corresponding Athila coverage matrices
DNAmeth_AthilaMats <- mclapply(seq_along(DNAmeth_AthilaMats), function(x) {
  if(length(chrName) > 1) {
    do.call(rbind, DNAmeth_AthilaMats[[x]])
  } else {
    DNAmeth_AthilaMats[[x]][[1]]
  }
}, mc.cores = length(DNAmeth_AthilaMats))

# DNAmeth
# Add column names
for(x in seq_along(DNAmeth_featureMats)) {
  colnames(DNAmeth_featureMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                           paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                           paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))
  colnames(DNAmeth_gapMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                       paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                       paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))
  colnames(DNAmeth_AthilaMats[[x]]) <- c(paste0("u", 1:(upstream/Athila_binSize)),
                                          paste0("t", ((upstream/Athila_binSize)+1):((upstream+Athila_bodyLength)/Athila_binSize)),
                                          paste0("d", (((upstream+Athila_bodyLength)/Athila_binSize)+1):(((upstream+Athila_bodyLength)/Athila_binSize)+(downstream/Athila_binSize))))
}

# Create list of lists in which each element in the enclosing list corresponds to a library
# and the two elements in the nested list correspond to coverage matrices for features and random loci 
DNAmeth_mats <- mclapply(seq_along(DNAmeth_featureMats), function(x) {
  list(
       # features
       DNAmeth_featureMats[[x]],
       # gaps
       DNAmeth_gapMats[[x]],
       # ranLocs
       DNAmeth_AthilaMats[[x]]
      ) 
}, mc.cores = length(DNAmeth_featureMats))

# Transpose matrix and convert into dataframe
# in which first column is window name
wideDFfeature_list_DNAmeth <- mclapply(seq_along(DNAmeth_mats), function(x) {
  lapply(seq_along(DNAmeth_mats[[x]]), function(y) {
    data.frame(window = colnames(DNAmeth_mats[[x]][[y]]),
               t(DNAmeth_mats[[x]][[y]]))
  })
}, mc.cores = length(DNAmeth_mats))

# Convert into tidy data.frame (long format)
tidyDFfeature_list_DNAmeth  <- mclapply(seq_along(wideDFfeature_list_DNAmeth), function(x) {
  lapply(seq_along(DNAmeth_mats[[x]]), function(y) {
    gather(data  = wideDFfeature_list_DNAmeth[[x]][[y]],
           key   = feature,
           value = coverage,
           -window)
  }) 
}, mc.cores = length(wideDFfeature_list_DNAmeth))

# Order levels of factor "window" so that sequential levels
# correspond to sequential windows
for(x in seq_along(tidyDFfeature_list_DNAmeth)) {
  for(y in seq_along(DNAmeth_mats[[x]])) {
    tidyDFfeature_list_DNAmeth[[x]][[y]]$window <- factor(tidyDFfeature_list_DNAmeth[[x]][[y]]$window,
                                                           levels = as.character(wideDFfeature_list_DNAmeth[[x]][[y]]$window))
  }
}

# Create summary data.frame in which each row corresponds to a window (Column 1),
# Column2 is the number of coverage values (features) per window,
# Column3 is the mean of coverage values per window,
# Column4 is the standard deviation of coverage values per window,
# Column5 is the standard error of the mean of coverage values per window,
# Column6 is the lower bound of the 95% confidence interval, and
# Column7 is the upper bound of the 95% confidence interval
summaryDFfeature_list_DNAmeth  <- mclapply(seq_along(tidyDFfeature_list_DNAmeth), function(x) {
  lapply(seq_along(DNAmeth_mats[[x]]), function(y) {
    data.frame(window = as.character(wideDFfeature_list_DNAmeth[[x]][[y]]$window),
               n      = tapply(X     = tidyDFfeature_list_DNAmeth[[x]][[y]]$coverage,
                               INDEX = tidyDFfeature_list_DNAmeth[[x]][[y]]$window,
                               FUN   = length),
               mean   = tapply(X     = tidyDFfeature_list_DNAmeth[[x]][[y]]$coverage,
                               INDEX = tidyDFfeature_list_DNAmeth[[x]][[y]]$window,
                               FUN   = mean,
                               na.rm = TRUE),
               sd     = tapply(X     = tidyDFfeature_list_DNAmeth[[x]][[y]]$coverage,
                               INDEX = tidyDFfeature_list_DNAmeth[[x]][[y]]$window,
                               FUN   = sd,
                               na.rm = TRUE))
  })
}, mc.cores = length(tidyDFfeature_list_DNAmeth))

for(x in seq_along(summaryDFfeature_list_DNAmeth)) {
  for(y in seq_along(DNAmeth_mats[[x]])) {
    summaryDFfeature_list_DNAmeth[[x]][[y]]$window <- factor(summaryDFfeature_list_DNAmeth[[x]][[y]]$window,
                                                              levels = as.character(wideDFfeature_list_DNAmeth[[x]][[y]]$window))
    summaryDFfeature_list_DNAmeth[[x]][[y]]$winNo <- factor(1:dim(summaryDFfeature_list_DNAmeth[[x]][[y]])[1])
    summaryDFfeature_list_DNAmeth[[x]][[y]]$sem <- summaryDFfeature_list_DNAmeth[[x]][[y]]$sd/sqrt(summaryDFfeature_list_DNAmeth[[x]][[y]]$n-1)
    summaryDFfeature_list_DNAmeth[[x]][[y]]$CI_lower <- summaryDFfeature_list_DNAmeth[[x]][[y]]$mean -
      qt(0.975, df = summaryDFfeature_list_DNAmeth[[x]][[y]]$n-1)*summaryDFfeature_list_DNAmeth[[x]][[y]]$sem
    summaryDFfeature_list_DNAmeth[[x]][[y]]$CI_upper <- summaryDFfeature_list_DNAmeth[[x]][[y]]$mean +
      qt(0.975, df = summaryDFfeature_list_DNAmeth[[x]][[y]]$n-1)*summaryDFfeature_list_DNAmeth[[x]][[y]]$sem
  }
}

summaryDFfeature_DNAmeth <- summaryDFfeature_list_DNAmeth

# Define y-axis limits
ymin_list_DNAmeth <- lapply(seq_along(summaryDFfeature_DNAmeth), function(x) {
  min(c(summaryDFfeature_DNAmeth[[x]][[1]]$CI_lower,
        summaryDFfeature_DNAmeth[[x]][[2]]$CI_lower,
        summaryDFfeature_DNAmeth[[x]][[3]]$CI_lower))
})
ymax_list_DNAmeth <- lapply(seq_along(summaryDFfeature_DNAmeth), function(x) {
  max(c(summaryDFfeature_DNAmeth[[x]][[1]]$CI_upper,
        summaryDFfeature_DNAmeth[[x]][[2]]$CI_upper,
        summaryDFfeature_DNAmeth[[x]][[3]]$CI_upper))
})

# Plot average profiles with 95% CI ribbon
## feature
ggObj1_combined_DNAmeth <- mclapply(seq_along(DNAmethNamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_DNAmeth[[x]][[1]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = 1)
        ) +
  geom_line(data = summaryDFfeature,
            colour = DNAmethColours[x],
            size = 1) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper),
              fill = DNAmethColours[x],
              alpha = 0.4) +
  scale_y_continuous(limits = c(ymin_list_DNAmeth[[x]], ymax_list_DNAmeth[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_DNAmeth[[x]][[1]])[1])-(downstream/binSize),
                              dim(summaryDFfeature_DNAmeth[[x]][[1]])[1]),
                   labels = c(paste0("-", flankName),
                              featureStartLab,
                              featureEndLab,
                              paste0("+", flankName))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_DNAmeth[[x]][[1]])[1])-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = DNAmethNamesPlot[x]) +
  theme_bw() +
  theme(
        axis.ticks = element_line(size = 1.0, colour = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.x = element_text(size = 22, colour = "black"),
        axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
        axis.title = element_text(size = 30, colour = DNAmethColours[x]),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_rect(size = 3.5, colour = "black"),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
        plot.title = element_text(hjust = 0.5, size = 30)) +
  ggtitle(bquote(.(featureNamePlot) ~ "(" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFfeature$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
}, mc.cores = length(DNAmethNamesPlot))

## gap
ggObj2_combined_DNAmeth <- mclapply(seq_along(DNAmethNamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_DNAmeth[[x]][[2]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = 1)
        ) +
  geom_line(data = summaryDFfeature,
            colour = DNAmethColours[x],
            size = 1) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper),
              fill = DNAmethColours[x],
              alpha = 0.4) +
  scale_y_continuous(limits = c(ymin_list_DNAmeth[[x]], ymax_list_DNAmeth[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_DNAmeth[[x]][[2]])[1])-(downstream/binSize),
                              dim(summaryDFfeature_DNAmeth[[x]][[2]])[1]),
                   labels = c(paste0("-", flankName),
                              featureStartLab,
                              featureEndLab,
                              paste0("+", flankName))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_DNAmeth[[x]][[2]])[1])-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = DNAmethNamesPlot[x]) +
  theme_bw() +
  theme(
        axis.ticks = element_line(size = 1.0, colour = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.x = element_text(size = 22, colour = "black"),
        axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
        axis.title = element_text(size = 30, colour = DNAmethColours[x]),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_rect(size = 3.5, colour = "black"),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
        plot.title = element_text(hjust = 0.5, size = 30)) +
  ggtitle(bquote(.(gapNamePlot) ~ "(" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFfeature$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
}, mc.cores = length(DNAmethNamesPlot))

## Athila
ggObj3_combined_DNAmeth <- mclapply(seq_along(DNAmethNamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_DNAmeth[[x]][[3]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = 1)
        ) +
  geom_line(data = summaryDFfeature,
            colour = DNAmethColours[x],
            size = 1) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper),
              fill = DNAmethColours[x],
              alpha = 0.4) +
  scale_y_continuous(limits = c(ymin_list_DNAmeth[[x]], ymax_list_DNAmeth[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/Athila_binSize)+1,
                              (dim(summaryDFfeature_DNAmeth[[x]][[3]])[1])-(downstream/Athila_binSize),
                              dim(summaryDFfeature_DNAmeth[[x]][[3]])[1]),
                   labels = c(paste0("-", flankName),
                              featureStartLab,
                              featureEndLab,
                              paste0("+", flankName))) +
  geom_vline(xintercept = c((upstream/Athila_binSize)+1,
                            (dim(summaryDFfeature_DNAmeth[[x]][[3]])[1])-(downstream/Athila_binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = DNAmethNamesPlot[x]) +
  theme_bw() +
  theme(
        axis.ticks = element_line(size = 1.0, colour = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.x = element_text(size = 22, colour = "black"),
        axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
        axis.title = element_text(size = 30, colour = DNAmethColours[x]),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_rect(size = 3.5, colour = "black"),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
        plot.title = element_text(hjust = 0.5, size = 30)) +
  ggtitle(bquote(.(AthilaNamePlot) ~ "(" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFfeature$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
}, mc.cores = length(DNAmethNamesPlot))

ggObjGA_combined <- grid.arrange(grobs = c(
                                           ggObj1_combined_DNAmeth,
                                           ggObj2_combined_DNAmeth,
                                           ggObj3_combined_DNAmeth
                                          ),
                                 layout_matrix = cbind(
                                                       1:length(c(DNAmethNamesPlot)),
                                                       (length(c(DNAmethNamesPlot))+1):(length(c(DNAmethNamesPlot))*2),
                                                       ((length(c(DNAmethNamesPlot))*2)+1):(length(c(DNAmethNamesPlot))*3)
                                                      ))
ggsave(paste0(plotDir,
              "DNAmeth_avgProfiles_around",
               "_CEN180_CENgap_CENAthila_in_RaGOO_v2.0_",
               paste0(chrName, collapse = "_"), "_", align, ".pdf"),
       plot = ggObjGA_combined,
       height = 6.5*length(c(DNAmethNamesPlot)), width = 30, limitsize = FALSE)

#### Free up memory by removing no longer required objects
rm(
   DNAmeth_featureMats, DNAmeth_AthilaMats,
   DNAmeth_mats,
   wideDFfeature_list_DNAmeth,
   tidyDFfeature_list_DNAmeth,
   summaryDFfeature_list_DNAmeth
#   summaryDFfeature_DNAmeth
  ) 
gc()
#####
