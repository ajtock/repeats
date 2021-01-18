#!/applications/R/R-4.0.0/bin/Rscript

# author: Andy Tock
# contact: ajt200@cam.ac.uk
# date: 18.12.2020

# Calculate and plot metaprofiles of ChIP-seq, MNase-seq, etc.
# (CEN180 windowed means and 95% confidence intervals, CIs)
# for all CEN180 sequences and randomly positioned loci

# Usage:
# /applications/R/R-4.0.0/bin/Rscript CEN180_CENAthila_CENgap_TEsuperfam_genes_6metaprofiles.R 'Chr1,Chr2,Chr3,Chr4,Chr5' both 180 11000 1000 2200 1000 1kb 10 200 10bp 10bp '0.02,0.40' 'WT_SPO11oligos_Rep1,WT_SPO11oligos_Rep2,WT_SPO11oligos_Rep3,met1_SPO11oligos_Rep1,met1_SPO11oligos_Rep2,met1_SPO11oligos_Rep3' '160518_Kyuha_SPO11oligos/WT/snakemake_SPO11oligos_T2T_Col,160518_Kyuha_SPO11oligos/WT/snakemake_SPO11oligos_T2T_Col,160518_Kyuha_SPO11oligos/WT/snakemake_SPO11oligos_T2T_Col,160518_Kyuha_SPO11oligos/met1/snakemake_SPO11oligos_T2T_Col,160518_Kyuha_SPO11oligos/met1/snakemake_SPO11oligos_T2T_Col,160518_Kyuha_SPO11oligos/met1/snakemake_SPO11oligos_T2T_Col' 'wt Rep1,wt Rep2,wt Rep3,met1 Rep1,met1 Rep2,met1 Rep3' 'grey0,grey20,grey40,magenta,purple,purple4' 'SPO11-1-oligos'

#chrName <- unlist(strsplit("Chr1,Chr2,Chr3,Chr4,Chr5",
#                           split = ","))
#align <- "both"
#bodyLength <- 180
#Athila_bodyLength <- 11000
#TEsf_bodyLength <- 1000
#gene_bodyLength <- 2200
#upstream <- 1000
#downstream <- 1000
#flankName <- "1kb"
#binSize <- 10
#Athila_binSize <- 200
#binName <- "10bp"
#Athila_binName <- "10bp"
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
#legendPos <- as.numeric(unlist(strsplit("0.02,0.40",
#                                        split = ",")))
#ChIPNames <- unlist(strsplit("WT_SPO11oligos_Rep1,WT_SPO11oligos_Rep2,WT_SPO11oligos_Rep3,met1_SPO11oligos_Rep1,met1_SPO11oligos_Rep2,met1_SPO11oligos_Rep3",
#                             split = ","))
#ChIPNamesDir <- unlist(strsplit("160518_Kyuha_SPO11oligos/WT/snakemake_SPO11oligos_T2T_Col,160518_Kyuha_SPO11oligos/WT/snakemake_SPO11oligos_T2T_Col,160518_Kyuha_SPO11oligos/WT/snakemake_SPO11oligos_T2T_Col,160518_Kyuha_SPO11oligos/met1/snakemake_SPO11oligos_T2T_Col,160518_Kyuha_SPO11oligos/met1/snakemake_SPO11oligos_T2T_Col,160518_Kyuha_SPO11oligos/met1/snakemake_SPO11oligos_T2T_Col",
#                                split = ","))
#log2ChIPNamesPlot <- unlist(strsplit("wt Rep1,wt Rep2,wt Rep3,met1 Rep1,met1 Rep2,met1 Rep3",
#                                     split = ","))
#log2ChIPColours <- unlist(strsplit("deepskyblue,blue,navy,tomato,red,red4",
#                                   split = ","))
#yLabPlot <- "SPO11-1-oligos"

args <- commandArgs(trailingOnly = T)
chrName <- unlist(strsplit(args[1],
                           split = ","))
align <- args[2]
bodyLength <- as.numeric(args[3])
Athila_bodyLength <- as.numeric(args[4])
TEsf_bodyLength <- as.numeric(args[5])
gene_bodyLength <- as.numeric(args[6])
upstream <- as.numeric(args[7])
downstream <- as.numeric(args[7])
flankName <- args[8]
binSize <- as.numeric(args[9])
Athila_binSize <- as.numeric(args[10])
binName <- args[11]
Athila_binName <- args[12]
legendPos <- as.numeric(unlist(strsplit(args[13],
                                        split = ",")))
ChIPNames <- unlist(strsplit(args[14],
                             split = ","))
ChIPNamesDir <- unlist(strsplit(args[15],
                                split = ","))
log2ChIPNamesPlot <- unlist(strsplit(args[16],
                                     split = ","))
log2ChIPColours <- unlist(strsplit(args[17],
                                   split = ","))
yLabPlot <- args[18]

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
  ranLocNamePlot <- "All ranLoc"
  gapNamePlot <- "All CENgap"
  AthilaNamePlot <- "All CENAthila"
  TEsfNamePlot <- "All Gypsy LTRs"
  geneNamePlot <- "All genes"
} else {
  featureNamePlot <- paste0(paste0(chrName, collapse = ","), " CEN180")
  ranLocNamePlot <- paste0(paste0(chrName, collapse = ","), " ranLoc")
  gapNamePlot <- paste0(paste0(chrName, collapse = ","), " CENgap")
  AthilaNamePlot <- paste0(paste0(chrName, collapse = ","), " CENAthila")
  TEsfNamePlot <- paste0(paste0(chrName, collapse = ","), " Gypsy LTRs")
  geneNamePlot <- paste0(paste0(chrName, collapse = ","), " genes")
}

# Define feature start and end labels for plotting
featureStartLab <- "Start"
featureEndLab <- "End"
geneStartLab <- "TSS"
geneEndLab <- "TTS"

# Load feature matrices for each chromatin dataset, calculate log2(ChIP/control),
ChIPNamesPlot <- log2ChIPNamesPlot
ChIPColours <- log2ChIPColours
# Load feature matrices for each chromatin dataset, calculate log2(ChIP/control),
ChIPDirs <- sapply(seq_along(ChIPNames), function(x) {
  paste0("/home/ajt200/analysis/",
         ChIPNamesDir[x],
         "/mapped/")
})

controlNames <- c(
#                  "WT_REC8_Myc_Rep1_input",
                  "WT_gDNA_Rep1_R1"
                 )
controlNamesDir <- c(
#                     "REC8_pooled/snakemake_ChIPseq_T2T_Col",
                     "150701_Natasha_gDNA/WT/R1/snakemake_SPO11oligos_T2T_Col"
                    )
controlNamesPlot <- c(
#                      "PE input",
                      "SE gDNA"
                     )
controlColours <- c(
                    rep("red", length(controlNamesPlot))
                   )
controlDirs <- sapply(seq_along(controlNames), function(x) {
  paste0("/home/ajt200/analysis/",
         controlNamesDir[x],
         "/mapped/")
})

## ChIP
# feature
ChIP_featureMats <- mclapply(seq_along(ChIPNames), function(x) {
  lapply(seq_along(chrName), function(y) {
    as.matrix(read.table(paste0(ChIPDirs[x], "CEN180profiles/matrices/",
                                ChIPNames[x],
                                "_MappedOn_T2T_Col_lowXM_", align, "_sort_norm_CEN180_in_",
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
    as.matrix(read.table(paste0(controlDirs[x], "CEN180profiles/matrices/",
                                controlNames[x],
                                "_MappedOn_T2T_Col_lowXM_", align, "_sort_norm_CEN180_in_",
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
    as.matrix(read.table(paste0(ChIPDirs[x], "CEN180profiles/matrices/",
                                ChIPNames[x],
                                "_MappedOn_T2T_Col_lowXM_", align, "_sort_norm_CEN180_in_",
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
    as.matrix(read.table(paste0(controlDirs[x], "CEN180profiles/matrices/",
                                controlNames[x],
                                "_MappedOn_T2T_Col_lowXM_", align, "_sort_norm_CEN180_in_",
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

## ChIP
# gap
ChIP_gapMats <- mclapply(seq_along(ChIPNames), function(x) {
  lapply(seq_along(chrName), function(y) {
    as.matrix(read.table(paste0(ChIPDirs[x], "CEN180profiles/matrices/",
                                ChIPNames[x],
                                "_MappedOn_T2T_Col_lowXM_", align, "_sort_norm_CENgap_in_",
                                chrName[y], "_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
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
    as.matrix(read.table(paste0(controlDirs[x], "CEN180profiles/matrices/",
                                controlNames[x],
                                "_MappedOn_T2T_Col_lowXM_", align, "_sort_norm_CENgap_in_",
                                chrName[y], "_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                         header = F, skip = 3))
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
    as.matrix(read.table(paste0(ChIPDirs[x], "CEN180profiles/matrices/",
                                ChIPNames[x],
                                "_MappedOn_T2T_Col_lowXM_", align, "_sort_norm_CENAthila_in_",
                                chrName[y], "_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
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
    as.matrix(read.table(paste0(controlDirs[x], "CEN180profiles/matrices/",
                                controlNames[x],
                                "_MappedOn_T2T_Col_lowXM_", align, "_sort_norm_CENAthila_in_",
                                chrName[y], "_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                         header = F, skip = 3))
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

## ChIP
# TEsf
ChIP_TEsfMats <- mclapply(seq_along(ChIPNames), function(x) {
  lapply(seq_along(chrName), function(y) {
    as.matrix(read.table(paste0(ChIPDirs[x], "TEprofiles/matrices/",
                                ChIPNames[x],
                                "_MappedOn_T2T_Col_lowXM_", align, "_sort_norm_TEs_Gypsy_LTR_in_",
                                chrName[y], "_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                         header = F, skip = 3))
  })
}, mc.cores = length(ChIPNames))
# If TEsfs from multiple chromosomes are to be analysed,
# concatenate the corresponding TEsf coverage matrices
ChIP_TEsfMats <- mclapply(seq_along(ChIP_TEsfMats), function(x) {
  if(length(chrName) > 1) {
    do.call(rbind, ChIP_TEsfMats[[x]])
  } else {
    ChIP_TEsfMats[[x]][[1]]
  }
}, mc.cores = length(ChIP_TEsfMats))

## control
# TEsf
control_TEsfMats <- mclapply(seq_along(controlNames), function(x) {
  lapply(seq_along(chrName), function(y) {
    as.matrix(read.table(paste0(controlDirs[x], "TEprofiles/matrices/",
                                controlNames[x],
                                "_MappedOn_T2T_Col_lowXM_", align, "_sort_norm_TEs_Gypsy_LTR_in_",
                                chrName[y], "_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                         header = F, skip = 3))
  })
}, mc.cores = length(controlNames))
# If TEsfs from multiple chromosomes are to be analysed,
# concatenate the corresponding TEsf coverage matrices
control_TEsfMats <- mclapply(seq_along(control_TEsfMats), function(x) {
  if(length(chrName) > 1) {
    do.call(rbind, control_TEsfMats[[x]])
  } else {
    control_TEsfMats[[x]][[1]]
  }
}, mc.cores = length(control_TEsfMats))

## ChIP
# gene
ChIP_geneMats <- mclapply(seq_along(ChIPNames), function(x) {
  lapply(seq_along(chrName), function(y) {
    as.matrix(read.table(paste0(ChIPDirs[x], "geneProfiles/matrices/",
                                ChIPNames[x],
                                "_MappedOn_T2T_Col_lowXM_", align, "_sort_norm_genes_in_",
                                chrName[y], "_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                         header = F, skip = 3))
  })
}, mc.cores = length(ChIPNames))
# If genes from multiple chromosomes are to be analysed,
# concatenate the corresponding gene coverage matrices
ChIP_geneMats <- mclapply(seq_along(ChIP_geneMats), function(x) {
  if(length(chrName) > 1) {
    do.call(rbind, ChIP_geneMats[[x]])
  } else {
    ChIP_geneMats[[x]][[1]]
  }
}, mc.cores = length(ChIP_geneMats))

## control
# gene
control_geneMats <- mclapply(seq_along(controlNames), function(x) {
  lapply(seq_along(chrName), function(y) {
    as.matrix(read.table(paste0(controlDirs[x], "geneProfiles/matrices/",
                                controlNames[x],
                                "_MappedOn_T2T_Col_lowXM_", align, "_sort_norm_genes_in_",
                                chrName[y], "_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                         header = F, skip = 3))
  })
}, mc.cores = length(controlNames))
# If genes from multiple chromosomes are to be analysed,
# concatenate the corresponding gene coverage matrices
control_geneMats <- mclapply(seq_along(control_geneMats), function(x) {
  if(length(chrName) > 1) {
    do.call(rbind, control_geneMats[[x]])
  } else {
    control_geneMats[[x]][[1]]
  }
}, mc.cores = length(control_geneMats))

## Conditionally calculate log2(ChIP/control)
## for each matrix depending on library
## feature
#log2ChIP_featureMats <- mclapply(seq_along(ChIP_featureMats), function(x) {
#  if ( grepl("MNase", ChIPNames[x]) ) {
#    print(paste0(ChIPNames[x], " library; using ", controlNames[1], " for log2((MNase+1)/(gDNA+1)) calculation"))
#    log2((ChIP_featureMats[[x]]+1)/(control_featureMats[[1]]+1))
#  } else if ( grepl("SPO11oligos", ChIPNames[x]) ) {
#    print(paste0(ChIPNames[x], " library; using ", controlNames[2], " for log2((SPO11-1-oligos+1)/(gDNA+1)) calculation"))
#    log2((ChIP_featureMats[[x]]+1)/(control_featureMats[[2]]+1))
#  } else {
#    print(paste0(ChIPNames[x], " library; using ", controlNames[1], " for log2((ChIP+1)/(input+1)) calculation"))
#    log2((ChIP_featureMats[[x]]+1)/(control_featureMats[[1]]+1))
#  }
#}, mc.cores = length(ChIP_featureMats))
#
## Conditionally calculate log2(ChIP/control)
## for each matrix depending on library
## ranLoc
#log2ChIP_ranLocMats <- mclapply(seq_along(ChIP_ranLocMats), function(x) {
#  if ( grepl("MNase", ChIPNames[x]) ) {
#    print(paste0(ChIPNames[x], " library; using ", controlNames[1], " for log2((MNase+1)/(gDNA+1)) calculation"))
#    log2((ChIP_ranLocMats[[x]]+1)/(control_ranLocMats[[1]]+1))
#  } else if ( grepl("SPO11oligos", ChIPNames[x]) ) {
#    print(paste0(ChIPNames[x], " library; using ", controlNames[2], " for log2((SPO11-1-oligos+1)/(gDNA+1)) calculation"))
#    log2((ChIP_ranLocMats[[x]]+1)/(control_ranLocMats[[2]]+1))
#  } else {
#    print(paste0(ChIPNames[x], " library; using ", controlNames[1], " for log2((ChIP+1)/(input+1)) calculation"))
#    log2((ChIP_ranLocMats[[x]]+1)/(control_ranLocMats[[1]]+1))
#  }
#}, mc.cores = length(ChIP_ranLocMats))
#
## Conditionally calculate log2(ChIP/control)
## for each matrix depending on library
## gap
#log2ChIP_gapMats <- mclapply(seq_along(ChIP_gapMats), function(x) {
#  if ( grepl("MNase", ChIPNames[x]) ) {
#    print(paste0(ChIPNames[x], " library; using ", controlNames[1], " for log2((MNase+1)/(gDNA+1)) calculation"))
#    log2((ChIP_gapMats[[x]]+1)/(control_gapMats[[1]]+1))
#  } else if ( grepl("SPO11oligos", ChIPNames[x]) ) {
#    print(paste0(ChIPNames[x], " library; using ", controlNames[2], " for log2((SPO11-1-oligos+1)/(gDNA+1)) calculation"))
#    log2((ChIP_gapMats[[x]]+1)/(control_gapMats[[2]]+1))
#  } else {
#    print(paste0(ChIPNames[x], " library; using ", controlNames[1], " for log2((ChIP+1)/(input+1)) calculation"))
#    log2((ChIP_gapMats[[x]]+1)/(control_gapMats[[1]]+1))
#  }
#}, mc.cores = length(ChIP_gapMats))
#
## Conditionally calculate log2(ChIP/control)
## for each matrix depending on library
## Athila
#log2ChIP_AthilaMats <- mclapply(seq_along(ChIP_AthilaMats), function(x) {
#  if ( grepl("MNase", ChIPNames[x]) ) {
#    print(paste0(ChIPNames[x], " library; using ", controlNames[1], " for log2((MNase+1)/(gDNA+1)) calculation"))
#    log2((ChIP_AthilaMats[[x]]+1)/(control_AthilaMats[[1]]+1))
#  } else if ( grepl("SPO11oligos", ChIPNames[x]) ) {
#    print(paste0(ChIPNames[x], " library; using ", controlNames[2], " for log2((SPO11-1-oligos+1)/(gDNA+1)) calculation"))
#    log2((ChIP_AthilaMats[[x]]+1)/(control_AthilaMats[[2]]+1))
#  } else {
#    print(paste0(ChIPNames[x], " library; using ", controlNames[1], " for log2((ChIP+1)/(input+1)) calculation"))
#    log2((ChIP_AthilaMats[[x]]+1)/(control_AthilaMats[[1]]+1))
#  }
#}, mc.cores = length(ChIP_AthilaMats))
#
#
## log2ChIP
## Add column names
#for(x in seq_along(log2ChIP_featureMats)) {
#  colnames(log2ChIP_featureMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
#                                           paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
#                                           paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))
#  colnames(log2ChIP_ranLocMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
#                                          paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
#                                          paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))
#  colnames(log2ChIP_gapMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
#                                       paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
#                                       paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))
#  colnames(log2ChIP_AthilaMats[[x]]) <- c(paste0("u", 1:(upstream/Athila_binSize)),
#                                           paste0("t", ((upstream/Athila_binSize)+1):((upstream+Athila_bodyLength)/Athila_binSize)),
#                                           paste0("d", (((upstream+Athila_bodyLength)/Athila_binSize)+1):(((upstream+Athila_bodyLength)/Athila_binSize)+(downstream/Athila_binSize))))
#}
#
## Create list of lists in which each element in the enclosing list corresponds to a library
## and the two elements in the nested list correspond to coverage matrices for features and random loci 
#log2ChIP_mats <- mclapply(seq_along(log2ChIP_featureMats), function(x) {
#  list(
#       # features
#       log2ChIP_featureMats[[x]],
#       # ranLocs
#       log2ChIP_ranLocMats[[x]],
#       # gaps
#       log2ChIP_gapMats[[x]],
#       # Athilas
#       log2ChIP_AthilaMats[[x]]
#      ) 
#}, mc.cores = length(log2ChIP_featureMats))
#
## Transpose matrix and convert into dataframe
## in which first column is window name
#wideDFfeature_list_log2ChIP <- mclapply(seq_along(log2ChIP_mats), function(x) {
#  lapply(seq_along(log2ChIP_mats[[x]]), function(y) {
#    data.frame(window = colnames(log2ChIP_mats[[x]][[y]]),
#               t(log2ChIP_mats[[x]][[y]]))
#  })
#}, mc.cores = length(log2ChIP_mats))
#
## Convert into tidy data.frame (long format)
#tidyDFfeature_list_log2ChIP  <- mclapply(seq_along(wideDFfeature_list_log2ChIP), function(x) {
#  lapply(seq_along(log2ChIP_mats[[x]]), function(y) {
#    gather(data  = wideDFfeature_list_log2ChIP[[x]][[y]],
#           key   = feature,
#           value = coverage,
#           -window)
#  }) 
#}, mc.cores = length(wideDFfeature_list_log2ChIP))
#
## Order levels of factor "window" so that sequential levels
## correspond to sequential windows
#for(x in seq_along(tidyDFfeature_list_log2ChIP)) {
#  for(y in seq_along(log2ChIP_mats[[x]])) {
#    tidyDFfeature_list_log2ChIP[[x]][[y]]$window <- factor(tidyDFfeature_list_log2ChIP[[x]][[y]]$window,
#                                                           levels = as.character(wideDFfeature_list_log2ChIP[[x]][[y]]$window))
#  }
#}
#
## Create summary data.frame in which each row corresponds to a window (Column 1),
## Column2 is the number of coverage values (features) per window,
## Column3 is the mean of coverage values per window,
## Column4 is the standard deviation of coverage values per window,
## Column5 is the standard error of the mean of coverage values per window,
## Column6 is the lower bound of the 95% confidence interval, and
## Column7 is the upper bound of the 95% confidence interval
#summaryDFfeature_list_log2ChIP  <- mclapply(seq_along(tidyDFfeature_list_log2ChIP), function(x) {
#  lapply(seq_along(log2ChIP_mats[[x]]), function(y) {
#    data.frame(window = as.character(wideDFfeature_list_log2ChIP[[x]][[y]]$window),
#               n      = tapply(X     = tidyDFfeature_list_log2ChIP[[x]][[y]]$coverage,
#                               INDEX = tidyDFfeature_list_log2ChIP[[x]][[y]]$window,
#                               FUN   = length),
#               mean   = tapply(X     = tidyDFfeature_list_log2ChIP[[x]][[y]]$coverage,
#                               INDEX = tidyDFfeature_list_log2ChIP[[x]][[y]]$window,
#                               FUN   = mean,
#                               na.rm = TRUE),
#               sd     = tapply(X     = tidyDFfeature_list_log2ChIP[[x]][[y]]$coverage,
#                               INDEX = tidyDFfeature_list_log2ChIP[[x]][[y]]$window,
#                               FUN   = sd,
#                               na.rm = TRUE))
#  })
#}, mc.cores = length(tidyDFfeature_list_log2ChIP))
#
#for(x in seq_along(summaryDFfeature_list_log2ChIP)) {
#  for(y in seq_along(log2ChIP_mats[[x]])) {
#    summaryDFfeature_list_log2ChIP[[x]][[y]]$window <- factor(summaryDFfeature_list_log2ChIP[[x]][[y]]$window,
#                                                              levels = as.character(wideDFfeature_list_log2ChIP[[x]][[y]]$window))
#    summaryDFfeature_list_log2ChIP[[x]][[y]]$winNo <- factor(1:dim(summaryDFfeature_list_log2ChIP[[x]][[y]])[1])
#    summaryDFfeature_list_log2ChIP[[x]][[y]]$sem <- summaryDFfeature_list_log2ChIP[[x]][[y]]$sd/sqrt(summaryDFfeature_list_log2ChIP[[x]][[y]]$n-1)
#    summaryDFfeature_list_log2ChIP[[x]][[y]]$CI_lower <- summaryDFfeature_list_log2ChIP[[x]][[y]]$mean -
#      qt(0.975, df = summaryDFfeature_list_log2ChIP[[x]][[y]]$n-1)*summaryDFfeature_list_log2ChIP[[x]][[y]]$sem
#    summaryDFfeature_list_log2ChIP[[x]][[y]]$CI_upper <- summaryDFfeature_list_log2ChIP[[x]][[y]]$mean +
#      qt(0.975, df = summaryDFfeature_list_log2ChIP[[x]][[y]]$n-1)*summaryDFfeature_list_log2ChIP[[x]][[y]]$sem
#  }
#}
#
## Convert list of lists summaryDFfeature_list_log2ChIP into
## a list of single data.frames containing all meta-profiles for plotting
#featureTmp <- lapply(seq_along(summaryDFfeature_list_log2ChIP), function(x) {
#  summaryDFfeature_list_log2ChIP[[x]][[1]]
#})
#ranLocTmp <- lapply(seq_along(summaryDFfeature_list_log2ChIP), function(x) {
#  summaryDFfeature_list_log2ChIP[[x]][[2]]
#})
#gapTmp <- lapply(seq_along(summaryDFfeature_list_log2ChIP), function(x) {
#  summaryDFfeature_list_log2ChIP[[x]][[3]]
#})
#AthilaTmp <- lapply(seq_along(summaryDFfeature_list_log2ChIP), function(x) {
#  summaryDFfeature_list_log2ChIP[[x]][[4]]
#})
#names(featureTmp) <- log2ChIPNamesPlot
#names(ranLocTmp) <- log2ChIPNamesPlot
#names(gapTmp) <- log2ChIPNamesPlot
#names(AthilaTmp) <- log2ChIPNamesPlot
#summaryDFfeature_log2ChIP <- list(
#  bind_rows(featureTmp, .id = "libName"),
#  bind_rows(ranLocTmp, .id = "libName"),
#  bind_rows(gapTmp, .id = "libName"),
#  bind_rows(AthilaTmp, .id = "libName")
#)  
#for(x in seq_along(summaryDFfeature_log2ChIP)) {
#  summaryDFfeature_log2ChIP[[x]]$libName <- factor(summaryDFfeature_log2ChIP[[x]]$libName,
#                                                   levels = log2ChIPNamesPlot)
#}
#
## Define y-axis limits
#ymin_log2ChIP <- min(c(summaryDFfeature_log2ChIP[[1]]$CI_lower,
#                       summaryDFfeature_log2ChIP[[2]]$CI_lower,
#                       summaryDFfeature_log2ChIP[[3]]$CI_lower,
#                       summaryDFfeature_log2ChIP[[4]]$CI_lower))
#ymax_log2ChIP <- max(c(summaryDFfeature_log2ChIP[[1]]$CI_upper,
#                       summaryDFfeature_log2ChIP[[2]]$CI_upper,
#                       summaryDFfeature_log2ChIP[[3]]$CI_upper,
#                       summaryDFfeature_log2ChIP[[4]]$CI_upper))
#
## Define legend labels
#legendLabs <- lapply(seq_along(log2ChIPNamesPlot), function(x) {
#  grobTree(textGrob(bquote(.(log2ChIPNamesPlot[x])),
#                    x = legendPos[1], y = legendPos[2]-((x-1)*0.06), just = "left",
#                    gp = gpar(col = log2ChIPColours[x], fontsize = 18)))
#})
#
## Plot average profiles with 95% CI ribbon
### feature
#summaryDFfeature <- summaryDFfeature_log2ChIP[[1]]
#ggObj1_combined_log2ChIP <- ggplot(data = summaryDFfeature,
#                                   mapping = aes(x = winNo,
#                                                 y = mean,
#                                                 group = libName)
#                                  ) +
#geom_line(data = summaryDFfeature,
#          mapping = aes(colour = libName),
#          size = 1) +
#scale_colour_manual(values = log2ChIPColours) +
#geom_ribbon(data = summaryDFfeature,
#            mapping = aes(ymin = CI_lower,
#                          ymax = CI_upper,
#                          fill = libName),
#            alpha = 0.4) +
#scale_fill_manual(values = log2ChIPColours) +
#scale_y_continuous(limits = c(ymin_log2ChIP, ymax_log2ChIP),
#                   labels = function(x) sprintf("%6.3f", x)) +
#scale_x_discrete(breaks = c(1,
#                            (upstream/binSize)+1,
#                            (dim(summaryDFfeature_log2ChIP[[1]])[1]/length(ChIPNames))-(downstream/binSize),
#                            dim(summaryDFfeature_log2ChIP[[1]])[1]/length(ChIPNames)),
#                 labels = c(paste0("-", flankName),
#                            featureStartLab,
#                            featureEndLab,
#                            paste0("+", flankName))) +
#geom_vline(xintercept = c((upstream/binSize)+1,
#                          (dim(summaryDFfeature_log2ChIP[[1]])[1]/length(ChIPNames))-(downstream/binSize)),
#           linetype = "dashed",
#           size = 1) +
#labs(x = "",
#     y = bquote("Log"[2] * "(" * .(yLabPlot) * "/control)")) +
#theme_bw() +
#theme(
#      axis.ticks = element_line(size = 1.0, colour = "black"),
#      axis.ticks.length = unit(0.25, "cm"),
#      axis.text.x = element_text(size = 22, colour = "black"),
#      axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
#      axis.title = element_text(size = 30, colour = "black"),
#      legend.position = "none",
#      panel.grid = element_blank(),
#      panel.border = element_rect(size = 3.5, colour = "black"),
#      panel.background = element_blank(),
#      plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
#      plot.title = element_text(hjust = 0.5, size = 30)) +
#ggtitle(bquote(.(featureNamePlot) ~ "(" * italic("n") ~ "=" ~
#               .(prettyNum(summaryDFfeature$n[1],
#                           big.mark = ",", trim = T)) *
#               ")"))
#
### ranLoc
#summaryDFfeature <- summaryDFfeature_log2ChIP[[2]]
#ggObj2_combined_log2ChIP <- ggplot(data = summaryDFfeature,
#                                   mapping = aes(x = winNo,
#                                                 y = mean,
#                                                 group = libName)
#                                  ) +
#geom_line(data = summaryDFfeature,
#          mapping = aes(colour = libName),
#          size = 1) +
#scale_colour_manual(values = log2ChIPColours) +
#geom_ribbon(data = summaryDFfeature,
#            mapping = aes(ymin = CI_lower,
#                          ymax = CI_upper,
#                          fill = libName),
#            alpha = 0.4) +
#scale_fill_manual(values = log2ChIPColours) +
#scale_y_continuous(limits = c(ymin_log2ChIP, ymax_log2ChIP),
#                   labels = function(x) sprintf("%6.3f", x)) +
#scale_x_discrete(breaks = c(1,
#                            (upstream/binSize)+1,
#                            (dim(summaryDFfeature_log2ChIP[[2]])[1]/length(ChIPNames))-(downstream/binSize),
#                            dim(summaryDFfeature_log2ChIP[[2]])[1]/length(ChIPNames)),
#                 labels = c(paste0("-", flankName),
#                            featureStartLab,
#                            featureEndLab,
#                            paste0("+", flankName))) +
#geom_vline(xintercept = c((upstream/binSize)+1,
#                          (dim(summaryDFfeature_log2ChIP[[2]])[1]/length(ChIPNames))-(downstream/binSize)),
#           linetype = "dashed",
#           size = 1) +
#labs(x = "",
#     y = bquote("Log"[2] * "(" * .(yLabPlot) * "/control)")) +
#annotation_custom(legendLabs[[1]]) +
#annotation_custom(legendLabs[[2]]) +
#annotation_custom(legendLabs[[3]]) +
#annotation_custom(legendLabs[[4]]) +
#annotation_custom(legendLabs[[5]]) +
#annotation_custom(legendLabs[[6]]) +
#theme_bw() +
#theme(
#      axis.ticks = element_line(size = 1.0, colour = "black"),
#      axis.ticks.length = unit(0.25, "cm"),
#      axis.text.x = element_text(size = 22, colour = "black"),
#      axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
#      axis.title = element_text(size = 30, colour = "black"),
#      legend.position = "none",
#      panel.grid = element_blank(),
#      panel.border = element_rect(size = 3.5, colour = "black"),
#      panel.background = element_blank(),
#      plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
#      plot.title = element_text(hjust = 0.5, size = 30)) +
#ggtitle(bquote(.(ranLocNamePlot) ~ "(" * italic("n") ~ "=" ~
#               .(prettyNum(summaryDFfeature$n[1],
#                           big.mark = ",", trim = T)) *
#               ")"))
#
### gap
#summaryDFfeature <- summaryDFfeature_log2ChIP[[3]]
#ggObj3_combined_log2ChIP <- ggplot(data = summaryDFfeature,
#                                   mapping = aes(x = winNo,
#                                                 y = mean,
#                                                 group = libName)
#                                  ) +
#geom_line(data = summaryDFfeature,
#          mapping = aes(colour = libName),
#          size = 1) +
#scale_colour_manual(values = log2ChIPColours) +
#geom_ribbon(data = summaryDFfeature,
#            mapping = aes(ymin = CI_lower,
#                          ymax = CI_upper,
#                          fill = libName),
#            alpha = 0.4) +
#scale_fill_manual(values = log2ChIPColours) +
#scale_y_continuous(limits = c(ymin_log2ChIP, ymax_log2ChIP),
#                   labels = function(x) sprintf("%6.3f", x)) +
#scale_x_discrete(breaks = c(1,
#                            (upstream/binSize)+1,
#                            (dim(summaryDFfeature_log2ChIP[[3]])[1]/length(ChIPNames))-(downstream/binSize),
#                            dim(summaryDFfeature_log2ChIP[[3]])[1]/length(ChIPNames)),
#                 labels = c(paste0("-", flankName),
#                            featureStartLab,
#                            featureEndLab,
#                            paste0("+", flankName))) +
#geom_vline(xintercept = c((upstream/binSize)+1,
#                          (dim(summaryDFfeature_log2ChIP[[3]])[1]/length(ChIPNames))-(downstream/binSize)),
#           linetype = "dashed",
#           size = 1) +
#labs(x = "",
#     y = bquote("Log"[2] * "(" * .(yLabPlot) * "/control)")) +
#theme_bw() +
#theme(
#      axis.ticks = element_line(size = 1.0, colour = "black"),
#      axis.ticks.length = unit(0.25, "cm"),
#      axis.text.x = element_text(size = 22, colour = "black"),
#      axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
#      axis.title = element_text(size = 30, colour = "black"),
#      legend.position = "none",
#      panel.grid = element_blank(),
#      panel.border = element_rect(size = 3.5, colour = "black"),
#      panel.background = element_blank(),
#      plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
#      plot.title = element_text(hjust = 0.5, size = 30)) +
#ggtitle(bquote(.(gapNamePlot) ~ "(" * italic("n") ~ "=" ~
#               .(prettyNum(summaryDFfeature$n[1],
#                           big.mark = ",", trim = T)) *
#               ")"))
#
### Athila
#summaryDFfeature <- summaryDFfeature_log2ChIP[[4]]
#ggObj4_combined_log2ChIP <- ggplot(data = summaryDFfeature,
#                                   mapping = aes(x = winNo,
#                                                 y = mean,
#                                                 group = libName)
#                                  ) +
#geom_line(data = summaryDFfeature,
#          mapping = aes(colour = libName),
#          size = 1) +
#scale_colour_manual(values = log2ChIPColours) +
#geom_ribbon(data = summaryDFfeature,
#            mapping = aes(ymin = CI_lower,
#                          ymax = CI_upper,
#                          fill = libName),
#            alpha = 0.4) +
#scale_fill_manual(values = log2ChIPColours) +
#scale_y_continuous(limits = c(ymin_log2ChIP, ymax_log2ChIP),
#                   labels = function(x) sprintf("%6.3f", x)) +
#scale_x_discrete(breaks = c(1,
#                            (upstream/Athila_binSize)+1,
#                            (dim(summaryDFfeature_log2ChIP[[4]])[1]/length(ChIPNames))-(downstream/Athila_binSize),
#                            dim(summaryDFfeature_log2ChIP[[4]])[1]/length(ChIPNames)),
#                 labels = c(paste0("-", flankName),
#                            featureStartLab,
#                            featureEndLab,
#                            paste0("+", flankName))) +
#geom_vline(xintercept = c((upstream/Athila_binSize)+1,
#                          (dim(summaryDFfeature_log2ChIP[[4]])[1]/length(ChIPNames))-(downstream/Athila_binSize)),
#           linetype = "dashed",
#           size = 1) +
#labs(x = "",
#     y = bquote("Log"[2] * "(" * .(yLabPlot) * "/control)")) +
#theme_bw() +
#theme(
#      axis.ticks = element_line(size = 1.0, colour = "black"),
#      axis.ticks.length = unit(0.25, "cm"),
#      axis.text.x = element_text(size = 22, colour = "black"),
#      axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
#      axis.title = element_text(size = 30, colour = "black"),
#      legend.position = "none",
#      panel.grid = element_blank(),
#      panel.border = element_rect(size = 3.5, colour = "black"),
#      panel.background = element_blank(),
#      plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
#      plot.title = element_text(hjust = 0.5, size = 30)) +
#ggtitle(bquote(.(AthilaNamePlot) ~ "(" * italic("n") ~ "=" ~
#               .(prettyNum(summaryDFfeature$n[1],
#                           big.mark = ",", trim = T)) *
#               ")"))
#
#ggObjGA_combined <- grid.arrange(grobs = list(
#                                              ggObj1_combined_log2ChIP,
#                                              ggObj2_combined_log2ChIP,
#                                              ggObj3_combined_log2ChIP,
#                                              ggObj4_combined_log2ChIP
#                                             ),
#                                 layout_matrix = cbind(
#                                                       1,
#                                                       2,
#                                                       3,
#                                                       4
#                                                      ))
#ggsave(paste0(plotDir,
#              "log2ChIPcontrol_",
#              paste0(ChIPNames, collapse = "_"),
#              "_avgProfiles_around",
#              "_CEN180_ranLoc_CENgap_CENAthila_in_T2T_Col_",
#              paste0(chrName, collapse = "_"), "_", align, ".pdf"),
#       plot = ggObjGA_combined,
#       height = 6.5, width = 28, limitsize = FALSE)


# ChIP
# Add column names
for(x in seq_along(ChIP_featureMats)) {
  colnames(ChIP_featureMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                       paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                       paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))
  colnames(ChIP_ranLocMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                      paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                      paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))
  colnames(ChIP_gapMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                   paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                   paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))
  colnames(ChIP_AthilaMats[[x]]) <- c(paste0("u", 1:(upstream/Athila_binSize)),
                                      paste0("t", ((upstream/Athila_binSize)+1):((upstream+Athila_bodyLength)/Athila_binSize)),
                                      paste0("d", (((upstream+Athila_bodyLength)/Athila_binSize)+1):(((upstream+Athila_bodyLength)/Athila_binSize)+(downstream/Athila_binSize))))
  colnames(ChIP_TEsfMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                    paste0("t", ((upstream/binSize)+1):((upstream+TEsf_bodyLength)/binSize)),
                                    paste0("d", (((upstream+TEsf_bodyLength)/binSize)+1):(((upstream+TEsf_bodyLength)/binSize)+(downstream/binSize))))
  colnames(ChIP_geneMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                    paste0("t", ((upstream/binSize)+1):((upstream+gene_bodyLength)/binSize)),
                                    paste0("d", (((upstream+gene_bodyLength)/binSize)+1):(((upstream+gene_bodyLength)/binSize)+(downstream/binSize))))
}

# Create list of lists in which each element in the enclosing list corresponds to a library
# and the two elements in the nested list correspond to coverage matrices for features and random loci
ChIP_mats <- mclapply(seq_along(ChIP_featureMats), function(x) {
  list(
       # features
       ChIP_featureMats[[x]],
       # ranLocs
       ChIP_ranLocMats[[x]],
       # gaps
       ChIP_gapMats[[x]],
       # Athilas
       ChIP_AthilaMats[[x]],
       # TEsfs
       ChIP_TEsfMats[[x]],
       # genes
       ChIP_geneMats[[x]]
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

# Convert list of lists summaryDFfeature_list_ChIP into
# a list of single data.frames containing all meta-profiles for plotting
featureTmp <- lapply(seq_along(summaryDFfeature_list_ChIP), function(x) {
  summaryDFfeature_list_ChIP[[x]][[1]]
})
ranLocTmp <- lapply(seq_along(summaryDFfeature_list_ChIP), function(x) {
  summaryDFfeature_list_ChIP[[x]][[2]]
})
gapTmp <- lapply(seq_along(summaryDFfeature_list_ChIP), function(x) {
  summaryDFfeature_list_ChIP[[x]][[3]]
})
AthilaTmp <- lapply(seq_along(summaryDFfeature_list_ChIP), function(x) {
  summaryDFfeature_list_ChIP[[x]][[4]]
})
TEsfTmp <- lapply(seq_along(summaryDFfeature_list_ChIP), function(x) {
  summaryDFfeature_list_ChIP[[x]][[5]]
})
geneTmp <- lapply(seq_along(summaryDFfeature_list_ChIP), function(x) {
  summaryDFfeature_list_ChIP[[x]][[6]]
})
names(featureTmp) <- ChIPNamesPlot
names(ranLocTmp) <- ChIPNamesPlot
names(gapTmp) <- ChIPNamesPlot
names(AthilaTmp) <- ChIPNamesPlot
names(TEsfTmp) <- ChIPNamesPlot
names(geneTmp) <- ChIPNamesPlot
summaryDFfeature_ChIP <- list(
  bind_rows(featureTmp, .id = "libName"),
  bind_rows(ranLocTmp, .id = "libName"),
  bind_rows(gapTmp, .id = "libName"),
  bind_rows(AthilaTmp, .id = "libName"),
  bind_rows(TEsfTmp, .id = "libName"),
  bind_rows(geneTmp, .id = "libName")
)
for(x in seq_along(summaryDFfeature_ChIP)) {
  summaryDFfeature_ChIP[[x]]$libName <- factor(summaryDFfeature_ChIP[[x]]$libName,
                                                   levels = ChIPNamesPlot)
}

# Define y-axis limits
ymin_ChIP <- min(c(summaryDFfeature_ChIP[[1]]$CI_lower,
                   summaryDFfeature_ChIP[[2]]$CI_lower,
                   summaryDFfeature_ChIP[[3]]$CI_lower,
                   summaryDFfeature_ChIP[[4]]$CI_lower,
                   summaryDFfeature_ChIP[[5]]$CI_lower,
                   summaryDFfeature_ChIP[[6]]$CI_lower),
                 na.rm = T)
ymax_ChIP <- max(c(summaryDFfeature_ChIP[[1]]$CI_upper,
                   summaryDFfeature_ChIP[[2]]$CI_upper,
                   summaryDFfeature_ChIP[[3]]$CI_upper,
                   summaryDFfeature_ChIP[[4]]$CI_upper,
                   summaryDFfeature_ChIP[[5]]$CI_upper,
                   summaryDFfeature_ChIP[[6]]$CI_upper),
                 na.rm = T)

# Define legend labels
legendLabs <- lapply(seq_along(ChIPNamesPlot), function(x) {
  grobTree(textGrob(bquote(.(ChIPNamesPlot[x])),
                    x = legendPos[1], y = legendPos[2]-((x-1)*0.06), just = "left",
                    gp = gpar(col = ChIPColours[x], fontsize = 18)))
})

# Plot average profiles with 95% CI ribbon
## feature
summaryDFfeature <- summaryDFfeature_ChIP[[1]]
ggObj1_combined_ChIP <- ggplot(data = summaryDFfeature,
                                   mapping = aes(x = winNo,
                                                 y = mean,
                                                 group = libName)
                                  ) +
geom_line(data = summaryDFfeature,
          mapping = aes(colour = libName),
          size = 1) +
scale_colour_manual(values = ChIPColours) +
geom_ribbon(data = summaryDFfeature,
            mapping = aes(ymin = CI_lower,
                          ymax = CI_upper,
                          fill = libName),
            alpha = 0.4) +
scale_fill_manual(values = ChIPColours) +
scale_y_continuous(limits = c(ymin_ChIP, ymax_ChIP),
                   labels = function(x) sprintf("%6.3f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/binSize)+1,
                            (dim(summaryDFfeature_ChIP[[1]])[1]/length(ChIPNames))-(downstream/binSize),
                            dim(summaryDFfeature_ChIP[[1]])[1]/length(ChIPNames)),
                 labels = c(paste0("-", flankName),
                            featureStartLab,
                            featureEndLab,
                            paste0("+", flankName))) +
geom_vline(xintercept = c((upstream/binSize)+1,
                          (dim(summaryDFfeature_ChIP[[1]])[1]/length(ChIPNames))-(downstream/binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = bquote(.(yLabPlot))) +
theme_bw() +
theme(
      axis.ticks = element_line(size = 1.0, colour = "black"),
      axis.ticks.length = unit(0.25, "cm"),
      axis.text.x = element_text(size = 22, colour = "black"),
      axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
      axis.title = element_text(size = 30, colour = "black"),
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

## ranLoc
summaryDFfeature <- summaryDFfeature_ChIP[[2]]
ggObj2_combined_ChIP <- ggplot(data = summaryDFfeature,
                                   mapping = aes(x = winNo,
                                                 y = mean,
                                                 group = libName)
                                  ) +
geom_line(data = summaryDFfeature,
          mapping = aes(colour = libName),
          size = 1) +
scale_colour_manual(values = ChIPColours) +
geom_ribbon(data = summaryDFfeature,
            mapping = aes(ymin = CI_lower,
                          ymax = CI_upper,
                          fill = libName),
            alpha = 0.4) +
scale_fill_manual(values = ChIPColours) +
scale_y_continuous(limits = c(ymin_ChIP, ymax_ChIP),
                   labels = function(x) sprintf("%6.3f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/binSize)+1,
                            (dim(summaryDFfeature_ChIP[[2]])[1]/length(ChIPNames))-(downstream/binSize),
                            dim(summaryDFfeature_ChIP[[2]])[1]/length(ChIPNames)),
                 labels = c(paste0("-", flankName),
                            featureStartLab,
                            featureEndLab,
                            paste0("+", flankName))) +
geom_vline(xintercept = c((upstream/binSize)+1,
                          (dim(summaryDFfeature_ChIP[[2]])[1]/length(ChIPNames))-(downstream/binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = bquote(.(yLabPlot))) +
annotation_custom(legendLabs[[1]]) +
annotation_custom(legendLabs[[2]]) +
annotation_custom(legendLabs[[3]]) +
annotation_custom(legendLabs[[4]]) +
annotation_custom(legendLabs[[5]]) +
annotation_custom(legendLabs[[6]]) +
theme_bw() +
theme(
      axis.ticks = element_line(size = 1.0, colour = "black"),
      axis.ticks.length = unit(0.25, "cm"),
      axis.text.x = element_text(size = 22, colour = "black"),
      axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
      axis.title = element_text(size = 30, colour = "black"),
      legend.position = "none",
      panel.grid = element_blank(),
      panel.border = element_rect(size = 3.5, colour = "black"),
      panel.background = element_blank(),
      plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
      plot.title = element_text(hjust = 0.5, size = 30)) +
ggtitle(bquote(.(ranLocNamePlot) ~ "(" * italic("n") ~ "=" ~
               .(prettyNum(summaryDFfeature$n[1],
                           big.mark = ",", trim = T)) *
               ")"))

## gap
summaryDFfeature <- summaryDFfeature_ChIP[[3]]
ggObj3_combined_ChIP <- ggplot(data = summaryDFfeature,
                                   mapping = aes(x = winNo,
                                                 y = mean,
                                                 group = libName)
                                  ) +
geom_line(data = summaryDFfeature,
          mapping = aes(colour = libName),
          size = 1) +
scale_colour_manual(values = ChIPColours) +
geom_ribbon(data = summaryDFfeature,
            mapping = aes(ymin = CI_lower,
                          ymax = CI_upper,
                          fill = libName),
            alpha = 0.4) +
scale_fill_manual(values = ChIPColours) +
scale_y_continuous(limits = c(ymin_ChIP, ymax_ChIP),
                   labels = function(x) sprintf("%6.3f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/binSize)+1,
                            (dim(summaryDFfeature_ChIP[[3]])[1]/length(ChIPNames))-(downstream/binSize),
                            dim(summaryDFfeature_ChIP[[3]])[1]/length(ChIPNames)),
                 labels = c(paste0("-", flankName),
                            featureStartLab,
                            featureEndLab,
                            paste0("+", flankName))) +
geom_vline(xintercept = c((upstream/binSize)+1,
                          (dim(summaryDFfeature_ChIP[[3]])[1]/length(ChIPNames))-(downstream/binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = bquote(.(yLabPlot))) +
theme_bw() +
theme(
      axis.ticks = element_line(size = 1.0, colour = "black"),
      axis.ticks.length = unit(0.25, "cm"),
      axis.text.x = element_text(size = 22, colour = "black"),
      axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
      axis.title = element_text(size = 30, colour = "black"),
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

## Athila
summaryDFfeature <- summaryDFfeature_ChIP[[4]]
ggObj4_combined_ChIP <- ggplot(data = summaryDFfeature,
                                   mapping = aes(x = winNo,
                                                 y = mean,
                                                 group = libName)
                                  ) +
geom_line(data = summaryDFfeature,
          mapping = aes(colour = libName),
          size = 1) +
scale_colour_manual(values = ChIPColours) +
geom_ribbon(data = summaryDFfeature,
            mapping = aes(ymin = CI_lower,
                          ymax = CI_upper,
                          fill = libName),
            alpha = 0.4) +
scale_fill_manual(values = ChIPColours) +
scale_y_continuous(limits = c(ymin_ChIP, ymax_ChIP),
                   labels = function(x) sprintf("%6.3f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/Athila_binSize)+1,
                            (dim(summaryDFfeature_ChIP[[4]])[1]/length(ChIPNames))-(downstream/Athila_binSize),
                            dim(summaryDFfeature_ChIP[[4]])[1]/length(ChIPNames)),
                 labels = c(paste0("-", flankName),
                            featureStartLab,
                            featureEndLab,
                            paste0("+", flankName))) +
geom_vline(xintercept = c((upstream/Athila_binSize)+1,
                          (dim(summaryDFfeature_ChIP[[4]])[1]/length(ChIPNames))-(downstream/Athila_binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = bquote(.(yLabPlot))) +
theme_bw() +
theme(
      axis.ticks = element_line(size = 1.0, colour = "black"),
      axis.ticks.length = unit(0.25, "cm"),
      axis.text.x = element_text(size = 22, colour = "black"),
      axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
      axis.title = element_text(size = 30, colour = "black"),
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

## TEsf
summaryDFfeature <- summaryDFfeature_ChIP[[5]]
ggObj5_combined_ChIP <- ggplot(data = summaryDFfeature,
                                   mapping = aes(x = winNo,
                                                 y = mean,
                                                 group = libName)
                                  ) +
geom_line(data = summaryDFfeature,
          mapping = aes(colour = libName),
          size = 1) +
scale_colour_manual(values = ChIPColours) +
geom_ribbon(data = summaryDFfeature,
            mapping = aes(ymin = CI_lower,
                          ymax = CI_upper,
                          fill = libName),
            alpha = 0.4) +
scale_fill_manual(values = ChIPColours) +
scale_y_continuous(limits = c(ymin_ChIP, ymax_ChIP),
                   labels = function(x) sprintf("%6.3f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/binSize)+1,
                            (dim(summaryDFfeature_ChIP[[5]])[1]/length(ChIPNames))-(downstream/binSize),
                            dim(summaryDFfeature_ChIP[[5]])[1]/length(ChIPNames)),
                 labels = c(paste0("-", flankName),
                            featureStartLab,
                            featureEndLab,
                            paste0("+", flankName))) +
geom_vline(xintercept = c((upstream/binSize)+1,
                          (dim(summaryDFfeature_ChIP[[5]])[1]/length(ChIPNames))-(downstream/binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = bquote(.(yLabPlot))) +
theme_bw() +
theme(
      axis.ticks = element_line(size = 1.0, colour = "black"),
      axis.ticks.length = unit(0.25, "cm"),
      axis.text.x = element_text(size = 22, colour = "black"),
      axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
      axis.title = element_text(size = 30, colour = "black"),
      legend.position = "none",
      panel.grid = element_blank(),
      panel.border = element_rect(size = 3.5, colour = "black"),
      panel.background = element_blank(),
      plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
      plot.title = element_text(hjust = 0.5, size = 30)) +
ggtitle(bquote(.(TEsfNamePlot) ~ "(" * italic("n") ~ "=" ~
               .(prettyNum(summaryDFfeature$n[1],
                           big.mark = ",", trim = T)) *
               ")"))

## gene
summaryDFfeature <- summaryDFfeature_ChIP[[6]]
ggObj6_combined_ChIP <- ggplot(data = summaryDFfeature,
                                   mapping = aes(x = winNo,
                                                 y = mean,
                                                 group = libName)
                                  ) +
geom_line(data = summaryDFfeature,
          mapping = aes(colour = libName),
          size = 1) +
scale_colour_manual(values = ChIPColours) +
geom_ribbon(data = summaryDFfeature,
            mapping = aes(ymin = CI_lower,
                          ymax = CI_upper,
                          fill = libName),
            alpha = 0.4) +
scale_fill_manual(values = ChIPColours) +
scale_y_continuous(limits = c(ymin_ChIP, ymax_ChIP),
                   labels = function(x) sprintf("%6.3f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/binSize)+1,
                            (dim(summaryDFfeature_ChIP[[6]])[1]/length(ChIPNames))-(downstream/binSize),
                            dim(summaryDFfeature_ChIP[[6]])[1]/length(ChIPNames)),
                 labels = c(paste0("-", flankName),
                            geneStartLab,
                            geneEndLab,
                            paste0("+", flankName))) +
geom_vline(xintercept = c((upstream/binSize)+1,
                          (dim(summaryDFfeature_ChIP[[6]])[1]/length(ChIPNames))-(downstream/binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = bquote(.(yLabPlot))) +
theme_bw() +
theme(
      axis.ticks = element_line(size = 1.0, colour = "black"),
      axis.ticks.length = unit(0.25, "cm"),
      axis.text.x = element_text(size = 22, colour = "black"),
      axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
      axis.title = element_text(size = 30, colour = "black"),
      legend.position = "none",
      panel.grid = element_blank(),
      panel.border = element_rect(size = 3.5, colour = "black"),
      panel.background = element_blank(),
      plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
      plot.title = element_text(hjust = 0.5, size = 30)) +
ggtitle(bquote(.(geneNamePlot) ~ "(" * italic("n") ~ "=" ~
               .(prettyNum(summaryDFfeature$n[1],
                           big.mark = ",", trim = T)) *
               ")"))

ggObjGA_combined <- grid.arrange(grobs = list(
                                              ggObj1_combined_ChIP,
                                              ggObj2_combined_ChIP,
                                              ggObj3_combined_ChIP,
                                              ggObj4_combined_ChIP,
                                              ggObj5_combined_ChIP,
                                              ggObj6_combined_ChIP
                                             ),
                                 layout_matrix = cbind(
                                                       1,
                                                       2,
                                                       3,
                                                       4,
                                                       5,
                                                       6
                                                      ))
ggsave(paste0(plotDir,
              "ChIP_",
              paste0(ChIPNames, collapse = "_"),
              "_avgProfiles_around",
              "_CEN180_ranLoc_CENgap_CENAthila_Gypsy_LTR_genes_in_T2T_Col_",
              paste0(chrName, collapse = "_"), "_", align, ".pdf"),
       plot = ggObjGA_combined,
       height = 6.5, width = 7*6, limitsize = FALSE)
