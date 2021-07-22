#!/applications/R/R-4.0.0/bin/Rscript

# author: Andy Tock
# contact: ajt200@cam.ac.uk
# date: 18.12.2020

# Calculate and plot metaprofiles of ChIP-seq, MNase-seq, etc.
# (CEN180 windowed means and 95% confidence intervals, CIs)
# for all CEN180 sequences and randomly positioned loci

# Usage:
# /applications/R/R-4.0.0/bin/Rscript gene_2metaprofiles_DAN.R 'Chr1,Chr2,Chr3,Chr4,Chr5' both 2200 2000 2kb 10 10bp '0.02,0.96' 'WT_MTOPVIB_HA_Rep1_ChIP,WT_MTOPVIB_HA_Rep2_ChIP' '20190819_dh580_Athaliana_ChIPseq_MTOPVIB/fastq_pooled/snakemake_ChIPseq_t2t-col.20210610,20190819_dh580_Athaliana_ChIPseq_MTOPVIB/fastq_pooled/snakemake_ChIPseq_t2t-col.20210610' 'WT MTOPVIB Rep1,WT MTOPVIB Rep2' 'dodgerblue1,red' 'ChIP'

#chrName <- unlist(strsplit("Chr1,Chr2,Chr3,Chr4,Chr5", split = ","))
#align <- "both"
#gene_bodyLength <- 2200
#upstream <- 2000
#downstream <- 2000
#flankName <- "2kb"
#binSize <- 10
#binName <- "10bp"
#legendPos <- as.numeric(unlist(strsplit("0.02,0.96", split = ",")))
#ChIPNames <- unlist(strsplit("WT_MTOPVIB_HA_Rep1_ChIP,WT_MTOPVIB_HA_Rep2_ChIP", split = ","))
#ChIPNamesDir <- unlist(strsplit("20190819_dh580_Athaliana_ChIPseq_MTOPVIB/fastq_pooled/snakemake_ChIPseq_t2t-col.20210610,20190819_dh580_Athaliana_ChIPseq_MTOPVIB/fastq_pooled/snakemake_ChIPseq_t2t-col.20210610", split = ","))
#log2ChIPNamesPlot <- unlist(strsplit("WT MTOPVIB Rep1,WT MTOPVIB Rep2", split = ","))
#log2ChIPColours <- unlist(strsplit("dodgerblue1,red", split = ","))
#yLabPlot <- "ChIP"

args <- commandArgs(trailingOnly = T)
chrName <- unlist(strsplit(args[1],
                           split = ","))
align <- args[2]
gene_bodyLength <- as.numeric(args[3])
upstream <- as.numeric(args[4])
downstream <- as.numeric(args[4])
flankName <- args[5]
binSize <- as.numeric(args[6])
binName <- args[7]
legendPos <- as.numeric(unlist(strsplit(args[8],
                                        split = ",")))
ChIPNames <- unlist(strsplit(args[9],
                             split = ","))
ChIPNamesDir <- unlist(strsplit(args[10],
                                split = ","))
log2ChIPNamesPlot <- unlist(strsplit(args[11],
                                     split = ","))
log2ChIPColours <- unlist(strsplit(args[12],
                                   split = ","))
yLabPlot <- args[13]

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
  geneNamePlot <- "All genes"
  ranLocNamePlot <- "All ranLoc"
} else {
  geneNamePlot <- paste0(paste0(chrName, collapse = ","), " genes")
  ranLocNamePlot <- paste0(paste0(chrName, collapse = ","), " ranLoc")
}

# Define feature start and end labels for plotting
featureStartLab <- "Start"
featureEndLab <- "End"
geneStartLab <- "TSS"
geneEndLab <- "TTS"

# Load feature matrices for each chromatin dataset, calculate log2(ChIP/control),
log2ChIPNames <- ChIPNames
ChIPNamesPlot <- log2ChIPNamesPlot
ChIPColours <- log2ChIPColours
# Load feature matrices for each chromatin dataset, calculate log2(ChIP/control),
ChIPDirs <- sapply(seq_along(ChIPNames), function(x) {
  paste0("/home/ajt200/analysis/",
         ChIPNamesDir[x],
         "/mapped/")
})

controlNames <- c(
  "WT_REC8_Myc_Rep1_input"
)
controlNamesDir <- c(
  "REC8_pooled/snakemake_ChIPseq_t2t-col.20210610"
)
controlNamesPlot <- c(
  "PE input"
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
# gene
ChIP_geneMats <- mclapply(seq_along(ChIPNames), function(x) {
  lapply(seq_along(chrName), function(y) {
    as.matrix(read.table(paste0(ChIPDirs[x], "geneProfiles/matrices/",
                                ChIPNames[x],
                                "_MappedOn_t2t-col.20210610_lowXM_", align, "_sort_norm_genes_in_",
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
                                "_MappedOn_t2t-col.20210610_lowXM_", align, "_sort_norm_genes_in_",
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


## ChIP
# ranLoc
ChIP_ranLocMats <- mclapply(seq_along(ChIPNames), function(x) {
  lapply(seq_along(chrName), function(y) {
    as.matrix(read.table(paste0(ChIPDirs[x], "geneProfiles/matrices/",
                                ChIPNames[x],
                                "_MappedOn_t2t-col.20210610_lowXM_", align, "_sort_norm_genes_in_",
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
    as.matrix(read.table(paste0(controlDirs[x], "geneProfiles/matrices/",
                                controlNames[x],
                                "_MappedOn_t2t-col.20210610_lowXM_", align, "_sort_norm_genes_in_",
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
# gene
log2ChIP_geneMats <- mclapply(seq_along(ChIP_geneMats), function(x) {
  if ( grepl("MNase", ChIPNames[x]) ) {
    print(paste0(ChIPNames[x], " library; using ", controlNames[1], " for log2((MNase+1)/(gDNA+1)) calculation"))
    log2((ChIP_geneMats[[x]]+1)/(control_geneMats[[1]]+1))
  } else if ( grepl("SPO11oligos", ChIPNames[x]) ) {
    print(paste0(ChIPNames[x], " library; using ", controlNames[2], " for log2((SPO11-1-oligos+1)/(gDNA+1)) calculation"))
    log2((ChIP_geneMats[[x]]+1)/(control_geneMats[[2]]+1))
  } else {
    print(paste0(ChIPNames[x], " library; using ", controlNames[1], " for log2((ChIP+1)/(input+1)) calculation"))
    log2((ChIP_geneMats[[x]]+1)/(control_geneMats[[1]]+1))
  }
}, mc.cores = length(ChIP_geneMats))


# ranLoc
log2ChIP_ranLocMats <- mclapply(seq_along(ChIP_ranLocMats), function(x) {
  if ( grepl("MNase", ChIPNames[x]) ) {
    print(paste0(ChIPNames[x], " library; using ", controlNames[1], " for log2((MNase+1)/(gDNA+1)) calculation"))
    log2((ChIP_ranLocMats[[x]]+1)/(control_ranLocMats[[1]]+1))
  } else if ( grepl("SPO11oligos", ChIPNames[x]) ) {
    print(paste0(ChIPNames[x], " library; using ", controlNames[2], " for log2((SPO11-1-oligos+1)/(gDNA+1)) calculation"))
    log2((ChIP_ranLocMats[[x]]+1)/(control_ranLocMats[[2]]+1))
  } else {
    print(paste0(ChIPNames[x], " library; using ", controlNames[1], " for log2((ChIP+1)/(input+1)) calculation"))
    log2((ChIP_ranLocMats[[x]]+1)/(control_ranLocMats[[1]]+1))
  }
}, mc.cores = length(ChIP_ranLocMats))
# log2ChIP
# Add column names
for(x in seq_along(log2ChIP_geneMats)) {
  colnames(log2ChIP_geneMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                        paste0("t", ((upstream/binSize)+1):((upstream+gene_bodyLength)/binSize)),
                                        paste0("d", (((upstream+gene_bodyLength)/binSize)+1):(((upstream+gene_bodyLength)/binSize)+(downstream/binSize))))
  colnames(log2ChIP_ranLocMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                        paste0("t", ((upstream/binSize)+1):((upstream+gene_bodyLength)/binSize)),
                                        paste0("d", (((upstream+gene_bodyLength)/binSize)+1):(((upstream+gene_bodyLength)/binSize)+(downstream/binSize))))
}

# Create list of lists in which each element in the enclosing list corresponds to a library
# and the two elements in the nested list correspond to coverage matrices for features and random loci
log2ChIP_mats <- mclapply(seq_along(log2ChIP_geneMats), function(x) {
  list(
    
    # genes
    log2ChIP_geneMats[[x]],
    # ranLoc
    log2ChIP_ranLocMats[[x]]
  )
}, mc.cores = length(log2ChIP_geneMats))

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

# Convert list of lists summaryDFfeature_list_log2ChIP into
# a list of single data.frames containing all meta-profiles for plotting

geneTmp <- lapply(seq_along(summaryDFfeature_list_log2ChIP), function(x) {
  summaryDFfeature_list_log2ChIP[[x]][[1]]
})

ranLocTmp <- lapply(seq_along(summaryDFfeature_list_log2ChIP), function(x) {
  summaryDFfeature_list_log2ChIP[[x]][[2]]
})
names(geneTmp) <- log2ChIPNamesPlot
names(ranLocTmp) <- log2ChIPNamesPlot
summaryDFfeature_log2ChIP <- list(
  bind_rows(geneTmp, .id = "libName"),
  bind_rows(ranLocTmp, .id = "libName")
)
for(x in seq_along(summaryDFfeature_log2ChIP)) {
  summaryDFfeature_log2ChIP[[x]]$libName <- factor(summaryDFfeature_log2ChIP[[x]]$libName,
                                                   levels = log2ChIPNamesPlot)
}

# Define y-axis limits
ymin_log2ChIP <- min(c(summaryDFfeature_log2ChIP[[1]]$CI_lower,
                       summaryDFfeature_log2ChIP[[2]]$CI_lower),
                     na.rm = T)
ymax_log2ChIP <- max(c(summaryDFfeature_log2ChIP[[1]]$CI_upper,
                       summaryDFfeature_log2ChIP[[2]]$CI_upper),
                     na.rm = T)

# Define legend labels
legendLabs <- lapply(seq_along(log2ChIPNamesPlot), function(x) {
  grobTree(textGrob(bquote(.(log2ChIPNamesPlot[x])),
                    x = legendPos[1], y = legendPos[2]-((x-1)*0.06), just = "left",
                    gp = gpar(col = log2ChIPColours[x], fontsize = 18)))
})

# Plot average profiles with 95% CI ribbon

## gene
summaryDFfeature <- summaryDFfeature_log2ChIP[[1]]
ggObj1_combined_log2ChIP <- ggplot(data = summaryDFfeature,
                                   mapping = aes(x = winNo,
                                                 y = mean,
                                                 group = libName)
) +
  geom_line(data = summaryDFfeature,
            mapping = aes(colour = libName),
            size = 1) +
  scale_colour_manual(values = log2ChIPColours) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = libName),
              alpha = 0.4) +
  scale_fill_manual(values = log2ChIPColours) +
  scale_y_continuous(limits = c(ymin_log2ChIP, ymax_log2ChIP),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_log2ChIP[[1]])[1]/length(log2ChIPNames))-(downstream/binSize),
                              dim(summaryDFfeature_log2ChIP[[1]])[1]/length(log2ChIPNames)),
                   labels = c(paste0("-", flankName),
                              geneStartLab,
                              geneEndLab,
                              paste0("+", flankName))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_log2ChIP[[1]])[1]/length(log2ChIPNames))-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = bquote("Log"[2] * "(" * .(yLabPlot) * "/control)")) +
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

## ranLoc
summaryDFfeature <- summaryDFfeature_log2ChIP[[2]]
ggObj2_combined_log2ChIP <- ggplot(data = summaryDFfeature,
                                   mapping = aes(x = winNo,
                                                 y = mean,
                                                 group = libName)
) +
  geom_line(data = summaryDFfeature,
            mapping = aes(colour = libName),
            size = 1) +
  scale_colour_manual(values = log2ChIPColours) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = libName),
              alpha = 0.4) +
  scale_fill_manual(values = log2ChIPColours) +
  scale_y_continuous(limits = c(ymin_log2ChIP, ymax_log2ChIP),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_log2ChIP[[2]])[1]/length(log2ChIPNames))-(downstream/binSize),
                              dim(summaryDFfeature_log2ChIP[[2]])[1]/length(log2ChIPNames)),
                   labels = c(paste0("-", flankName),
                              "Start",
                              "End",
                              paste0("+", flankName))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_log2ChIP[[2]])[1]/length(log2ChIPNames))-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = bquote("Log"[2] * "(" * .(yLabPlot) * "/control)")) +
  annotation_custom(legendLabs[[1]]) + 
  annotation_custom(legendLabs[[2]]) +
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

ggObjGA_combined <- grid.arrange(grobs = list(
  ggObj1_combined_log2ChIP,
  ggObj2_combined_log2ChIP
),
layout_matrix = cbind(
  1,
  2
))
ggsave(paste0(plotDir,
              "log2ChIPcontrol_",
              paste0(log2ChIPNames, collapse = "_"),
              "_avgProfiles_around",
              "_genes_ranLoc_in_t2t-col.20210610_",
              paste0(chrName, collapse = "_"), "_", align, ".pdf"),
       plot = ggObjGA_combined,
       height = 6.5, width = 7*2, limitsize = FALSE)
