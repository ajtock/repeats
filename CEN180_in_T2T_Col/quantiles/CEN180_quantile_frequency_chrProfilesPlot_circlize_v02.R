#!/applications/R/R-4.0.0/bin/Rscript

# author: Andy Tock
# contact: ajt200@cam.ac.uk
# date: 16.02.2021

# Plot windowed CEN180 frequencies for CEN180 sequences within orderingFactor quantiles

# Usage:
# /applications/R/R-4.0.0/bin/Rscript CEN180_quantile_frequency_chrProfilesPlot_circlize_v02.R 'Chr1,Chr2,Chr3,Chr4,Chr5' 4 10000 perchrom 260221

#chrName <- unlist(strsplit("Chr1,Chr2,Chr3,Chr4,Chr5",
#                           split = ","))
#quantiles <- 4
#genomeBinSize <- 10000
#quantileDef <- "genomewide"
#quantileDef <- "perchrom"
#date <- "260221"

args <- commandArgs(trailingOnly = T)
chrName <- unlist(strsplit(args[1],
                           split = ","))
quantiles <- as.integer(args[2])
genomeBinSize <- as.integer(args[3])
quantileDef <- args[4]
date <- as.character(args[5])

if(floor(log10(genomeBinSize)) + 1 < 4) {
  genomeBinName <- paste0(genomeBinSize, "bp")
} else if(floor(log10(genomeBinSize)) + 1 >= 4 &
          floor(log10(genomeBinSize)) + 1 <= 6) {
  genomeBinName <- paste0(genomeBinSize/1e3, "kb")
} else if(floor(log10(genomeBinSize)) + 1 >= 7) {
  genomeBinName <- paste0(genomeBinSize/1e6, "Mb")
}

options(stringsAsFactors = F)
options(scipen=999)
library(parallel)
library(GenomicRanges)
library(circlize)
library(ComplexHeatmap)
library(gridBase)
library(viridis)

plotDir <- paste0("plots/")
system(paste0("[ -d ", plotDir, " ] || mkdir -p ", plotDir))

# Define quantile colours
quantileColours <- c("red", "orange", "dodgerblue", "navy")
makeTransparent <- function(thisColour, alpha = 210)
{
  newColour <- col2rgb(thisColour)
  apply(newColour, 2, function(x) {
    rgb(red = x[1], green = x[2], blue = x[3],
        alpha = alpha, maxColorValue = 255)
  })
}
quantileColours <- makeTransparent(quantileColours)

# Genomic definitions
fai <- read.table("/home/ajt200/analysis/nanopore/T2T_Col/T2T_Col.fa.fai", header = F)
chrs <- fai$V1[which(fai$V1 %in% chrName)]
chrLens <- fai$V2[which(fai$V1 %in% chrName)]
regionGR <- GRanges(seqnames = chrs,
                    ranges = IRanges(start = 1,
                                     end = chrLens),
                    strand = "*")
CENstart <- c(14840750, 3724530, 13597090, 4203495, 11783990)[which(fai$V1 %in% chrName)]
CENend <- c(17558182, 5946091, 15733029, 6977107, 14551874)[which(fai$V1 %in% chrName)]
CENGR <- GRanges(seqnames = chrs,
                 ranges = IRanges(start = CENstart,
                                  end = CENend),
                 strand = "*")

if(quantileDef == "genomewide") {
  # Load table of windowed feature frequencies for each quantile
  CENH3_in_bodies_CEN180 <- read.table(paste0("quantiles_by_CENH3_in_bodies/", paste0(chrName, collapse = "_"), "/",
                                              "CEN180_frequency_per_", genomeBinName,
                                              "_", quantiles, "quantiles_",
                                              "_by_CENH3_in_bodies",
                                              "_of_CEN180_in_T2T_Col_",
                                              paste0(chrName, collapse = "_"), "_smoothed.tsv"),
                                       header = T, sep = "\t")
  HORlengthsSum_CEN180 <- read.table(paste0("quantiles_by_HORlengthsSum/", paste0(chrName, collapse = "_"), "/",
                                            "CEN180_frequency_per_", genomeBinName,
                                            "_", quantiles, "quantiles_",
                                            "_by_HORlengthsSum",
                                            "_of_CEN180_in_T2T_Col_",
                                            paste0(chrName, collapse = "_"), "_smoothed.tsv"),
                                     header = T, sep = "\t")
  wSNV_CEN180 <- read.table(paste0("quantiles_by_wSNV/", paste0(chrName, collapse = "_"), "/",
                                   "CEN180_frequency_per_", genomeBinName,
                                   "_", quantiles, "quantiles_",
                                   "_by_wSNV",
                                   "_of_CEN180_in_T2T_Col_",
                                   paste0(chrName, collapse = "_"), "_smoothed.tsv"),
                            header = T, sep = "\t")
} else if(quantileDef == "perchrom") {
  # Load table of windowed feature frequencies for each quantile
  CENH3_in_bodies_CEN180List <- lapply(seq_along(chrName), function(x) {  
    read.table(paste0("quantiles_by_CENH3_in_bodies/", chrName[x], "/",
                      "CEN180_frequency_per_", genomeBinName,
                      "_", quantiles, "quantiles_",
                      "_by_CENH3_in_bodies",
                      "_of_CEN180_in_T2T_Col_",
                      chrName[x], "_smoothed.tsv"),
               header = T, sep = "\t")
  })
  if(length(chrName) > 1) {
    CENH3_in_bodies_CEN180 <- do.call(rbind, CENH3_in_bodies_CEN180List)
  } else {
    CENH3_in_bodies_CEN180 <- CENH3_in_bodies_CEN180List[[1]]
  }
  HORlengthsSum_CEN180List <- lapply(seq_along(chrName), function(x) {  
    read.table(paste0("quantiles_by_HORlengthsSum/", chrName[x], "/",
                      "CEN180_frequency_per_", genomeBinName,
                      "_", quantiles, "quantiles_",
                      "_by_HORlengthsSum",
                      "_of_CEN180_in_T2T_Col_",
                      chrName[x], "_smoothed.tsv"),
               header = T, sep = "\t")
  })
  if(length(chrName) > 1) {
    HORlengthsSum_CEN180 <- do.call(rbind, HORlengthsSum_CEN180List)
  } else {
    HORlengthsSum_CEN180 <- HORlengthsSum_CEN180List[[1]]
  }
  wSNV_CEN180List <- lapply(seq_along(chrName), function(x) {  
    read.table(paste0("quantiles_by_wSNV/", chrName[x], "/",
                      "CEN180_frequency_per_", genomeBinName,
                      "_", quantiles, "quantiles_",
                      "_by_wSNV",
                      "_of_CEN180_in_T2T_Col_",
                      chrName[x], "_smoothed.tsv"),
               header = T, sep = "\t")
  })
  if(length(chrName) > 1) {
    wSNV_CEN180 <- do.call(rbind, wSNV_CEN180List)
  } else {
    wSNV_CEN180 <- wSNV_CEN180List[[1]]
  }
}

# Load table of windowed log2(ChIP/input) coverage for CENH3 ChIP-seq
CENH3 <- read.table(paste0("/home/ajt200/analysis/CENH3_seedlings_Maheshwari_Comai_2017_GenomeRes/",
                           "snakemake_ChIPseq_T2T_Col/mapped/both/tsv/log2ChIPcontrol/",
                           "WT_CENH3_Rep1_ChIP_SRR4430537_WT_REC8_Myc_Rep1_input_MappedOn_T2T_Col_lowXM_both_sort_norm_binSize",
                           genomeBinName, "_smoothed.tsv"),
                    header = T, sep = "\t")
# Load table of windowed superfam TE frequencies
Gypsy <- read.table(paste0("/home/ajt200/analysis/nanopore/T2T_Col/annotation/TEs/T2T_Col_TEs_Gypsy_LTR_frequency_per_",
                           genomeBinName, "_smoothed.tsv"),
                    header = T, sep = "\t")
# Load table of windowed CEN180 frequencies
CEN180freq <- read.table(paste0("/home/ajt200/analysis/repeats/CEN180_in_T2T_Col/T2T_Col_CEN180_frequency_per_",
                                genomeBinName, "_smoothed.tsv"),
                         header = T, sep = "\t")

# Convert to BED-like format for use with circlize
CENH3_in_bodies_CEN180_bed <- data.frame(chr = CENH3_in_bodies_CEN180$chr,
                                         start = CENH3_in_bodies_CEN180$window-1,
                                         end = CENH3_in_bodies_CEN180$window-1+genomeBinSize,
                                         value1 = CENH3_in_bodies_CEN180[,3],
                                         value2 = CENH3_in_bodies_CEN180[,4],
                                         value3 = CENH3_in_bodies_CEN180[,5],
                                         value4 = CENH3_in_bodies_CEN180[,6])
HORlengthsSum_CEN180_bed <- data.frame(chr = HORlengthsSum_CEN180$chr,
                                       start = HORlengthsSum_CEN180$window-1,
                                       end = HORlengthsSum_CEN180$window-1+genomeBinSize,
                                       value1 = HORlengthsSum_CEN180[,3],
                                       value2 = HORlengthsSum_CEN180[,4],
                                       value3 = HORlengthsSum_CEN180[,5],
                                       value4 = HORlengthsSum_CEN180[,6])
wSNV_CEN180_bed <- data.frame(chr = wSNV_CEN180$chr,
                              start = wSNV_CEN180$window-1,
                              end = wSNV_CEN180$window-1+genomeBinSize,
                              value1 = wSNV_CEN180[,3],
                              value2 = wSNV_CEN180[,4],
                              value3 = wSNV_CEN180[,5],
                              value4 = wSNV_CEN180[,6])
CENH3_bed <- data.frame(chr = CENH3$chr,
                        start = CENH3$window-1,
                        end = CENH3$window-1+genomeBinSize,
                        value1 = CENH3[,6])
Gypsy_bed <- data.frame(chr = Gypsy$chr,
                        start = Gypsy$window-1,
                        end = Gypsy$window-1+genomeBinSize,
                        value1 = Gypsy[,4])
CEN180freq_bed <- data.frame(chr = CEN180freq$chr,
                             start = CEN180freq$window-1,
                             end = CEN180freq$window-1+genomeBinSize,
                             value1 = CEN180freq[,4])

# Redefine end coordinate of last window to match chrLens for each chromosome
for(x in seq_along(chrs)) {
  CENH3_in_bodies_CEN180_bed[CENH3_in_bodies_CEN180_bed$chr == chrs[x],][dim(CENH3_in_bodies_CEN180_bed[CENH3_in_bodies_CEN180_bed$chr == chrs[x],])[1],]$end <- chrLens[x]
  HORlengthsSum_CEN180_bed[HORlengthsSum_CEN180_bed$chr == chrs[x],][dim(HORlengthsSum_CEN180_bed[HORlengthsSum_CEN180_bed$chr == chrs[x],])[1],]$end <- chrLens[x]
  wSNV_CEN180_bed[wSNV_CEN180_bed$chr == chrs[x],][dim(wSNV_CEN180_bed[wSNV_CEN180_bed$chr == chrs[x],])[1],]$end <- chrLens[x]
  CENH3_bed[CENH3_bed$chr == chrs[x],][dim(CENH3_bed[CENH3_bed$chr == chrs[x],])[1],]$end <- chrLens[x]
  Gypsy_bed[Gypsy_bed$chr == chrs[x],][dim(Gypsy_bed[Gypsy_bed$chr == chrs[x],])[1],]$end <- chrLens[x]
  CEN180freq_bed[CEN180freq_bed$chr == chrs[x],][dim(CEN180freq_bed[CEN180freq_bed$chr == chrs[x],])[1],]$end <- chrLens[x]
}

# Define feature density heatmap colour scale interpolation
#rich8to6equal <- c("#000041", "#0000CB", "#0081FF", "#FDEE02", "#FFAB00", "#FF3300")
Gypsy_col_fun <- colorRamp2(quantile(Gypsy_bed$value1,
                                     c(0.5, 0.6, 0.7, 0.8, 0.9, 1.0),
                                     na.rm = T),
                            viridis_pal()(6))
CEN180freq_col_fun <- colorRamp2(quantile(CEN180freq_bed$value1,
                                          c(0.90, 0.92, 0.94, 0.96, 0.98, 1.0),
                                          na.rm = T),
                                 viridis_pal()(6))
# Define corresponding heatmap legends
lgd_Gypsy <- Legend(col_fun = Gypsy_col_fun, title = "Gypsy", title_gp = gpar(fontface = "bold.italic"), title_position = "leftcenter-rot")
lgd_CEN180freq <- Legend(col_fun = CEN180freq_col_fun, title = "CEN180", title_gp = gpar(fontface = "bold.italic"), title_position = "leftcenter-rot")
lgd_list1 <- packLegend(lgd_Gypsy, lgd_CEN180freq)

lgd_CENH3_in_bodies_CEN180 <- Legend(at = c("Q1", "Q2", "Q3", "Q4"), type = "lines", legend_gp = gpar(col = quantileColours, lwd = 1.5),
                                     background = NULL, title = "CENH3", title_position = "leftcenter-rot")
lgd_HORlengthsSum_CEN180 <- Legend(at = c("Q1", "Q2", "Q3", "Q4"), type = "lines", legend_gp = gpar(col = quantileColours, lwd = 1.5),
                                   background = NULL, title = "Repeats", title_position = "leftcenter-rot")
lgd_wSNV_CEN180 <- Legend(at = c("Q1", "Q2", "Q3", "Q4"), type = "lines", legend_gp = gpar(col = quantileColours, lwd = 1.5),
                          background = NULL, title = "SNVs", title_position = "leftcenter-rot")
lgd_list2 <- packLegend(lgd_CENH3_in_bodies_CEN180, lgd_HORlengthsSum_CEN180, lgd_wSNV_CEN180)


CEN180_bed <- read.table(paste0("/home/ajt200/analysis/repeats/CEN180_in_T2T_Col/CEN180_in_T2T_Col_",
                                paste0(chrName, collapse = "_"), ".bed"),
                         header = F, colClasses = c(rep(NA, 3), rep("NULL", 5)))
colnames(CEN180_bed) <- c("chr", "start", "end")
CENAthila_bed <- read.table(paste0("/home/ajt200/analysis/repeats/CEN180_in_T2T_Col/CENAthila_in_T2T_Col_",
                                   paste0(chrName, collapse = "_"), ".bed"),
                            header = F, colClasses = c(rep(NA, 3), rep("NULL", 2)))
colnames(CENAthila_bed) <- c("chr", "start", "end")


## circlize

genomeDF <- data.frame(chr = chrs,
                       start = c(12e6, 2e6, 11e6, 2e6, 10e6),
                       end = c(20e6, 8e6, 18e6, 10e6, 17e6))

# Initialize circular layout
circlize_plot <- function() {
  gapDegree <- 6
  circos.par(track.height = 0.15,
             canvas.xlim = c(-1.1, 1.1),
             canvas.ylim = c(-1.1, 1.1),
             gap.degree = c(rep(1, length(chrs)-1), gapDegree),
             start.degree = 90-(gapDegree/2))
  circos.genomicInitialize(data = genomeDF,
                           plotType = NULL,
                           tickLabelsStartFromZero = FALSE)
  circos.track(ylim = c(0, 1),
               bg.col = "grey70",
               bg.border = NA,
               track.height = 0.05,
               panel.fun = function(x, y) {
                 circos.genomicAxis(h = "top",
                                    labels.facing = "clockwise",
                                    tickLabelsStartFromZero = FALSE)
               })
  sapply(seq_along(chrs), function(x) {
    circos.text(x = (CENstart[x]+CENend[x])/2,
                y = 0.5,
                sector.index = chrs[x],
                track.index = 1,
                labels = paste0("CEN", x),
                niceFacing = TRUE,
                cex = 0.8,
                col = "white",
                font = 4)
  })
  
  # Gypsy heatmap
  circos.genomicHeatmap(bed = do.call(rbind, lapply(seq_along(chrs), function(x) {
    Gypsy_bed[Gypsy_bed$chr == chrs[x] &
              Gypsy_bed$start >= genomeDF$start[x] &
              Gypsy_bed$end <= genomeDF$end[x],] } )),
    col = Gypsy_col_fun,
    border = NA,
    side = "inside",
    heatmap_height = 0.05,
    connection_height = NULL)
  # CEN180 heatmap
  circos.genomicHeatmap(bed = do.call(rbind, lapply(seq_along(chrs), function(x) {
    CEN180freq_bed[CEN180freq_bed$chr == chrs[x] &
                   CEN180freq_bed$start >= genomeDF$start[x] &
                   CEN180freq_bed$end <= genomeDF$end[x],] } )),
    col = CEN180freq_col_fun,
    border = NA,
    side = "inside",
    heatmap_height = 0.05,
    connection_height = NULL)
  # CENAthila rainfall plot
  circos.genomicRainfall(data = do.call(rbind, lapply(seq_along(chrs), function(x) {
    CENAthila_bed[CENAthila_bed$chr == chrs[x] &
                  CENAthila_bed$start >= genomeDF$start[x] &
                  CENAthila_bed$end <= genomeDF$end[x],] } )),
                         bg.border = NA,
                         track.height = 0.05,
                         pch = 16,
                         cex = 0.4,
                         col = c("#0000FF80"))
  
  set_track_gap(gap = 0.005)

  # Plot windowed CEN180 frequency for each quantile
  circos.genomicTrack(data = do.call(rbind, lapply(seq_along(chrs), function(x) {
    CENH3_in_bodies_CEN180_bed[CENH3_in_bodies_CEN180_bed$chr == chrs[x] &
                               CENH3_in_bodies_CEN180_bed$start >= genomeDF$start[x] &
                               CENH3_in_bodies_CEN180_bed$end <= genomeDF$end[x],] } )),
                      panel.fun = function(region, value, ...) {
                        circos.genomicLines(region,
                                            value,
                                            col = quantileColours,
                                            lwd = 1.5,
                                            lty = 1,
                                            area = FALSE,
                                            ...)
                      }, bg.border = NA)
  circos.yaxis(side = "left", sector.index = get.all.sector.index()[1], labels.cex = 0.3, tick.length = convert_x(0.5, "mm"), lwd = 0.5)

  set_track_gap(gap = 0.005)

  circos.genomicTrack(data = do.call(rbind, lapply(seq_along(chrs), function(x) {
    HORlengthsSum_CEN180_bed[HORlengthsSum_CEN180_bed$chr == chrs[x] &
                             HORlengthsSum_CEN180_bed$start >= genomeDF$start[x] &
                             HORlengthsSum_CEN180_bed$end <= genomeDF$end[x],] } )),
                      panel.fun = function(region, value, ...) {
                        circos.genomicLines(region,
                                            value,
                                            col = quantileColours,
                                            lwd = 1.5,
                                            lty = 1,
                                            area = FALSE,
                                            ...)
                      }, bg.border = NA)
  circos.yaxis(side = "left", sector.index = get.all.sector.index()[1], labels.cex = 0.3, tick.length = convert_x(0.5, "mm"), lwd = 0.5)

  set_track_gap(gap = 0.005)

  circos.genomicTrack(data = do.call(rbind, lapply(seq_along(chrs), function(x) {
    wSNV_CEN180_bed[wSNV_CEN180_bed$chr == chrs[x] &
                    wSNV_CEN180_bed$start >= genomeDF$start[x] &
                    wSNV_CEN180_bed$end <= genomeDF$end[x],] } )),
                      panel.fun = function(region, value, ...) {
                        circos.genomicLines(region,
                                            value,
                                            col = quantileColours,
                                            lwd = 1.5,
                                            lty = 1,
                                            area = FALSE,
                                            ...)
                      }, bg.border = NA)
  circos.yaxis(side = "left", sector.index = get.all.sector.index()[1], labels.cex = 0.3, tick.length = convert_x(0.5, "mm"), lwd = 0.5)

  set_track_gap(gap = 0.005)

  circos.genomicTrack(data = do.call(rbind, lapply(seq_along(chrs), function(x) {
    CENH3_bed[CENH3_bed$chr == chrs[x] &
              CENH3_bed$start >= genomeDF$start[x] &
              CENH3_bed$end <= genomeDF$end[x],] } )),
                      panel.fun = function(region, value, ...) {
                        circos.genomicLines(region,
                                            value,
                                            col = "purple",
                                            area = TRUE,
                                            baseline = 0,
                                            border = NA,
                                            ...)
                      }, bg.border = NA)
  circos.yaxis(side = "left", at = seq(0, 4, by = 2), sector.index = get.all.sector.index()[1], labels.cex = 0.3, tick.length = convert_x(0.5, "mm"), lwd = 0.5)

  # Reset graphic parameters and internal variables
  circos.clear()
}

pdf(paste0(plotDir,
           "CEN180_frequency_per_", genomeBinName,
           "_", quantileDef, "_", quantiles, "quantiles",
           "_of_CEN180_in_T2T_Col_",
           paste0(chrName, collapse = "_"), "_circlize_zoom_v", date, ".pdf"))
circlize_plot()
draw(lgd_list2, x = unit(4, "mm"), y = unit(4, "mm"), just = c("left", "bottom"))
draw(lgd_list1, x = unit(1, "npc") - unit(2, "mm"), y = unit(4, "mm"), just = c("right", "bottom"))
dev.off()
