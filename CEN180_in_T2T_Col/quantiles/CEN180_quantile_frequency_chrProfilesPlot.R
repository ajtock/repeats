#!/applications/R/R-4.0.0/bin/Rscript

# author: Andy Tock
# contact: ajt200@cam.ac.uk
# date: 14.01.2021

# Plot windowed CEN180 frequencies for CEN180 sequences within orderingFactor quantiles

# Usage:
# /applications/R/R-4.0.0/bin/Rscript CEN180_quantile_frequency_chrProfilesPlot.R 'Chr1,Chr2,Chr3,Chr4,Chr5' wSNV 4 10000

#chrName <- unlist(strsplit("Chr1,Chr2,Chr3,Chr4,Chr5",
#                           split = ","))
#orderingFactor <- "wSNV"
#quantiles <- 4
#genomeBinSize <- 10000

args <- commandArgs(trailingOnly = T)
chrName <- unlist(strsplit(args[1],
                           split = ","))
orderingFactor <- args[2]
quantiles <- as.integer(args[3])
genomeBinSize <- as.integer(args[4])

if(floor(log10(genomeBinSize)) + 1 < 4) {
  genomeBinName <- paste0(genomeBinSize, "bp")
} else if(floor(log10(genomeBinSize)) + 1 >= 4 &
          floor(log10(genomeBinSize)) + 1 <= 6) {
  genomeBinName <- paste0(genomeBinSize/1e3, "kb")
} else if(floor(log10(genomeBinSize)) + 1 >= 7) {
  genomeBinName <- paste0(genomeBinSize/1e6, "Mb")
}

options(stringsAsFactors = F)
library(parallel)
library(GenomicRanges)

outDir <- paste0("quantiles_by_", orderingFactor, "/",
                 paste0(chrName, collapse = "_"), "/")
plotDir <- paste0(outDir, "plots/")
system(paste0("[ -d ", outDir, " ] || mkdir -p ", outDir))
system(paste0("[ -d ", plotDir, " ] || mkdir -p ", plotDir))

# Define plot titles
if(length(chrName) == 5) {
  if(grepl("_in_", orderingFactor)) {
    featureNamePlot <- paste0("Genome-wide", " ",
                              sub("_in_\\w+", "", orderingFactor), " quantiles")
  } else if(grepl("SNV", orderingFactor)) {
    featureNamePlot <- paste0("Genome-wide", " ",
                              "SNV quantiles")
  } else if(orderingFactor == "array_size") {
    featureNamePlot <- paste0("Genome-wide", " ",
                              "array-size quantiles")
  } else if(orderingFactor == "HORlengthsSum") {
    featureNamePlot <- paste0("Genome-wide", " ",
                              "activity quantiles")
  } else if(orderingFactor == "HORcount") {
    featureNamePlot <- paste0("Genome-wide", " ",
                              "HORcount quantiles")
  }
} else {
  if(grepl("_in_", orderingFactor)) {
    featureNamePlot <- paste0(paste0(chrName, collapse = "_"), " ",
                              sub("_in_\\w+", "", orderingFactor), " quantiles")
  } else if(grepl("SNV", orderingFactor)) {
    featureNamePlot <- paste0(paste0(chrName, collapse = "_"), " ",
                              "SNV quantiles")
  } else if(orderingFactor == "array_size") {
    featureNamePlot <- paste0(paste0(chrName, collapse = "_"), " ",
                              "array-size quantiles")
  } else if(orderingFactor == "HORlengthsSum") {
    featureNamePlot <- paste0(paste0(chrName, collapse = "_"), " ",
                              "activity quantiles")
  } else if(orderingFactor == "HORcount") {
    featureNamePlot <- paste0(paste0(chrName, collapse = "_"), " ",
                              "HORcount quantiles")
  }
}

# Define quantile colours
quantileColours <- c("red", "purple", "blue", "navy")

# Genomic definitions
fai <- read.table("T2T_Col.fa.fai", header = F)
chrs <- fai$V1[1:5]
chrLens <- fai$V2[1:5]

# Load table of windowed feature frequencies for each quantile
filt_quantileProfiles <- read.table(paste0(outDir,
                                           "CEN180_frequency_per_", genomeBinName,
                                           "_", quantiles, "quantiles_",
                                           "_by_", orderingFactor,
                                           "_of_CEN180_in_T2T_Col_",
                                           paste0(chrName, collapse = "_"), "_smoothed.tsv"),
                                    header = T, sep = "\t")

filt_quantileProfilesChr <- mclapply(seq_along(chrs), function(x) {
  filt_quantileProfiles[filt_quantileProfiles$chr == chrs[x],]
}, mc.cores = length(chrs))

# Load table of windowed log2(ChIP/input) coverage for CENH3 ChIP-seq
filt_covProfiles <- read.table(paste0("/home/ajt200/analysis/CENH3_seedlings_Maheshwari_Comai_2017_GenomeRes/",
                                      "snakemake_ChIPseq_T2T_Col/mapped/both/tsv/log2ChIPcontrol/",
                                      "WT_CENH3_Rep1_ChIP_SRR4430537_WT_REC8_Myc_Rep1_input_MappedOn_T2T_Col_lowXM_both_sort_norm_binSize",
                                      genomeBinName, "_smoothed.tsv"),
                               header = T, sep = "\t")
filt_covProfilesChr <- mclapply(seq_along(chrs), function(x) {
  filt_covProfiles[filt_covProfiles$chr == chrs[x],]
}, mc.cores = length(chrs))

# Get minimum and maximum levels across chrName
min_quantileProfilesChr <- min(unlist(mclapply(which(chrs %in% chrName),
  function(x) {
    min(c(filt_quantileProfilesChr[[x]][,3],
          filt_quantileProfilesChr[[x]][,4],
          filt_quantileProfilesChr[[x]][,5],
          filt_quantileProfilesChr[[x]][,6]), na.rm = T)
}, mc.cores = length(filt_quantileProfilesChr))))
max_quantileProfilesChr <- max(unlist(mclapply(which(chrs %in% chrName),
  function(x) {
    max(c(filt_quantileProfilesChr[[x]][,3],
          filt_quantileProfilesChr[[x]][,4],
          filt_quantileProfilesChr[[x]][,5],
          filt_quantileProfilesChr[[x]][,6]), na.rm = T)
}, mc.cores = length(filt_quantileProfilesChr))))

min_covProfilesChr <- min(unlist(mclapply(which(chrs %in% chrName),
  function(x) {
    min(c(filt_covProfilesChr[[x]][,6]), na.rm = T)
}, mc.cores = length(filt_covProfilesChr))))
max_covProfilesChr <- max(unlist(mclapply(which(chrs %in% chrName),
  function(x) {
    max(c(filt_covProfilesChr[[x]][,6]), na.rm = T)
}, mc.cores = length(filt_covProfilesChr))))

pdf(paste0(plotDir,
           "CEN180_frequency_per_", genomeBinName,
           "_", quantiles, "quantiles_",
           "_by_", orderingFactor,
           "_of_CEN180_in_T2T_Col_",
           paste0(chrName, collapse = "_"), ".pdf"),
           height = 4*length(chrName), width = 16)
par(mfrow = c(length(chrName), 1))
par(mar = c(5.0, 9.0, 4.0, 9.0))
for(x in which(chrs %in% chrName)) {
  par(mgp = c(3, 1, 0))
  plot(x = filt_covProfilesChr[[x]][,2],
       y = filt_covProfilesChr[[x]][,6],
       col = "grey50",
       type = "h", lwd = 1.0,
       xlim = c(0, max(chrLens[chrs %in% chrName])),
       ylim = c(-max((min_covProfilesChr*-1), max_covProfilesChr),
                max((min_covProfilesChr*-1), max_covProfilesChr)),
       xlab = "", ylab = "",
       xaxt = "n", yaxt = "n",
       main = bquote(.(featureNamePlot)), font.main = 1, cex.main = 3.0)
  axis(side = 2, cex.axis = 2.0, lwd = 2.0, lwd.tick = 2.0, col = "black", col.axis = "black", line = 0.0)
  mtext(side = 2, line = 3.0, cex = 2.0, text = bquote("Log"[2]*"(ChIP/input)"), col = "black")
  par(mgp = c(3, 1.25, 0))
  axis(side = 1, cex.axis = 2.0, lwd = 2.0, lwd.tick = 2.0,
       labels = as.character(seq(0, round(max(chrLens[chrs %in% chrName])/1e+06), by = 5)),
       at = seq(0, round(max(chrLens[chrs %in% chrName])/1e+06), by = 5)*1e+06)
  mtext(side = 1, line = 3.5, cex = 2.0, text = paste0(chrs[x], " coordinates (Mb)"), col = "black")
  par(new = T, mgp = c(3, 1.5, 0))
  plot(x = filt_quantileProfilesChr[[x]][,2],
       y = filt_quantileProfilesChr[[x]][,3],
       col = quantileColours[1],
       type = "l", lwd = 2,
       xlim = c(0, max(chrLens[chrs %in% chrName])),
       ylim = c(0-max_quantileProfilesChr, max_quantileProfilesChr),
       xlab = "", ylab = "",
       xaxt = "n", yaxt = "n")
  axis(side = 4, cex.axis = 2.0, lwd = 2.0, lwd.tick = 2.0, col = "black", col.axis = "black", line = 0.0)
  p <- par('usr')
  text(p[2], mean(p[3:4]), cex = 3.0, adj = c(0.5, -2.0), labels = "CEN180 frequency", xpd = NA, srt = -90, col = "black")
  lines(x = filt_quantileProfilesChr[[x]][,2], y = filt_quantileProfilesChr[[x]][,4], col = quantileColours[2], type = "l", lwd = 2.0)
  lines(x = filt_quantileProfilesChr[[x]][,2], y = filt_quantileProfilesChr[[x]][,5], col = quantileColours[3], type = "l", lwd = 2.0)
  lines(x = filt_quantileProfilesChr[[x]][,2], y = filt_quantileProfilesChr[[x]][,6], col = quantileColours[4], type = "l", lwd = 2.0)
  legend("topright",
         inset = c(0.15, 0.08),
         legend = c(paste0("Quantile ", 1:quantiles), "CENH3"),
         col = c("white"),
         text.col = c(quantileColours, "grey50"),
         text.font = c(1, 1),
         ncol = 1, cex = 1.5, lwd = 1.5, bty = "n")
   box(lwd = 2.0)
}
dev.off()

