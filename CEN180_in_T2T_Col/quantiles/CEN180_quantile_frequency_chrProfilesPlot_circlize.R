#!/applications/R/R-4.0.0/bin/Rscript

# author: Andy Tock
# contact: ajt200@cam.ac.uk
# date: 16.02.2021

# Plot windowed CEN180 frequencies for CEN180 sequences within orderingFactor quantiles

# Usage:
# /applications/R/R-4.0.0/bin/Rscript CEN180_quantile_frequency_chrProfilesPlot_circlize.R 'Chr1,Chr2,Chr3,Chr4,Chr5' wSNV 4 10000

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
library(circlize)

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
fai <- read.table("/home/ajt200/analysis/nanopore/T2T_Col/T2T_Col.fa.fai", header = F)
chrs <- fai$V1[1:5]
chrLens <- fai$V2[1:5]
regionGR <- GRanges(seqnames = chrs,
                    ranges = IRanges(start = rep(1, length(chrs)),
                                     end = chrLens),
                    strand = "*")
CENstart <- c(14840750, 3724530, 13597090, 4203495, 11783990)
CENend <- c(17558182, 5946091, 15733029, 6977107, 14551874)
CENGR <- GRanges(seqnames = chrs,
                 ranges = IRanges(start = CENstart,
                                  end = CENend),
                 strand = "*")

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

# Convert to BED-like format for use with circlize
filt_quantileProfiles_bed <- data.frame(chr = filt_quantileProfiles$chr,
                                        start = filt_quantileProfiles$window-1,
                                        end = filt_quantileProfiles$window-1+genomeBinSize,
                                        value1 = filt_quantileProfiles[,3],
                                        value2 = filt_quantileProfiles[,4],
                                        value3 = filt_quantileProfiles[,5],
                                        value4 = filt_quantileProfiles[,6])
filt_covProfiles_bed <- data.frame(chr = filt_covProfiles$chr,
                                   start = filt_covProfiles$window-1,
                                   end = filt_covProfiles$window-1+genomeBinSize,
                                   value1 = filt_covProfiles$filt_log2ChIPcontrol)
# Redefine end coordinate of last window to match chrLens for each chromosome
for(x in seq_along(chrs)) {
  filt_quantileProfiles_bed[filt_quantileProfiles_bed$chr == chrs[x],][dim(filt_quantileProfiles_bed[filt_quantileProfiles_bed$chr == chrs[x],])[1],]$end <- chrLens[x]
  filt_covProfiles_bed[filt_covProfiles_bed$chr == chrs[x],][dim(filt_covProfiles_bed[filt_covProfiles_bed$chr == chrs[x],])[1],]$end <- chrLens[x]
}

filt_Profiles_bedList <- list(filt_quantileProfiles_bed,
                              filt_covProfiles_bed)


## circlize

# Create df of randomly positioned centromeric loci with randomly generated values
set.seed(999)
n <- 200
df <- do.call(rbind, lapply(seq_along(chrs), function(x) {
  data.frame(sectors = sample(chrs[x], size = n, replace = T),
             x = sample(CENstart[x] : CENend[x], size = n, replace = F),
             y = runif(n))
}))

# Initialize circular layout in PDF
pdf(paste0(plotDir,
           "CEN180_frequency_per_", genomeBinName,
           "_", quantiles, "quantiles",
           "_by_", orderingFactor,
           "_of_CEN180_in_T2T_Col_",
           paste0(chrName, collapse = "_"), "_circlize.pdf"))
circos.par(track.height = 0.1,
           canvas.xlim = c(-1.1, 1.1),
           canvas.ylim = c(-1.1, 1.1),
           start.degree = 90)
#circos.initialize(sectors = df$sectors,
#                  x = df$x)
circos.initialize(sectors = rep(chrs, 2),
                  x = c(c(rep(0, 5)), chrLens))

# Add graphics in track-by-track manner using circos.trackPlotRegion() or circos.track() for short
# Similar to the "base R graphic engine, [where] you need [to] first call plot(),
# then you can use functions such as points() and lines() to add graphics."
circos.track(sectors = df$sectors,
             y = df$y,
             panel.fun = function(x, y) {
               circos.text(x = CELL_META$xcenter,
                           y = CELL_META$cell.ylim[2] + mm_y(8),
                           labels = CELL_META$sector.index)
               circos.axis(h = 1.2,
                           labels.cex = 0.6,
                           labels.niceFacing = FALSE)
             })
col = c("dodgerblue2", "orange2", "green2", "magenta2", "purple4")
circos.trackPoints(sectors = df$sectors,
                   x = df$x,
                   y = df$y,
                   col = col,
                   pch = 16,
                   cex = 0.5)
sapply(seq_along(chrs), function(x) {
  circos.text(x = (CENstart[x]+CENend[x])/2,
              y = 0.5,
              labels = paste0("CEN", x),
              sector.index = chrs[x],
              track.index = 1,
              col = "black",
              font = 4)
})
# Plot windowed CEN180 frequency for each orderingFactor quantiles
circos.genomicTrack(data = filt_quantileProfiles_bed,
                    panel.fun = function(region, value, ...) {
                      circos.genomicLines(region,
                                          value,
                                          col = quantileColours,
                                          lwd = 1,
                                          lty = 1,
                                          area = FALSE,
                                          ...)
                    })
circos.genomicTrack(data = filt_covProfiles_bed,
                    panel.fun = function(region, value, ...) {
                      circos.genomicLines(region,
                                          value,
                                          col = "grey50",
                                          area = TRUE,
                                          baseline = 0,
                                          border = NA,
                                          ...)
                    })
dev.off()
# Reset graphic parameters and internal variables
circos.clear()



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

