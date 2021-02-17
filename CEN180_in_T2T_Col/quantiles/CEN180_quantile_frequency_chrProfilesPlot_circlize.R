#!/applications/R/R-4.0.0/bin/Rscript

# author: Andy Tock
# contact: ajt200@cam.ac.uk
# date: 16.02.2021

# Plot windowed CEN180 frequencies for CEN180 sequences within orderingFactor quantiles

# Usage:
# /applications/R/R-4.0.0/bin/Rscript CEN180_quantile_frequency_chrProfilesPlot_circlize.R 'Chr1,Chr2,Chr3,Chr4,Chr5' 4 10000 perchrom

#chrName <- unlist(strsplit("Chr1,Chr2,Chr3,Chr4,Chr5",
#                           split = ","))
#quantiles <- 4
#genomeBinSize <- 10000
#quantileDef <- "genomewide"
#quantileDef <- "perchrom"

args <- commandArgs(trailingOnly = T)
chrName <- unlist(strsplit(args[1],
                           split = ","))
quantiles <- as.integer(args[2])
genomeBinSize <- as.integer(args[3])
quantileDef <- args[4]

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

plotDir <- paste0("plots/")
system(paste0("[ -d ", plotDir, " ] || mkdir -p ", plotDir))

# Define quantile colours
quantileColours <- c("red", "purple", "blue", "navy")

# Genomic definitions
fai <- read.table("/home/ajt200/analysis/nanopore/T2T_Col/T2T_Col.fa.fai", header = F)
chrs <- fai$V1[seq_along(chrName)]
chrLens <- fai$V2[seq_along(chrName)]
regionGR <- GRanges(seqnames = chrs,
                    ranges = IRanges(start = rep(1, length(chrs)),
                                     end = chrLens),
                    strand = "*")
CENstart <- c(14840750, 3724530, 13597090, 4203495, 11783990)[seq_along(chrName)]
CENend <- c(17558182, 5946091, 15733029, 6977107, 14551874)[seq_along(chrName)]
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
DNAmethPE <- read.table(paste0("/home/ajt200/analysis/BSseq_seedling_Yang_Zhu_2016_CellRes/",
                               "snakemake_BSseq_T2T_Col/coverage/tsv/",
                               "DNAmeth_Col0_BSseq_Rep1_MappedOn_T2T_Col_dedup_binSize",
                               "200kb_smoothed.tsv"),
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
DNAmethPE_bed <- data.frame(chr = DNAmethPE$chr,
                            start = DNAmethPE$window-1,
                            end = DNAmethPE$window-1+200000,
                            value1 = DNAmethPE[,4],
                            value2 = DNAmethPE[,5],
                            value3 = DNAmethPE[,6])

# Redefine end coordinate of last window to match chrLens for each chromosome
for(x in seq_along(chrs)) {
  CENH3_in_bodies_CEN180_bed[CENH3_in_bodies_CEN180_bed$chr == chrs[x],][dim(CENH3_in_bodies_CEN180_bed[CENH3_in_bodies_CEN180_bed$chr == chrs[x],])[1],]$end <- chrLens[x]
  HORlengthsSum_CEN180_bed[HORlengthsSum_CEN180_bed$chr == chrs[x],][dim(HORlengthsSum_CEN180_bed[HORlengthsSum_CEN180_bed$chr == chrs[x],])[1],]$end <- chrLens[x]
  wSNV_CEN180_bed[wSNV_CEN180_bed$chr == chrs[x],][dim(wSNV_CEN180_bed[wSNV_CEN180_bed$chr == chrs[x],])[1],]$end <- chrLens[x]
  CENH3_bed[CENH3_bed$chr == chrs[x],][dim(CENH3_bed[CENH3_bed$chr == chrs[x],])[1],]$end <- chrLens[x]
  DNAmethPE_bed[DNAmethPE_bed$chr == chrs[x],][dim(DNAmethPE_bed[DNAmethPE_bed$chr == chrs[x],])[1],]$end <- chrLens[x]
}


## circlize

# Create df of randomly positioned centromeric loci with randomly generated values
set.seed(999)
n <- 200
df <- do.call(rbind, lapply(seq_along(chrs), function(x) {
  data.frame(sectors = sample(chrs[x], size = n, replace = T),
             x = sample(CENstart[x] : CENend[x], size = n, replace = F),
             y = runif(n))
}))

genomeDF <- data.frame(chr = chrs,
                       start = rep(0, length(chrs)),

                       end = chrLens)
#circos.initialize(sectors = df$sectors,
#                  x = df$x)
#circos.initialize(sectors = rep(chrs, 2),
#                  x = c(c(rep(0, 5)), chrLens))

# Initialize circular layout in PDF
pdf(paste0(plotDir,
           "CEN180_frequency_per_", genomeBinName,
           "_", quantiles, "quantiles",
           "_of_CEN180_in_T2T_Col_",
           paste0(chrName, collapse = "_"), "_circlize_", quantileDef, "_test.pdf"))
circos.par(track.height = 0.1,
           canvas.xlim = c(-1.1, 1.1),
           canvas.ylim = c(-1.1, 1.1),
           start.degree = 90)
circos.genomicInitialize(data = genomeDF,
                         plotType = NULL,
                         tickLabelsStartFromZero = TRUE)
circos.track(ylim = c(0, 1),
             bg.col = "grey80",
             bg.border = NA,
             track.height = 0.05,
             panel.fun = function(x, y) {
               circos.genomicAxis(h = "top",
                                  labels.facing = "clockwise")
             })
sapply(seq_along(chrs), function(x) {
  circos.text(x = (CENstart[x]+CENend[x])/2,
              y = 0.5,
              labels = paste0("CEN", x),
              sector.index = chrs[x],
              track.index = 1,
              cex = 0.8,
              col = "white",
              font = 4)
})

# Add graphics in track-by-track manner using circos.trackPlotRegion() or circos.track() for short
# Similar to the "base R graphic engine, [where] you need [to] first call plot(),
# then you can use functions such as points() and lines() to add graphics."
#circos.track(sectors = df$sectors,
#             y = df$y,
#             panel.fun = function(x, y) {
#               circos.text(x = CELL_META$xcenter,
#                           y = CELL_META$cell.ylim[2] + mm_y(8),
#                           labels = CELL_META$sector.index)
#               circos.axis(h = 1.3,
#                           labels.cex = 0.6,
#                           labels.niceFacing = FALSE)
#             })
#col = c("dodgerblue2", "orange2", "green2", "magenta2", "purple4")
#circos.trackPoints(sectors = df$sectors,
#                   x = df$x,
#                   y = df$y,
#                   col = col,
#                   pch = 16,
#                   cex = 0.5)
#sapply(seq_along(chrs), function(x) {
#  circos.text(x = (CENstart[x]+CENend[x])/2,
#              y = 0.5,
#              labels = paste0("CEN", x),
#              sector.index = chrs[x],
#              track.index = 1,
#              col = "white",
#              font = 4)
#})

# Plot windowed CEN180 frequency for each quantile
circos.genomicTrack(data = CENH3_in_bodies_CEN180_bed,
                    panel.fun = function(region, value, ...) {
                      circos.genomicLines(region,
                                          value,
                                          col = quantileColours,
                                          lwd = 1,
                                          lty = 1,
                                          area = FALSE,
                                          ...)
                    })
circos.genomicTrack(data = HORlengthsSum_CEN180_bed,
                    panel.fun = function(region, value, ...) {
                      circos.genomicLines(region,
                                          value,
                                          col = quantileColours,
                                          lwd = 1,
                                          lty = 1,
                                          area = FALSE,
                                          ...)
                    })
circos.genomicTrack(data = wSNV_CEN180_bed,
                    panel.fun = function(region, value, ...) {
                      circos.genomicLines(region,
                                          value,
                                          col = quantileColours,
                                          lwd = 1,
                                          lty = 1,
                                          area = FALSE,
                                          ...)
                    })
circos.genomicTrack(data = CENH3_bed,
                    panel.fun = function(region, value, ...) {
                      circos.genomicLines(region,
                                          value,
                                          col = "darkorange1",
                                          area = TRUE,
                                          baseline = 0,
                                          border = NA,
                                          ...)
                    })
circos.genomicTrack(data = DNAmethPE_bed, numeric.column = 4, 
                    panel.fun = function(region, value, ...) {
                      circos.genomicLines(region,
                                          value,
                                          col = "dodgerblue4",
                                          area = TRUE,
                                          baseline = 0,
                                          border = NA,
                                          ...)
                    })
circos.genomicTrack(data = DNAmethPE_bed, numeric.column = 5, 
                    panel.fun = function(region, value, ...) {
                      circos.genomicLines(region,
                                          value,
                                          col = "dodgerblue1",
                                          area = TRUE,
                                          baseline = 0,
                                          border = NA,
                                          ...)
                    })
circos.genomicTrack(data = DNAmethPE_bed, numeric.column = 6, 
                    panel.fun = function(region, value, ...) {
                      circos.genomicLines(region,
                                          value,
                                          col = "cyan2",
                                          area = TRUE,
                                          baseline = 0,
                                          border = NA,
                                          ...)
                    })
dev.off()
# Reset graphic parameters and internal variables
circos.clear()



# Get minimum and maximum levels across chrName
min_wSNV_CEN180Chr <- min(unlist(mclapply(which(chrs %in% chrName),
  function(x) {
    min(c(wSNV_CEN180Chr[[x]][,3],
          wSNV_CEN180Chr[[x]][,4],
          wSNV_CEN180Chr[[x]][,5],
          wSNV_CEN180Chr[[x]][,6]), na.rm = T)
}, mc.cores = length(wSNV_CEN180Chr))))
max_wSNV_CEN180Chr <- max(unlist(mclapply(which(chrs %in% chrName),
  function(x) {
    max(c(wSNV_CEN180Chr[[x]][,3],
          wSNV_CEN180Chr[[x]][,4],
          wSNV_CEN180Chr[[x]][,5],
          wSNV_CEN180Chr[[x]][,6]), na.rm = T)
}, mc.cores = length(wSNV_CEN180Chr))))

min_CENH3Chr <- min(unlist(mclapply(which(chrs %in% chrName),
  function(x) {
    min(c(CENH3Chr[[x]][,6]), na.rm = T)
}, mc.cores = length(CENH3Chr))))
max_CENH3Chr <- max(unlist(mclapply(which(chrs %in% chrName),
  function(x) {
    max(c(CENH3Chr[[x]][,6]), na.rm = T)
}, mc.cores = length(CENH3Chr))))

pdf(paste0(plotDir,
           "CEN180_frequency_per_", genomeBinName,
           "_", quantiles, "quantiles_",
           "_by_wSNV",
           "_of_CEN180_in_T2T_Col_",
           paste0(chrName, collapse = "_"), ".pdf"),
           height = 4*length(chrName), width = 16)
par(mfrow = c(length(chrName), 1))
par(mar = c(5.0, 9.0, 4.0, 9.0))
for(x in which(chrs %in% chrName)) {
  par(mgp = c(3, 1, 0))
  plot(x = CENH3Chr[[x]][,2],
       y = CENH3Chr[[x]][,6],
       col = "grey50",
       type = "h", lwd = 1.0,
       xlim = c(0, max(chrLens[chrs %in% chrName])),
       ylim = c(-max((min_CENH3Chr*-1), max_CENH3Chr),
                max((min_CENH3Chr*-1), max_CENH3Chr)),
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
  plot(x = wSNV_CEN180Chr[[x]][,2],
       y = wSNV_CEN180Chr[[x]][,3],
       col = quantileColours[1],
       type = "l", lwd = 2,
       xlim = c(0, max(chrLens[chrs %in% chrName])),
       ylim = c(0-max_wSNV_CEN180Chr, max_wSNV_CEN180Chr),
       xlab = "", ylab = "",
       xaxt = "n", yaxt = "n")
  axis(side = 4, cex.axis = 2.0, lwd = 2.0, lwd.tick = 2.0, col = "black", col.axis = "black", line = 0.0)
  p <- par('usr')
  text(p[2], mean(p[3:4]), cex = 3.0, adj = c(0.5, -2.0), labels = "CEN180 frequency", xpd = NA, srt = -90, col = "black")
  lines(x = wSNV_CEN180Chr[[x]][,2], y = wSNV_CEN180Chr[[x]][,4], col = quantileColours[2], type = "l", lwd = 2.0)
  lines(x = wSNV_CEN180Chr[[x]][,2], y = wSNV_CEN180Chr[[x]][,5], col = quantileColours[3], type = "l", lwd = 2.0)
  lines(x = wSNV_CEN180Chr[[x]][,2], y = wSNV_CEN180Chr[[x]][,6], col = quantileColours[4], type = "l", lwd = 2.0)
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

