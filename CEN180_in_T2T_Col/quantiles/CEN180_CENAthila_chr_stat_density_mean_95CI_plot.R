#!/applications/R/R-4.0.0/bin/Rscript

# Plot density and means with 95% confidence intervals (CIs)
# of given statistic for each group of CEN180 sequences

# Usage:
# /applications/R/R-4.0.0/bin/Rscript CEN180_CENAthila_chr_stat_density_mean_95CI_plot.R 'Chr1,Chr2,Chr3,Chr4,Chr5' wSNV 4 1.00 1.0 '%1.1f' '%1.1f' '%1.1f' 0.80

#chrName <- unlist(strsplit("Chr1,Chr2,Chr3,Chr4,Chr5",
#                           split = ","))
#orderingFactor <- "wSNV"
#quantiles <- 4
#densityProp <- 1.00
#maxDensityPlus <- 1.0
#xDec <- "%1.1f"
#yDec <- "%1.1f"
#yDec2 <- "%1.1f"
#legendLabX <- 0.80

args <- commandArgs(trailingOnly = T)
chrName <- unlist(strsplit(args[1],
                           split = ","))
orderingFactor <- args[2]
quantiles <- as.numeric(args[3])
densityProp <- as.numeric(args[4])
maxDensityPlus <- as.numeric(args[5])
xDec <- as.character(args[6])
yDec <- as.character(args[7])
yDec2 <- as.character(args[8])
legendLabX <- as.numeric(args[9])

options(stringsAsFactors = F)
library(GenomicRanges)
library(parallel)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggbeeswarm)
library(ggthemes)
library(grid)
library(gridExtra)
library(extrafont)
library(scales)
library(viridis)

# Genomic definitions
fai <- read.table("/home/ajt200/analysis/nanopore/T2T_Col/T2T_Col.fa.fai", header = F)
chrs <- fai$V1[1:5]
chrLens <- fai$V2[1:5]
regionGR <- GRanges(seqnames = chrs,
                    ranges = IRanges(start = 1,
                                     end = chrLens),
                    strand = "*")
CENstart <- c(14840750, 3724530, 13597090, 4203495, 11783990)
CENend <- c(17558182, 5946091, 15733029, 6977107, 14551874)
CENGR <- GRanges(seqnames = chrs,
                 ranges = IRanges(start = CENstart,
                                  end = CENend),
                 strand = "*")

# Define chr colours
chrColours <- viridis_pal()(length(chrName))
chrColours[length(chrColours)] <- "orange"
makeTransparent <- function(thisColour, alpha = 250)
{
  newColour <- col2rgb(thisColour)
  apply(newColour, 2, function(x) {
    rgb(red = x[1], green = x[2], blue = x[3],
        alpha = alpha, maxColorValue = 255)
  })
}
chrColours <- makeTransparent(chrColours)

# Load table of features grouped into quantiles
inDir <- paste0("quantiles_by_", orderingFactor, "/",
                paste0(chrName, collapse = "_"), "/")
featuresDF <- read.table(paste0(inDir,
                                "features_", quantiles, "quantiles",
                                "_by_", orderingFactor,
                                "_of_CEN180_in_T2T_Col_",
                                paste0(chrName, collapse = "_"), ".tsv"),
                         header = T, sep = "\t", row.names = NULL)

# Load features to confirm feature (row) ordering in "featuresDF" is the same
# as in "features" (which was used for generating the coverage matrices)
features <- lapply(seq_along(chrName), function(y) {
  read.table(paste0("/home/ajt200/analysis/repeats/CEN180_in_T2T_Col/CEN180_in_T2T_Col_",
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

# Load table of centromeric Athila Gypsy LTRs
CENAthila <- lapply(seq_along(chrName), function(y) {
  read.table(paste0("/home/ajt200/analysis/repeats/CEN180_in_T2T_Col/CENAthila_in_T2T_Col_",
                    chrName[y], ".bed"),
             header = F)
})
# If features from multiple subgenomes and/or compartments are to be analysed,
# concatenate the corresponding feature data.frames
if(length(chrName) > 1) {
  CENAthila <- do.call(rbind, CENAthila)
} else {
  CENAthila <- CENAthila[[1]]
}
colnames(CENAthila) <- c("chr", "start", "end", "featureID", "score", "strand")
CENAthilaGR <- GRanges(seqnames = CENAthila$chr,
                       ranges = IRanges(CENAthila$start+1, CENAthila$end),
                       strand = CENAthila$strand,
                       featureID = CENAthila$featureID)

# Get distance between each CEN180 and the
# CENAthila that is closest or furthest away
featuresDF_dists <- data.frame()
for(i in seq_along(chrName)) {
  featuresDFchr <- featuresDF[featuresDF$chr == chrName[i],]
  CENAthilaGRchr <- CENAthilaGR[seqnames(CENAthilaGR) == chrName[i]] 

  # Calculate distances between start and end coordinates of CEN180 and CENAthila
  featuresStart_vs_CENAthilaStart <- mclapply(seq_along(featuresDFchr$start), function(x) {
    abs(featuresDFchr[x,]$start-start(CENAthilaGRchr))
  }, mc.cores = detectCores(), mc.preschedule = T)
  featuresEnd_vs_CENAthilaEnd <- mclapply(seq_along(featuresDFchr$end), function(x) {
    abs(featuresDFchr[x,]$end-end(CENAthilaGRchr))
  }, mc.cores = detectCores(), mc.preschedule = T)
  featuresStart_vs_CENAthilaEnd <- mclapply(seq_along(featuresDFchr$start), function(x) {
    abs(featuresDFchr[x,]$start-end(CENAthilaGRchr))
  }, mc.cores = detectCores(), mc.preschedule = T)
  featuresEnd_vs_CENAthilaStart <- mclapply(seq_along(featuresDFchr$end), function(x) {
    abs(featuresDFchr[x,]$end-start(CENAthilaGRchr))
  }, mc.cores = detectCores(), mc.preschedule = T)
  
  # Get distance between each CEN180 and the
  # CENAthila that is closest or furthest away
  minDists <- unlist(mclapply(seq_along(featuresDFchr$start), function(x) {
    min(c(featuresStart_vs_CENAthilaStart[[x]],
          featuresEnd_vs_CENAthilaEnd[[x]],
          featuresStart_vs_CENAthilaEnd[[x]],
          featuresEnd_vs_CENAthilaStart[[x]]))
  }, mc.cores = detectCores(), mc.preschedule = T))
  maxDists <- unlist(mclapply(seq_along(featuresDFchr$start), function(x) {
    max(c(featuresStart_vs_CENAthilaStart[[x]],
          featuresEnd_vs_CENAthilaEnd[[x]],
          featuresStart_vs_CENAthilaEnd[[x]],
          featuresEnd_vs_CENAthilaStart[[x]]))
  }, mc.cores = detectCores(), mc.preschedule = T))

  featuresDFchr <- data.frame(featuresDFchr,
                              minDistToCENAthila = minDists,
                              maxDistToCENAthila = maxDists)
  featuresDF_dists <- rbind(featuresDF_dists, featuresDFchr)
}
featuresDF <- featuresDF_dists
rm(featuresDF_dists); gc()

outDir <- paste0("correlations/", paste0(chrName, collapse = "_"), "/")

# Make scatter plots with trend lines
# minDistToCENAthila vs wSNV
ggTrend1 <- ggplot(data = featuresDF,
                   mapping = aes(x = minDistToCENAthila,
                                 y = wSNV)) +
  geom_hex(bins = 60) +
  scale_x_continuous(trans = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) +
  scale_y_continuous(trans = log2_trans(),
                     breaks = trans_breaks("log2", function(x) 2^x),
                     labels = trans_format("log2", math_format(2^.x))) +
  annotation_logticks(sides = "lb") +
  scale_fill_viridis() +
  geom_smooth(colour = "red", fill = "grey70", alpha = 0.9,
              method = "gam", formula = y~s(x)) +
  labs(x = "Distance to nearest CENAthila (bp)",
       y = "SNVs") +
  theme_bw() +
  theme(
        axis.ticks = element_line(size = 0.5, colour = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.x = element_text(size = 16, colour = "black"),
        axis.text.y = element_text(size = 16, colour = "black"),
        axis.title = element_text(size = 18, colour = "black"),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(size = 1.0, colour = "black"),
        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
        plot.title = element_text(hjust = 0.5, size = 18)) +
  ggtitle(bquote(italic(r[s]) ~ "=" ~
                 .(round(cor.test(featuresDF$minDistToCENAthila, featuresDF$wSNV, method = "spearman", use = "pairwise.complete.obs")$estimate[[1]],
                         digits = 2)) *
                 ";" ~ italic(P) ~ "=" ~
                 .(round(min(0.5, cor.test(featuresDF$minDistToCENAthila, featuresDF$wSNV, method = "spearman", use = "pairwise.complete.obs")$p.value * sqrt( (dim(featuresDF)[1]/100) )),
                         digits = 5)) ~
                 "(CEN180 in" ~ .(paste0(chrName, collapse = ",")) * ")"))

# minDistToCENAthila vs HORlengthsSum
ggTrend2 <- ggplot(data = featuresDF,
                   mapping = aes(x = minDistToCENAthila,
                                 y = HORlengthsSum+1)) +
  geom_hex(bins = 60) +
  scale_x_continuous(trans = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) +
  scale_y_continuous(trans = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks(sides = "lb") +
  scale_fill_viridis() +
  geom_smooth(colour = "red", fill = "grey70", alpha = 0.9,
              method = "gam", formula = y~s(x)) +
  labs(x = "Distance to nearest CENAthila (bp)",
       y = "Activity") +
  theme_bw() +
  theme(
        axis.ticks = element_line(size = 0.5, colour = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.x = element_text(size = 16, colour = "black"),
        axis.text.y = element_text(size = 16, colour = "black"),
        axis.title = element_text(size = 18, colour = "black"),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(size = 1.0, colour = "black"),
        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
        plot.title = element_text(hjust = 0.5, size = 18)) +
  ggtitle(bquote(italic(r[s]) ~ "=" ~
                 .(round(cor.test(featuresDF$minDistToCENAthila, featuresDF$HORlengthsSum, method = "spearman", use = "pairwise.complete.obs")$estimate[[1]],
                         digits = 2)) *
                 ";" ~ italic(P) ~ "=" ~
                 .(round(min(0.5, cor.test(featuresDF$minDistToCENAthila, featuresDF$HORlengthsSum, method = "spearman", use = "pairwise.complete.obs")$p.value * sqrt( (dim(featuresDF)[1]/100) )),
                         digits = 5)) ~
                 "(CEN180 in" ~ .(paste0(chrName, collapse = ",")) * ")"))

# minDistToCENAthila vs HORcount
ggTrend3 <- ggplot(data = featuresDF,
                   mapping = aes(x = minDistToCENAthila,
                                 y = HORcount+1)) +
  geom_hex(bins = 60) +
  scale_x_continuous(trans = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) +
  scale_y_continuous(trans = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks(sides = "lb") +
  scale_fill_viridis() +
  geom_smooth(colour = "red", fill = "grey70", alpha = 0.9,
              method = "gam", formula = y~s(x)) +
  labs(x = "Distance to nearest CENAthila (bp)",
       y = "HOR count") +
  theme_bw() +
  theme(
        axis.ticks = element_line(size = 0.5, colour = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.x = element_text(size = 16, colour = "black"),
        axis.text.y = element_text(size = 16, colour = "black"),
        axis.title = element_text(size = 18, colour = "black"),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(size = 1.0, colour = "black"),
        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
        plot.title = element_text(hjust = 0.5, size = 18)) +
  ggtitle(bquote(italic(r[s]) ~ "=" ~
                 .(round(cor.test(featuresDF$minDistToCENAthila, featuresDF$HORcount, method = "spearman", use = "pairwise.complete.obs")$estimate[[1]],
                         digits = 2)) *
                 ";" ~ italic(P) ~ "=" ~
                 .(round(min(0.5, cor.test(featuresDF$minDistToCENAthila, featuresDF$HORcount, method = "spearman", use = "pairwise.complete.obs")$p.value * sqrt( (dim(featuresDF)[1]/100) )),
                         digits = 5)) ~
                 "(CEN180 in" ~ .(paste0(chrName, collapse = ",")) * ")"))

# minDistToCENAthila vs CENH3_in_bodies
ggTrend4 <- ggplot(data = featuresDF,
                   mapping = aes(x = minDistToCENAthila,
                                 y = CENH3_in_bodies)) +
  geom_hex(bins = 60) +
  scale_x_continuous(trans = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks(sides = "b") +
  scale_fill_viridis() +
  geom_smooth(colour = "red", fill = "grey70", alpha = 0.9,
              method = "gam", formula = y~s(x)) +
  labs(x = "Distance to nearest CENAthila (bp)",
       y = "CENH3") +
  theme_bw() +
  theme(
        axis.ticks = element_line(size = 0.5, colour = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.x = element_text(size = 16, colour = "black"),
        axis.text.y = element_text(size = 16, colour = "black"),
        axis.title = element_text(size = 18, colour = "black"),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(size = 1.0, colour = "black"),
        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
        plot.title = element_text(hjust = 0.5, size = 18)) +
  ggtitle(bquote(italic(r[s]) ~ "=" ~
                 .(round(cor.test(featuresDF$minDistToCENAthila, featuresDF$CENH3_in_bodies, method = "spearman", use = "pairwise.complete.obs")$estimate[[1]],
                         digits = 2)) *
                 ";" ~ italic(P) ~ "=" ~
                 .(round(min(0.5, cor.test(featuresDF$minDistToCENAthila, featuresDF$CENH3_in_bodies, method = "spearman", use = "pairwise.complete.obs")$p.value * sqrt( (dim(featuresDF)[1]/100) )),
                         digits = 5)) ~
                 "(CEN180 in" ~ .(paste0(chrName, collapse = ",")) * ")"))

ggTrend_combined <- grid.arrange(grobs = list(
                                              ggTrend1,
                                              ggTrend2,
                                              ggTrend3,
                                              ggTrend4
                                             ),
                                 layout_matrix = rbind(
                                                       1:4
                                                      ))
ggsave(paste0(outDir,
              "trends_for_CEN180_distance_to_nearest_CENAthila_in_T2T_Col_",
              paste0(chrName, collapse = "_"), ".pdf"),
       plot = ggTrend_combined, height = 7*1, width = 8*4)

if(length(chrName) > 1) {

# Get row indices for each feature chr
chrIndices <- lapply(seq_along(chrName), function(k) {
  which(featuresDF$chr == chrName[k])
})

# Calculate means, SDs, SEMs and 95% CIs
# and create dataframe of summary statistics for plotting
featureNamePlot <- "CEN180"
parName <- colnames(featuresDF)[c(5, 7, 8, 12, 46)]
print(parName)
parNamePlot <- c("SNVs", "Activity", "HOR count", "CENH3", "Distance to nearest CENAthila (bp)")
for(i in seq_along(parName)) {
  featuresDF_chrMean <- sapply(seq_along(chrName), function(k) {
    mean(featuresDF[featuresDF$chr == chrName[k],][,colnames(featuresDF) == parName[i],], na.rm = T)
  })
  featuresDF_chrSD <- sapply(seq_along(chrName), function(k) {
    sd(featuresDF[featuresDF$chr == chrName[k],][,colnames(featuresDF) == parName[i],], na.rm = T)
  })
  featuresDF_chrSEM <- sapply(seq_along(chrName), function(k) {
    featuresDF_chrSD[k] / sqrt( (dim(featuresDF[featuresDF$chr == chrName[k],])[1] - 1) )
  })
  featuresDF_chrCIlower <- sapply(seq_along(chrName), function(k) {
    featuresDF_chrMean[k] -
      ( qt(0.975, df = dim(featuresDF[featuresDF$chr == chrName[k],])[1]-1 ) *
        featuresDF_chrSEM[k] )
  })
  featuresDF_chrCIupper <- sapply(seq_along(chrName), function(k) {
    featuresDF_chrMean[k] +
      ( qt(0.975, df = dim(featuresDF[featuresDF$chr == chrName[k],])[1]-1 ) *
        featuresDF_chrSEM[k] )
  })
  featuresDF_summary_stats <- data.frame(chr = chrName,
                                         Mean = featuresDF_chrMean,
                                         SD = featuresDF_chrSD,
                                         SEM = featuresDF_chrSEM,
                                         CIlower = featuresDF_chrCIlower,
                                         CIupper = featuresDF_chrCIupper)
  
  summary_stats_min <- min(c(featuresDF_summary_stats$CIlower), na.rm = T)
  summary_stats_max <- max(c(featuresDF_summary_stats$CIupper), na.rm = T)
  
  featuresDFq <- featuresDF[which(featuresDF[,which(colnames(featuresDF) == parName[i])] <=
                                  quantile(featuresDF[,which(colnames(featuresDF) == parName[i])],
                                           probs = densityProp, na.rm = T)),]
  if(parName[i] %in% c("HORlengthsSum", "HORcount")) { 
    featuresDFq[,which(colnames(featuresDFq) == parName[i])] <- featuresDFq[,which(colnames(featuresDFq) == parName[i])]+1
  }
  xmin <- min(c(featuresDFq[,which(colnames(featuresDFq) == parName[i])]),
                na.rm = T)
  xmax <- max(c(featuresDFq[,which(colnames(featuresDFq) == parName[i])]),
                na.rm = T)
  minDensity <- 0
  maxDensity <- max(
    c(
      sapply(seq_along(chrName), function(k) {
        max(density(featuresDFq[featuresDFq$chr == chrName[k],][,which(colnames(featuresDFq) == parName[i])],
            na.rm = T)$y)
      })
     )
  )+maxDensityPlus
  
  # Define legend labels
  legendLabs_feature <- lapply(seq_along(chrName), function(x) {
    grobTree(textGrob(bquote(.(chrName[x])),
                      x = legendLabX, y = 0.90-((x-1)*0.07), just = "left",
                      gp = gpar(col = chrColours[x], fontsize = 22)))
  })
  
  # Parameter density plot function
  popgen_stats_plotFun <- function(lociDF,
                                   parameter,
                                   parameterLab,
                                   featureGroup,
                                   featureNamePlot,
                                   legendLabs,
                                   chrColours) {
    ggplot(data = lociDF,
           mapping = aes(x = get(parameter),
                         colour = reorder(x = get(featureGroup), X = desc(get(featureGroup))),
                         group = reorder(x = get(featureGroup), X = desc(get(featureGroup))))) +
    scale_colour_manual(values = rev(chrColours)) +
    geom_density(size = 1.5) +
#    scale_y_continuous(limits = c(minDensity, maxDensity),
    scale_y_continuous(
                       labels = function(x) sprintf(yDec, x)) +
    labs(x = parameterLab,
         y = "Density") +
    annotation_custom(legendLabs[[1]]) +
    annotation_custom(legendLabs[[2]]) +
    annotation_custom(legendLabs[[3]]) +
    annotation_custom(legendLabs[[4]]) +
    annotation_custom(legendLabs[[5]]) +
    theme_bw() +
    theme(axis.line.y = element_line(size = 1.0, colour = "black"),
          axis.line.x = element_line(size = 1.0, colour = "black"),
          axis.ticks.y = element_line(size = 0.5, colour = "black"),
          axis.ticks.x = element_line(size = 0.5, colour = "black"),
          axis.ticks.length = unit(0.25, "cm"),
          axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
          axis.text.x = element_text(size = 18, colour = "black", family = "Luxi Mono"),
          axis.title = element_text(size = 26, colour = "black"),
          legend.position = "none",
          panel.grid = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          plot.margin = unit(c(0.3,1.2,0.3,0.3),"cm"),
          plot.title = element_text(hjust = 0.5, size = 30)) +
    ggtitle(bquote(.(featureNamePlot)))
  }
  
  # Plot means and 95% confidence intervals
  popgen_stats_meanCIs <- function(dataFrame,
                                   parameterLab,
                                   featureGroup,
                                   featureNamePlot,
                                   chrColours) {
    ggplot(data = dataFrame,
           mapping = aes(x = get(featureGroup),
                         y = Mean,
                         colour = get(featureGroup))) +
    labs(colour = "") +
    geom_point(shape = 19, size = 6, position = position_dodge(width = 0.2)) +
    geom_errorbar(mapping = aes(ymin = CIlower,
                                ymax = CIupper),
                  width = 0.2, size = 2, position = position_dodge(width = 0.2)) +
    scale_colour_manual(values = chrColours) +
#    scale_y_continuous(trans = log10_trans(),
#                       limits = c(summary_stats_min, summary_stats_max),
#                       breaks = trans_breaks("log10", function(x) 10^x),
#                       labels = trans_format("log10", math_format(10^.x))) +
#    annotation_logticks(sides = "l") +
    scale_y_continuous(limits = c(summary_stats_min, summary_stats_max),
                       labels = function(x) sprintf(yDec2, x)) +
    labs(x = "",
         y = parameterLab) +
    theme_bw() +
    theme(axis.line.y = element_line(size = 1.0, colour = "black"),
          axis.ticks.y = element_line(size = 0.5, colour = "black"),
          axis.ticks.x = element_blank(),
          axis.ticks.length = unit(0.25, "cm"),
          axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
          axis.text.x = element_text(size = 22, colour = chrColours, hjust = 1.0, vjust = 1.0, angle = 45),
          axis.title = element_text(size = 26, colour = "black"),
          legend.position = "none",
          panel.grid = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          plot.margin = unit(c(0.3,1.2,0.1,0.3),"cm"),
          plot.title = element_text(hjust = 0.5, size = 30)) +
    ggtitle(bquote(.(featureNamePlot)))
  }
  
  ggObjGA_feature <- popgen_stats_plotFun(lociDF = featuresDFq,
                                          parameter = parName[i],
                                          parameterLab = bquote(.(parNamePlot[i])),
                                          featureGroup = "chr", 
                                          featureNamePlot = featureNamePlot,
                                          legendLabs = legendLabs_feature,
                                          chrColours = chrColours
                                         )
  if(parName[i] %in% c("HORlengthsSum", "HORcount", "minDistToCENAthila")) {
    ggObjGA_feature <- ggObjGA_feature + scale_x_continuous(trans = log10_trans(),
                                                            breaks = trans_breaks("log10", function(x) 10^x),
                                                            labels = trans_format("log10", math_format(10^.x))) +
                                                            annotation_logticks(sides = "b")
  } else if(parName[i] %in% c("wSNV")) {
    ggObjGA_feature <- ggObjGA_feature + scale_x_continuous(trans = log2_trans(),
                                                            breaks = trans_breaks("log2", function(x) 2^x),
                                                            labels = trans_format("log2", math_format(2^.x))) +
                                                            annotation_logticks(sides = "b")
 } else {
   ggObjGA_feature <- ggObjGA_feature + scale_x_continuous(limits = c(xmin, xmax),
                                                           labels = function(x) sprintf(xDec, x))
 }
  ggObjGA_feature_mean <- popgen_stats_meanCIs(dataFrame = featuresDF_summary_stats,
                                               parameterLab = bquote(.(parNamePlot[i])),
                                               featureGroup = "chr",
                                               featureNamePlot = featureNamePlot,
                                               chrColours = chrColours
                                              )
  ggObjGA_combined <- grid.arrange(ggObjGA_feature,
                                   ggObjGA_feature_mean,
                                   ncol = 1, as.table = F)
  ggsave(paste0(outDir,
                parName[i], "_densityProp", densityProp, "_at_",
                "_CEN180_in_T2T_Col_",
                paste0(chrName, collapse = "_"), ".pdf"),
         plot = ggObjGA_combined,
         height = 13, width = 7)
}

}
