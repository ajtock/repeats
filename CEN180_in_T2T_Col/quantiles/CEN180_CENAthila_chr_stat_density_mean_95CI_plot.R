#!/applications/R/R-4.0.0/bin/Rscript

# Plot density and means with 95% confidence intervals (CIs)
# of given statistic for each group of CEN180 sequences

# Usage:
# /applications/R/R-3.5.0/bin/Rscript CEN180_CENAthila_chr_stat_density_mean_95CI_plot.R 'Chr1,Chr2,Chr3,Chr4,Chr5' wSNV 4 1.00 0.2 '%2.0f' '%3.1f' '%2.0f' 0.65

#chrName <- unlist(strsplit("Chr1,Chr2,Chr3,Chr4,Chr5",
#                           split = ","))
#orderingFactor <- "wSNV"
#quantiles <- 4
#densityProp <- 0.99
#maxDensityPlus <- 0.2
#xDec <- "%1.1f"
#yDec <- "%1.1f"
#yDec2 <- "%1.2f"
#legendLabX <- 0.65

args <- commandArgs(trailingOnly = T)
chrName <- unlist(strsplit(args[3],
                           split = ","))
orderingFactor <- args[4]
quantiles <- as.numeric(args[5])
densityProp <- as.numeric(args[6])
maxDensityPlus <- as.numeric(args[7])
xDec <- as.character(args[8])
yDec <- as.character(args[9])
yDec2 <- as.character(args[10])
legendLabX <- as.numeric(args[11])

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
chrColours <- viridis_pal()(length(chrs))
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
for(i in seq_along(chrs)) {
  featuresDFchr <- featuresDF[featuresDF$chr == chrs[i],]
  CENAthilaGRchr <- CENAthilaGR[seqnames(CENAthilaGR) == chrs[i]] 

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
  geom_hex(bins = 40) +
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
  labs(x = "Distance to nearest CENAthila",
       y = "SNVs") +
  theme_bw() +
  theme(
        axis.ticks = element_line(size = 1.0, colour = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.x = element_text(size = 16, colour = "black"),
        axis.text.y = element_text(size = 16, colour = "black"),
        axis.title = element_text(size = 18, colour = "black"),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(size = 3.5, colour = "black"),
        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
        plot.title = element_text(hjust = 0.5, size = 18)) +
  ggtitle(bquote(italic(r[s]) ~ "=" ~
                 .(round(cor.test(featuresDF$minDistToCENAthila, featuresDF$wSNV, method = "spearman", use = "pairwise.complete.obs")$estimate[[1]],
                         digits = 2)) *
                 ";" ~ italic(P) ~ "=" ~
                 .(round(min(0.5, cor.test(featuresDF$minDistToCENAthila, featuresDF$wSNV, method = "spearman", use = "pairwise.complete.obs")$p.value * sqrt( (dim(featuresDF)[1]/100) )),
                         digits = 5)) ~
                 "(T2T_Col" ~ .(paste0(chrName, collapse = ",")) ~ ")"))

# minDistToCENAthila vs HORlengthsSum
ggTrend2 <- ggplot(data = featuresDF,
                   mapping = aes(x = minDistToCENAthila,
                                 y = HORlengthsSum)) +
  geom_hex(bins = 40) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks(sides = "lb") +
  scale_fill_viridis() +
  geom_smooth(colour = "red", fill = "grey70", alpha = 0.9,
              method = "gam", formula = y~s(x)) +
  labs(x = "Distance to nearest CENAthila",
       y = "Activity") +
  theme_bw() +
  theme(
        axis.ticks = element_line(size = 1.0, colour = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.x = element_text(size = 16, colour = "black"),
        axis.text.y = element_text(size = 16, colour = "black"),
        axis.title = element_text(size = 18, colour = "black"),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(size = 3.5, colour = "black"),
        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
        plot.title = element_text(hjust = 0.5, size = 18)) +
  ggtitle(bquote(italic(r[s]) ~ "=" ~
                 .(round(cor.test(featuresDF$minDistToCENAthila, featuresDF$HORlengthsSum, method = "spearman", use = "pairwise.complete.obs")$estimate[[1]],
                         digits = 2)) *
                 ";" ~ italic(P) ~ "=" ~
                 .(round(min(0.5, cor.test(featuresDF$minDistToCENAthila, featuresDF$HORlengthsSum, method = "spearman", use = "pairwise.complete.obs")$p.value * sqrt( (dim(featuresDF)[1]/100) )),
                         digits = 5)) ~
                 "(T2T_Col" ~ .(paste0(chrName, collapse = ",")) ~ ")"))

# minDistToCENAthila vs HORcount
ggTrend3 <- ggplot(data = featuresDF,
                   mapping = aes(x = minDistToCENAthila,
                                 y = HORcount)) +
  geom_hex(bins = 40) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks(sides = "lb") +
  scale_fill_viridis() +
  geom_smooth(colour = "red", fill = "grey70", alpha = 0.9,
              method = "gam", formula = y~s(x)) +
  labs(x = "Distance to nearest CENAthila",
       y = "Activity") +
  theme_bw() +
  theme(
        axis.ticks = element_line(size = 1.0, colour = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.x = element_text(size = 16, colour = "black"),
        axis.text.y = element_text(size = 16, colour = "black"),
        axis.title = element_text(size = 18, colour = "black"),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(size = 3.5, colour = "black"),
        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
        plot.title = element_text(hjust = 0.5, size = 18)) +
  ggtitle(bquote(italic(r[s]) ~ "=" ~
                 .(round(cor.test(featuresDF$minDistToCENAthila, featuresDF$HORcount, method = "spearman", use = "pairwise.complete.obs")$estimate[[1]],
                         digits = 2)) *
                 ";" ~ italic(P) ~ "=" ~
                 .(round(min(0.5, cor.test(featuresDF$minDistToCENAthila, featuresDF$HORcount, method = "spearman", use = "pairwise.complete.obs")$p.value * sqrt( (dim(featuresDF)[1]/100) )),
                         digits = 5)) ~
                 "(T2T_Col" ~ .(paste0(chrName, collapse = ",")) ~ ")"))

# minDistToCENAthila vs CENH3_in_bodies
ggTrend4 <- ggplot(data = featuresDF,
                   mapping = aes(x = minDistToCENAthila,
                                 y = CENH3_in_bodies)) +
  geom_hex(bins = 40) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks(sides = "lb") +
  scale_fill_viridis() +
  geom_smooth(colour = "red", fill = "grey70", alpha = 0.9,
              method = "gam", formula = y~s(x)) +
  labs(x = "Distance to nearest CENAthila",
       y = "CENH3") +
  theme_bw() +
  theme(
        axis.ticks = element_line(size = 1.0, colour = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.x = element_text(size = 16, colour = "black"),
        axis.text.y = element_text(size = 16, colour = "black"),
        axis.title = element_text(size = 18, colour = "black"),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(size = 3.5, colour = "black"),
        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
        plot.title = element_text(hjust = 0.5, size = 18)) +
  ggtitle(bquote(italic(r[s]) ~ "=" ~
                 .(round(cor.test(featuresDF$minDistToCENAthila, featuresDF$CENH3_in_bodies, method = "spearman", use = "pairwise.complete.obs")$estimate[[1]],
                         digits = 2)) *
                 ";" ~ italic(P) ~ "=" ~
                 .(round(min(0.5, cor.test(featuresDF$minDistToCENAthila, featuresDF$CENH3_in_bodies, method = "spearman", use = "pairwise.complete.obs")$p.value * sqrt( (dim(featuresDF)[1]/100) )),
                         digits = 5)) ~
                 "(T2T_Col" ~ .(paste0(chrName, collapse = ",")) ~ ")"))

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
              "trends_distance_to_nearest_CENAthila_for_CEN180_in_T2T_Col_",
              paste0(chrName, collapse = "_"), ".pdf"),
       plot = ggTrend_combined, height = 7*1, width = 8*4)


# Get row indices for each feature quantile
chrIndices <- lapply(1:quantiles, function(k) {
  which(featuresDF$chr == paste0("Quantile ", k))
})


# Calculate means, SDs, SEMs and 95% CIs
# and create dataframe of summary statistics for plotting
featuresDF_quantileMean <- sapply(1:quantiles, function(k) {
  mean(featuresDF[featuresDF$quantile == paste0("Quantile ", k),][,colnames(featuresDF) == parName,], na.rm = T)
})
featuresDF_quantileSD <- sapply(1:quantiles, function(k) {
  sd(featuresDF[featuresDF$quantile == paste0("Quantile ", k),][,colnames(featuresDF) == parName,], na.rm = T)
})
featuresDF_quantileSEM <- sapply(1:quantiles, function(k) {
  featuresDF_quantileSD[k] / sqrt( (dim(featuresDF[featuresDF$quantile == paste0("Quantile ", k),])[1] - 1) )
})
featuresDF_quantileCIlower <- sapply(1:quantiles, function(k) {
  featuresDF_quantileMean[k] -
    ( qt(0.975, df = dim(featuresDF[featuresDF$quantile == paste0("Quantile ", k),])[1]-1 ) *
      featuresDF_quantileSEM[k] )
})
featuresDF_quantileCIupper <- sapply(1:quantiles, function(k) {
  featuresDF_quantileMean[k] +
    ( qt(0.975, df = dim(featuresDF[featuresDF$quantile == paste0("Quantile ", k),])[1]-1 ) *
      featuresDF_quantileSEM[k] )
})
featuresDF_summary_stats <- data.frame(quantile = paste0("Quantile ", 1:quantiles),
                                       Mean = featuresDF_quantileMean,
                                       SD = featuresDF_quantileSD,
                                       SEM = featuresDF_quantileSEM,
                                       CIlower = featuresDF_quantileCIlower,
                                       CIupper = featuresDF_quantileCIupper,
                                       stringsAsFactors = F)

ranFeatsDF_randomMean <- sapply(1:quantiles, function(k) {
  mean(ranFeatsDF[ranFeatsDF$random == paste0("Random ", k),][,colnames(ranFeatsDF) == parName,], na.rm = T)
})
ranFeatsDF_randomSD <- sapply(1:quantiles, function(k) {
  sd(ranFeatsDF[ranFeatsDF$random == paste0("Random ", k),][,colnames(ranFeatsDF) == parName,], na.rm = T)
})
ranFeatsDF_randomSEM <- sapply(1:quantiles, function(k) {
  ranFeatsDF_randomSD[k] / sqrt( (dim(ranFeatsDF[ranFeatsDF$random == paste0("Random ", k),])[1] - 1) )
})
ranFeatsDF_randomCIlower <- sapply(1:quantiles, function(k) {
  ranFeatsDF_randomMean[k] -
    ( qt(0.975, df = dim(ranFeatsDF[ranFeatsDF$random == paste0("Random ", k),])[1]-1 ) *
      ranFeatsDF_randomSEM[k] )
})
ranFeatsDF_randomCIupper <- sapply(1:quantiles, function(k) {
  ranFeatsDF_randomMean[k] +
    ( qt(0.975, df = dim(ranFeatsDF[ranFeatsDF$random == paste0("Random ", k),])[1]-1 ) *
      ranFeatsDF_randomSEM[k] )
})
ranFeatsDF_summary_stats <- data.frame(random = paste0("Random ", 1:quantiles),
                                       Mean = ranFeatsDF_randomMean,
                                       SD = ranFeatsDF_randomSD,
                                       SEM = ranFeatsDF_randomSEM,
                                       CIlower = ranFeatsDF_randomCIlower,
                                       CIupper = ranFeatsDF_randomCIupper,
                                       stringsAsFactors = F)
summary_stats_min <- min(c(featuresDF_summary_stats$CIlower, ranFeatsDF_summary_stats$CIlower), na.rm = T)
summary_stats_max <- max(c(featuresDF_summary_stats$CIupper, ranFeatsDF_summary_stats$CIupper), na.rm = T)

featuresDF <- featuresDF[which(featuresDF[,which(colnames(featuresDF) == parName)] <=
                               quantile(featuresDF[,which(colnames(featuresDF) == parName)],
                                        probs = densityProp, na.rm = T)),]
ranFeatsDF <- ranFeatsDF[which(ranFeatsDF[,which(colnames(ranFeatsDF) == parName)] <=
                               quantile(ranFeatsDF[,which(colnames(ranFeatsDF) == parName)],
                                        probs = densityProp, na.rm = T)),]
xmin <- min(c(featuresDF[,which(colnames(featuresDF) == parName)]),
              na.rm = T)
xmax <- max(c(featuresDF[,which(colnames(featuresDF) == parName)]),
              na.rm = T)
minDensity <- 0
maxDensity <- max(density(featuresDF[featuresDF$quantile == "Quantile 4",][,which(colnames(featuresDF) == parName)],
                          na.rm = T)$y)+maxDensityPlus
maxDensity <- max(
  c(
    sapply(1:quantiles, function(k) {
      max(c(max(density(featuresDF[featuresDF$quantile == paste0("Quantile ", k),][,which(colnames(featuresDF) == parName)],
                        na.rm = T)$y),
            max(density(ranFeatsDF[ranFeatsDF$random == paste0("Random ", k),][,which(colnames(featuresDF) == parName)],
                        na.rm = T)$y)))
     })
   )
)+maxDensityPlus

# Define legend labels
legendLabs_feature <- lapply(1:quantiles, function(x) {
  grobTree(textGrob(bquote(.(paste0("Quantile ", 1:quantiles)[x])),
                    x = legendLabX, y = 0.90-((x-1)*0.07), just = "left",
                    gp = gpar(col = chrColours[x], fontsize = 22)))
})
legendLabs_ranFeat <- lapply(1:quantiles, function(x) {
  grobTree(textGrob(bquote(.(paste0("Random ", 1:quantiles)[x])),
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
  scale_x_continuous(limits = c(xmin, xmax),
                     labels = function(x) sprintf(xDec, x)) +
  scale_y_continuous(limits = c(minDensity, maxDensity),
                     labels = function(x) sprintf(yDec, x)) +
  labs(x = parameterLab,
       y = "Density") +
  annotation_custom(legendLabs[[1]]) +
  annotation_custom(legendLabs[[2]]) +
  annotation_custom(legendLabs[[3]]) +
  annotation_custom(legendLabs[[4]]) +
  theme_bw() +
  theme(axis.line.y = element_line(size = 2.0, colour = "black"),
        axis.line.x = element_line(size = 2.0, colour = "black"),
        axis.ticks.y = element_line(size = 2.0, colour = "black"),
        axis.ticks.x = element_line(size = 2.0, colour = "black"),
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
  scale_y_continuous(limits = c(summary_stats_min, summary_stats_max),
                     labels = function(x) sprintf(yDec2, x)) +
  labs(x = "",
       y = parameterLab) +
  theme_bw() +
  theme(axis.line.y = element_line(size = 2.0, colour = "black"),
        axis.ticks.y = element_line(size = 2.0, colour = "black"),
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

ggObjGA_feature <- popgen_stats_plotFun(lociDF = featuresDF[grepl("Quantile ", featuresDF$quantile),],
                                        parameter = parName,
                                        parameterLab = bquote(.(parNamePlot)),
                                        featureGroup = "quantile", 
                                        featureNamePlot = featureNamePlot,
                                        legendLabs = legendLabs_feature,
                                        chrColours = chrColours
                                       )
ggObjGA_ranFeat <- popgen_stats_plotFun(lociDF = ranFeatsDF[grepl("Random ", ranFeatsDF$random),],
                                        parameter = parName,
                                        parameterLab = bquote(.(parNamePlot)),
                                        featureGroup = "random", 
                                        featureNamePlot = ranFeatNamePlot,
                                        legendLabs = legendLabs_ranFeat,
                                        chrColours = chrColours
                                       )
ggObjGA_feature_mean <- popgen_stats_meanCIs(dataFrame = featuresDF_summary_stats,
                                             parameterLab = bquote(.(parNamePlot)),
                                             featureGroup = "quantile",
                                             featureNamePlot = featureNamePlot,
                                             chrColours = chrColours
                                            )
ggObjGA_ranFeat_mean <- popgen_stats_meanCIs(dataFrame = ranFeatsDF_summary_stats,
                                             parameterLab = bquote(.(parNamePlot)),
                                             featureGroup = "random",
                                             featureNamePlot = ranFeatNamePlot,
                                             chrColours = chrColours
                                            )
ggObjGA_combined <- grid.arrange(ggObjGA_feature,
                                 ggObjGA_feature_mean,
                                 ggObjGA_ranFeat,
                                 ggObjGA_ranFeat_mean,
                                 ncol = 2, as.table = F)
ggsave(paste0(plotDir,
              parName, "_densityProp", densityProp, "_around_", quantiles, "quantiles",
              "_by_", orderingFactor,
              "_of_CEN180_in_T2T_Col_",
              paste0(chrName, collapse = "_"), ".pdf"),
       plot = ggObjGA_combined,
       height = 13, width = 14)
