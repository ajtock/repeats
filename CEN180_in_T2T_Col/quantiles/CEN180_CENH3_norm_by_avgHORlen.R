#!/applications/R/R-3.5.0/bin/Rscript

# Plot density and means with 95% confidence intervals (CIs)
# of given statistic for each group of CEN180 sequences

# Usage:
# /applications/R/R-3.5.0/bin/Rscript CEN180_CENH3_norm_by_avgHORlen.R wSNV SNVs 'Chr1,Chr2,Chr3,Chr4,Chr5' wSNV 4 1.00 0.2 '%2.0f' '%3.1f' '%2.0f' 0.65
# /applications/R/R-3.5.0/bin/Rscript CEN180_CENH3_norm_by_avgHORlen.R map_K150_E4_in_bodies 'Mappability (k=150 e=4)' 'Chr1,Chr2,Chr3,Chr4,Chr5' wSNV 4 1.00 0.2 '%1.2f' '%3.1f' '%1.2f' 0.65
# /applications/R/R-3.5.0/bin/Rscript CEN180_CENH3_norm_by_avgHORlen.R CENH3_in_bodies 'CENH3' 'Chr1,Chr2,Chr3,Chr4,Chr5' wSNV 4 1.00 0.2 '%1.1f' '%3.1f' '%1.1f' 0.65
# /applications/R/R-3.5.0/bin/Rscript CEN180_CENH3_norm_by_avgHORlen.R H3K9me2_in_bodies 'H3K9me2' 'Chr1,Chr2,Chr3,Chr4,Chr5' wSNV 4 1.00 0.2 '%1.1f' '%3.1f' '%1.1f' 0.65
# /applications/R/R-3.5.0/bin/Rscript CEN180_CENH3_norm_by_avgHORlen.R H3K27me1_in_bodies 'H3K27me1' 'Chr1,Chr2,Chr3,Chr4,Chr5' wSNV 4 1.00 0.2 '%1.1f' '%3.1f' '%1.2f' 0.65
# /applications/R/R-3.5.0/bin/Rscript CEN180_CENH3_norm_by_avgHORlen.R HORlengthsSum 'Activity' 'Chr1,Chr2,Chr3,Chr4,Chr5' wSNV 4 1.00 0.2 '%1.0f' '%1.2f' '%1.0f' 0.65
# /applications/R/R-3.5.0/bin/Rscript CEN180_CENH3_norm_by_avgHORlen.R array_size 'Array size' 'Chr1,Chr2,Chr3,Chr4,Chr5' wSNV 4 1.00 0.2 '%1.0f' '%1.2f' '%1.0f' 0.65

#parName <- "wSNV"
#parNamePlot <- "SNVs"
#chrName <- unlist(strsplit("Chr1,Chr2,Chr3,Chr4,Chr5",
#                           split = ","))
#orderingFactor <- "wSNV"
#quantiles <- 4
#densityProp <- 1.0
#maxDensityPlus <- 0.2
#xDec <- "%1.1f"
#yDec <- "%1.1f"
#yDec2 <- "%1.2f"
#legendLabX <- 0.65

args <- commandArgs(trailingOnly = T)
parName <- args[1]
parNamePlot <- args[2]
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

library(parallel)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggbeeswarm)
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
  ranFeatNamePlot <- paste0("Genome-wide", " ",
                            "random quantiles")
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
  ranFeatNamePlot <- paste0(paste0(chrName, collapse = "_"), " ",
                            "random quantiles")
}

# Define quantile colours
quantileColours <- c("red", "purple", "blue", "navy")
makeTransparent <- function(thisColour, alpha = 250)
{
  newColour <- col2rgb(thisColour)
  apply(newColour, 2, function(x) {
    rgb(red = x[1], green = x[2], blue = x[3],
        alpha = alpha, maxColorValue = 255)
  })
}
quantileColours <- makeTransparent(quantileColours)

# Load table of features grouped into quantiles
featuresDF <- read.table(paste0(outDir,
                                "features_", quantiles, "quantiles",
                                "_by_", orderingFactor,
                                "_of_CEN180_in_T2T_Col_",
                                paste0(chrName, collapse = "_"), ".tsv"),
                         header = T, sep = "\t", row.names = NULL)

# Load features to confirm feature (row) ordering in "featuresDF" is the same
# as in "features" (which was used for generating the coverage matrices)
features <- lapply(seq_along(chrName), function(y) {
  read.table(paste0("CEN180_in_T2T_Col_",
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

featuresDFtmp <- data.frame(featuresDF,
                            random = as.character(""),
                            stringsAsFactors = F)
ranFeatsDF <- data.frame()
for(k in 1:quantiles) {
  featuresDFtmp[randomPCIndices[[k]],]$random <- paste0("Random ", k)
  ranFeatsDFk <- featuresDFtmp[featuresDFtmp$random == paste0("Random ", k),]
  ranFeatsDF <- rbind(ranFeatsDF, ranFeatsDFk)
}

featuresDF_avgHORlen <- (featuresDF$HORlengthsSum+1) / (featuresDF$HORcount+1)
featuresDF[,colnames(featuresDF) == parName,] <- featuresDF[,colnames(featuresDF) == parName,] /
                                                 featuresDF_avgHORlen 

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
                    gp = gpar(col = quantileColours[x], fontsize = 22)))
})
legendLabs_ranFeat <- lapply(1:quantiles, function(x) {
  grobTree(textGrob(bquote(.(paste0("Random ", 1:quantiles)[x])),
                    x = legendLabX, y = 0.90-((x-1)*0.07), just = "left",
                    gp = gpar(col = quantileColours[x], fontsize = 22)))
})

# Parameter density plot function
popgen_stats_plotFun <- function(lociDF,
                                 parameter,
                                 parameterLab,
                                 featureGroup,
                                 featureNamePlot,
                                 legendLabs,
                                 quantileColours) {
  ggplot(data = lociDF,
         mapping = aes(x = get(parameter),
                       colour = reorder(x = get(featureGroup), X = desc(get(featureGroup))),
                       group = reorder(x = get(featureGroup), X = desc(get(featureGroup))))) +
  scale_colour_manual(values = rev(quantileColours)) +
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
                                 quantileColours) {
  ggplot(data = dataFrame,
         mapping = aes(x = get(featureGroup),
                       y = Mean,
                       colour = get(featureGroup))) +
  labs(colour = "") +
  geom_point(shape = 19, size = 6, position = position_dodge(width = 0.2)) +
  geom_errorbar(mapping = aes(ymin = CIlower,
                              ymax = CIupper),
                width = 0.2, size = 2, position = position_dodge(width = 0.2)) +
  scale_colour_manual(values = quantileColours) +
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
        axis.text.x = element_text(size = 22, colour = quantileColours, hjust = 1.0, vjust = 1.0, angle = 45),
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
                                        quantileColours = quantileColours
                                       )
ggObjGA_ranFeat <- popgen_stats_plotFun(lociDF = ranFeatsDF[grepl("Random ", ranFeatsDF$random),],
                                        parameter = parName,
                                        parameterLab = bquote(.(parNamePlot)),
                                        featureGroup = "random", 
                                        featureNamePlot = ranFeatNamePlot,
                                        legendLabs = legendLabs_ranFeat,
                                        quantileColours = quantileColours
                                       )
ggObjGA_feature_mean <- popgen_stats_meanCIs(dataFrame = featuresDF_summary_stats,
                                             parameterLab = bquote(.(parNamePlot)),
                                             featureGroup = "quantile",
                                             featureNamePlot = featureNamePlot,
                                             quantileColours = quantileColours
                                            )
ggObjGA_ranFeat_mean <- popgen_stats_meanCIs(dataFrame = ranFeatsDF_summary_stats,
                                             parameterLab = bquote(.(parNamePlot)),
                                             featureGroup = "random",
                                             featureNamePlot = ranFeatNamePlot,
                                             quantileColours = quantileColours
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
