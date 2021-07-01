#!/applications/R/R-4.0.0/bin/Rscript

# author: Andy Tock
# contact: ajt200@cam.ac.uk
# date: 17.06.2021

# Calculate windowed CEN180 frequencies for CEN180 sequences within orderingFactor quantiles

# Usage:
# /applications/R/R-4.0.0/bin/Rscript CEN180_quantile_frequency_chrProfiles.R 'Chr1,Chr2,Chr3,Chr4,Chr5' wSNV 4 10000 101

#chrName <- unlist(strsplit("Chr1,Chr2,Chr3,Chr4,Chr5",
#                           split = ","))
#orderingFactor <- "wSNV"
#quantiles <- 4
#genomeBinSize <- 10000
#maPeriod <- 101

args <- commandArgs(trailingOnly = T)
chrName <- unlist(strsplit(args[1],
                           split = ","))
orderingFactor <- args[2]
quantiles <- as.integer(args[3])
genomeBinSize <- as.integer(args[4])
maPeriod <- as.integer(args[5])

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

# Load table of features grouped into quantiles
featuresDF <- read.table(paste0(outDir,
                                "features_", quantiles, "quantiles",
                                "_by_", orderingFactor,
                                "_of_CEN180_in_t2t-col.20210610_",
                                paste0(chrName, collapse = "_"), ".tsv"),
                         header = T, sep = "\t", row.names = NULL)

# Define list of data.frames in which each element corresponds to a CEN180 quantile
quantilesDFlist <- mclapply(1:quantiles, function(k) {
  featuresDF[featuresDF$quantile == paste0("Quantile ", k),]
}, mc.cores = quantiles)

# Genomic definitions
fai <- read.table("/home/ajt200/analysis/nanopore/t2t-col.20210610/t2t-col.20210610.fa.fai", header = F)
chrs <- fai$V1[1:5]
chrLens <- fai$V2[1:5]
chrLens <- chrLens[chrs %in% chrName]
chrs <- chrs[chrs %in% chrName]

# Define windows as GRanges object
windowsGR <- GRanges()
for(i in 1:length(chrs)) {
  seqWindows <- seq(1, chrLens[i], by = genomeBinSize)
  windowsIR <- IRanges(start = seqWindows,
                       width = genomeBinSize)
  windowsIR <- windowsIR[-length(windowsIR)]
  windowsIR <- append(windowsIR,
                      IRanges(start = seqWindows[length(seqWindows)],
                              end = chrLens[i]))
  chrWindowsGR <- GRanges(seqnames = chrs[i],
                          ranges = windowsIR,
                          strand = "*")
  print(chrWindowsGR)
  windowsGR <- append(windowsGR, chrWindowsGR)
}

quantileProfiles <- NULL
for(k in 1:quantiles) {
  quantilekProfile <- NULL
  for(i in 1:length(chrs)) {
    # Count quantile features within adjacent windows
    featuresChr <- quantilesDFlist[[k]][quantilesDFlist[[k]]$chr == chrs[i],]
    windowsChrGR <- windowsGR[seqnames(windowsGR) == chrs[i]]
    if(dim(featuresChr)[1] > 0) {
      featuresChrGR <- GRanges(seqnames = chrs[i],
                               ranges = IRanges(start = featuresChr$start,
                                                end = featuresChr$end),
                               strand = featuresChr$strand)
    } else {
      featuresChrGR <- GRanges()
    }
    winfeatures <- countOverlaps(query = windowsChrGR,
                                 subject = featuresChrGR,
                                 type = "any",
                                 ignore.strand = T)
    profileChr <- data.frame(chr = as.character(chrs[i]),
                             window = as.integer(start(windowsChrGR)),
                             features = as.integer(winfeatures),
                             stringsAsFactors = F)
    quantilekProfile <- rbind(quantilekProfile, profileChr)
  }
  quantileProfiles <- cbind(quantileProfiles, quantilekProfile[,3])
}
quantileProfiles <- data.frame(quantilekProfile[,1:2],
                               quantileProfiles,
                               stringsAsFactors = F)
colnames(quantileProfiles) <- c("chr", "window",
                                paste0("quantile", 1:quantiles))
write.table(quantileProfiles,
            file = paste0(outDir,
                          "CEN180_frequency_per_", genomeBinName,
                          "_", quantiles, "quantiles_",
                          "_by_", orderingFactor,
                          "_of_CEN180_in_t2t-col.20210610_",
                          paste0(chrName, collapse = "_"), "_unsmoothed.tsv"),
            row.names = F, col.names = T, quote = F, sep = "\t")

# Calculate moving average of current window,
#### (maPeriod/2) previous windows (where maPeriod is even) OR
# (maPeriod/2)-0.5 previous windows (where maPeriod is odd),
# and
#### (maPeriod/2) subsequent windows (where maPeriod is even) OR
# (maPeriod/2)-0.5 subsequent windows (where maPeriod is odd)
# (the higher maPeriod is, the greater the smoothing)
stopifnot(maPeriod %% 2 != 0)
flank <- (maPeriod/2)-0.5
# Define MA filter coefficients
f <- rep(1/maPeriod, maPeriod)

chrProfiles <- mclapply(seq_along(chrs), function(x) {
  quantileProfiles[quantileProfiles$chr == chrs[x],]
}, mc.cores = length(chrs))

filt_chrProfiles <- mclapply(seq_along(chrProfiles), function(x) {
  filt_quantileProfiles <- NULL
  for(k in 1:quantiles) {
    filt_chrProfilek <- stats::filter(x = chrProfiles[[x]][,k+2],
                                      filter = f,
                                      sides = 2)
    filt_chrProfilek[1:flank] <- filt_chrProfilek[flank+1]
    filt_chrProfilek[(length(filt_chrProfilek)-flank+1):length(filt_chrProfilek)] <- filt_chrProfilek[(length(filt_chrProfilek)-flank)]
    filt_quantilekProfile <- data.frame(chr = as.character(chrProfiles[[x]]$chr),
                                        window = as.integer(chrProfiles[[x]]$window),
                                        filt_features = as.numeric(filt_chrProfilek),
                                        stringsAsFactors = F)
    filt_quantileProfiles <- cbind(filt_quantileProfiles, filt_quantilekProfile[,3])
  }
  filt_quantileProfiles <- data.frame(filt_quantilekProfile[,1:2],
                                      filt_quantileProfiles,
                                      stringsAsFactors = F)
  colnames(filt_quantileProfiles) <- c("chr", "window",
                                       paste0("filt_quantile", 1:quantiles))
  return(filt_quantileProfiles)
}, mc.cores = length(chrProfiles))

# Combine list of 1 data.frame per chromosome into one data.frame
filt_quantileProfiles <- do.call(rbind, filt_chrProfiles)
write.table(filt_quantileProfiles,
            file = paste0(outDir,
                          "CEN180_frequency_per_", genomeBinName,
                          "_", quantiles, "quantiles_",
                          "_by_", orderingFactor,
                          "_of_CEN180_in_t2t-col.20210610_",
                          paste0(chrName, collapse = "_"), "_smoothed.tsv"),
            row.names = F, col.names = T, quote = F, sep = "\t")
