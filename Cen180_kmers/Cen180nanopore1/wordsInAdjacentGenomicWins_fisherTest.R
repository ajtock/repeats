#!/applications/R/R-3.5.0/bin/Rscript

# Perform a Fisher's exact test on the counts of each word of wordSize nt
# in each genomic window and shuffled window of winSize nt

# Usage:
# wordsInAdjacentGenomicWins_fisherTest.R 10 1000 1kb 0.1

wordSize <- 10
wordName <- as.character(wordSize)
winSize <- 1000
winName <- "1kb"
FDR <- 0.1
FDRname <- as.character(FDR)

args <- commandArgs(trailingOnly = T)
wordSize <- as.numeric(args[1])
wordName <- as.character(wordSize)
winSize <- as.numeric(args[2])
winName <- as.character(args[3])
FDR <- as.numeric(args[4])
FDRname <- as.character(FDR)

## REMOVE
#setwd("/home/ajt200/analysis/repeats/Cen180_kmers/Cen180nanopore1/")
dir0 <- "./plots/"
dir1 <- paste0(dir0, "wordSize", wordName)
dir2 <- paste0(dir1, "/winSize", winName)
plotDir <- paste0(dir2, "/FDR", FDRname, "/")
system(paste0("[ -d ", dir0, " ] || mkdir ", dir0))
system(paste0("[ -d ", dir1, " ] || mkdir ", dir1))
system(paste0("[ -d ", dir2, " ] || mkdir ", dir2))
system(paste0("[ -d ", plotDir, " ] || mkdir ", plotDir))

library(parallel)
library(doParallel)
registerDoParallel(cores = 5)
print("Currently registered parallel backend name, version and cores")
print(getDoParName())
print(getDoParVersion())
print(getDoParWorkers())
library(BSgenome.Athaliana.Cambridge.Cen180nanopore1)
library(BSgenome.Athaliana.TAIR.TAIR9)
ls("package:BSgenome.Athaliana.Cambridge.Cen180nanopore1")
ls("package:BSgenome.Athaliana.TAIR.TAIR9")
Cen180 <- BSgenome.Athaliana.Cambridge.Cen180nanopore1
genome <- BSgenome.Athaliana.TAIR.TAIR9
print(Cen180)
print(genome)

# Load windowed word counts
# (objects named winCounts and ranCounts)
load(paste0("./wordsInAdjacentCen180nanopore1Wins",
            "_wordSize", wordName,
            "_winSize", winName,
            ".RData"))
load(paste0("./wordsInAdjacentCen180nanopore1ShuffledWins",
            "_wordSize", wordName,
            "_winSize", winName,
            ".RData"))
# Load windowed word counts (winCounts and ranCounts) for
# reverse-complemented sequences in adjacent windows
# (objects named winCountsRC and ranCountsRC)
load(paste0("./wordsInAdjacentCen180nanopore1WinsRC",
            "_wordSize", wordName,
            "_winSize", winName,
            ".RData"))
load(paste0("./wordsInAdjacentCen180nanopore1ShuffledWinsRC",
            "_wordSize", wordName,
            "_winSize", winName,
            ".RData"))

# Function to create a 2x2 contingency table of counts for a word of wordSize nt
# inside and outside a genomic window and a shuffled window of winSize nt
contingencyTable <- function(inWin, outWin, inRan, outRan) {
  conTab <- matrix(c(inWin, outWin, inRan, outRan), ncol = 2)
  rownames(conTab) <- c("in", "out")
  colnames(conTab) <- c("win", "ran")
  conTab
}

# Perform Fisher's exact test on each contingency table
fisherPvalsList <- mclapply(seq_along(winCounts), function(x) {
    sapply(seq_along(winCounts[[x]]), function(y) {
      fisher.test(contingencyTable(inWin = winCounts[[x]][y],
                                   outWin = sum(winCounts[[x]])-winCounts[[x]][y],
                                   inRan = ranCounts[[x]][y],
                                   outRan = sum(ranCounts[[x]])-ranCounts[[x]][y]),
                  alternative = "greater")$p.value
    })
}, mc.cores = detectCores())
# And for reverse-complemented sequences in adjacent windows
fisherPvalsListRC <- mclapply(seq_along(winCountsRC), function(x) {
    sapply(seq_along(winCountsRC[[x]]), function(y) {
      fisher.test(contingencyTable(inWin = winCountsRC[[x]][y],
                                   outWin = sum(winCountsRC[[x]])-winCountsRC[[x]][y],
                                   inRan = ranCountsRC[[x]][y],
                                   outRan = sum(ranCountsRC[[x]])-ranCountsRC[[x]][y]),
                  alternative = "greater")$p.value
    })
}, mc.cores = detectCores())

# Assign word names to elements in list of P-values
names(fisherPvalsList) <- names(winCounts)
names(fisherPvalsListRC) <- names(winCountsRC)

# Adjust P-values by applying the Benjamini-Hochberg multiple testing correction
fisherAdjPvalsList <- mclapply(seq_along(fisherPvalsList), function(x) {
  p.adjust(p = fisherPvalsList[[x]], method = "BH")
}, mc.cores = detectCores())
# And for reverse-complemented sequences in adjacent windows
fisherAdjPvalsListRC <- mclapply(seq_along(fisherPvalsListRC), function(x) {
  p.adjust(p = fisherPvalsListRC[[x]], method = "BH")
}, mc.cores = detectCores())
# Assign word names to elements in list of adjusted P-values
names(fisherAdjPvalsList) <- names(fisherPvalsList)
names(fisherAdjPvalsListRC) <- names(fisherPvalsListRC)

# Determine which words are over-represented in at least one window
overRepWordsBool <- mclapply(seq_along(fisherAdjPvalsList), function(x) {
  min(fisherAdjPvalsList[[x]]) < FDR
}, mc.cores = detectCores())
names(overRepWordsBool) <- names(fisherAdjPvalsList)
overRepWords <- names(overRepWordsBool)[overRepWordsBool == TRUE]

# Confirm that the sum of "TRUE" elements in list overRepWordsBool is equal
# to the number of words in overRepWords
stopifnot(identical(sum(unlist(overRepWordsBool[names(overRepWordsBool) %in% overRepWords])),
                    length(overRepWords)))

# Retain words that are over-represented in at least one window
if(length(overRepWords) > 0) {
  sigFisherAdjPvalsList <- fisherAdjPvalsList[names(fisherAdjPvalsList) %in% overRepWords]
  save(sigFisherAdjPvalsList,
       file = paste0("./wordsInAdjacentCen180nanopore1Wins",
                     "_wordSize", wordName,
                     "_winSize", winName,
                     "_FDR", FDRname,
                     "_sigFisherAdjPvals.RData"))
  sigFisherPvalsList <- fisherPvalsList[names(fisherPvalsList) %in% overRepWords]
  save(sigFisherPvalsList,
       file = paste0("./wordsInAdjacentCen180nanopore1Wins",
                     "_wordSize", wordName,
                     "_winSize", winName,
                     "_FDR", FDRname,
                     "_sigFisherPvalsList.RData"))
  
  sigWinCounts <- winCounts[names(winCounts) %in% overRepWords]
  save(sigWinCounts,
       file = paste0("./wordsInAdjacentCen180nanopore1Wins",
                     "_wordSize", wordName,
                     "_winSize", winName,
                     "_FDR", FDRname,
                     "_sigWinCounts.RData"))
  
  sigRanCounts <- ranCounts[names(ranCounts) %in% overRepWords]
  save(sigRanCounts,
       file = paste0("./wordsInAdjacentCen180nanopore1ShuffledWins",
                     "_wordSize", wordName,
                     "_winSize", winName,
                     "_FDR", FDRname,
                     "_sigRanCounts.RData"))
} else {
  print("There were no over-represented words")
}


## And for reverse-complemented sequences in adjacent windows
# Determine which words are over-represented in at least one window
overRepWordsBoolRC <- mclapply(seq_along(fisherAdjPvalsListRC), function(x) {
  min(fisherAdjPvalsListRC[[x]]) < FDR
}, mc.cores = detectCores())
names(overRepWordsBoolRC) <- names(fisherAdjPvalsListRC)
overRepWordsRC <- names(overRepWordsBoolRC)[overRepWordsBoolRC == TRUE]

# Confirm that the sum of "TRUE" elements in list overRepWordsBoolRC is equal
# to the number of words in overRepWordsRC
stopifnot(identical(sum(unlist(overRepWordsBoolRC[names(overRepWordsBoolRC) %in% overRepWordsRC])),
                    length(overRepWordsRC)))

# Retain words that are over-represented in at least one window
if(length(overRepWordsRC) > 0) {
  sigFisherAdjPvalsListRC <- fisherAdjPvalsListRC[names(fisherAdjPvalsListRC) %in% overRepWordsRC]
  save(sigFisherAdjPvalsListRC,
       file = paste0("./wordsInAdjacentCen180nanopore1WinsRC",
                     "_wordSize", wordName,
                     "_winSize", winName,
                     "_FDR", FDRname,
                     "_sigFisherAdjPvals.RData"))
  sigFisherPvalsListRC <- fisherPvalsListRC[names(fisherPvalsListRC) %in% overRepWordsRC]
  save(sigFisherPvalsListRC,
       file = paste0("./wordsInAdjacentCen180nanopore1WinsRC",
                     "_wordSize", wordName,
                     "_winSize", winName,
                     "_FDR", FDRname,
                     "_sigFisherPvalsListRC.RData"))
  
  sigWinCountsRC <- winCountsRC[names(winCountsRC) %in% overRepWordsRC]
  save(sigWinCountsRC,
       file = paste0("./wordsInAdjacentCen180nanopore1WinsRC",
                     "_wordSize", wordName,
                     "_winSize", winName,
                     "_FDR", FDRname,
                     "_sigWinCountsRC.RData"))
  
  sigRanCountsRC <- ranCountsRC[names(ranCountsRC) %in% overRepWordsRC]
  save(sigRanCountsRC,
       file = paste0("./wordsInAdjacentCen180nanopore1ShuffledWinsRC",
                     "_wordSize", wordName,
                     "_winSize", winName,
                     "_FDR", FDRname,
                     "_sigRanCountsRC.RData"))
} else {
  print("There were no over-represented words")
}

stopifnot(length(overRepWords) + length(overRepWordsRC) > 0)

# Function to plot -log10-transformed P-values derived from Fisher's exact tests
minusLog10PvalPlot <- function(xplot, word,
                               Pvals, PvalsCol, minPval, maxPval,
                               AdjPvals, AdjPvalsCol, minAdjPval, maxAdjPval, RC) {
  plot(x = xplot, y = -log10(Pvals), col = PvalsCol, type = "l", lwd = 1.5,
       ylim = c(minPval,
                maxPval),
       xlab = "", ylab = "",
       xaxt = "n", yaxt = "n",
       main = "")
  mtext(side = 2, line = 2.25, cex = 1, col = PvalsCol,
        text = bquote("-Log"[10]*"("*italic("P")*"-value) "*.(winName)^-1))
  axis(side = 2, cex.axis = 1.5, lwd.tick = 1.5)
  abline(h = -log10(0.05), lty = 5, lwd = 1, col = PvalsCol) 

  par(new = T)
  plot(x = xplot, y = -log10(AdjPvals), col = AdjPvalsCol, type = "l", lwd = 1.5,
       ylim = c(minAdjPval,
                maxAdjPval),
       xlab = "", ylab = "",
       xaxt = "n", yaxt = "n",
       main = paste0("Fisher's exact tests on ", word,
                     " (", wordName, " nt) counts along Cen180nanopore1 sequences ", RC),
       cex.main = 1.0)
  p <- par('usr')
  text(p[2], mean(p[3:4]), cex = 1, adj = c(0.5, -1.0), xpd = NA, srt = -90, col = AdjPvalsCol,
       labels = bquote(atop("-Log"[10]*"(BH-adjusted ", italic("P")*"-value) "*.(winName)^-1)))
  axis(side = 4, cex.axis = 1.5, lwd.tick = 1.5)
  axis(side = 1, lwd.tick = 1.5, cex.axis = 1.5,
       at = c(0, 2e5, 4e5, 6e5, 8e5, 1e6),
       labels = c("0", "0.2", "0.4", "0.6", "0.8", "1.0"))
  abline(h = -log10(FDR), lty = 5, lwd = 1, col = AdjPvalsCol) 

  abline(v = c(1, xplot[length(xplot)]), lty = 1, lwd = 1, col = "black")
  box(lwd = 1.5) 
}

# Function to plot word counts in windows containing genomic and shuffled sequence
countPlot <- function(xplot, word, minCount, maxCount,
                      winCounts, winCountsCol,
                      ranCounts, ranCountsCol, RC) {
  plot(x = xplot, y = winCounts, col = winCountsCol, type = "l", lwd = 1.5,
       ylim = c(minCount, maxCount),
       xlab = "", ylab = "",
       xaxt = "n", yaxt = "n",
       main = paste0(word,
                     " (", wordName, " nt) counts along Cen180nanopore1 sequences ", RC),
       cex.main = 1.0)
  lines(x = xplot, y = ranCounts, col = ranCountsCol, type = "l", lwd = 1.5)
  axis(side = 1, lwd.tick = 1.5, cex.axis = 1.5,
       at = c(0, 2e5, 4e5, 6e5, 8e5, 1e6),
       labels = c("0", "0.2", "0.4", "0.6", "0.8", "1.0"))
  axis(side = 2, cex.axis = 1.5, lwd.tick = 1.5)
  mtext(side = 2, line = 2.25, cex = 1, col = "black",
        text = bquote(.(word)*" counts "*.(winName)^-1))
  legend("topright",
         legend = c(paste0("Cen180nanopore1 ", RC), paste0("Shuffled ", RC)),
         col = c(winCountsCol, ranCountsCol),
         text.col = c(winCountsCol, ranCountsCol),
         text.font = c(1, 1),
         ncol = 1, cex = 0.7, lwd = 1.5, bty = "n")

  abline(v = c(1, xplot[length(xplot)]), lty = 1, lwd = 1, col = "black")
  box(lwd = 1.5) 
}

# Plot chromosome profiles of P-values and counts for
# words that are over-represented in at least one genomic window
mclapply(seq_along(sigFisherAdjPvalsList), function(x) {
  pdf(paste0(plotDir,
             names(sigFisherAdjPvalsList)[x],
             "_inAdjacentWins",
             "_wordSize", wordName,
             "_winSize", winName,
             "_FDR", FDRname, ".pdf"),
      height = 2.5, width = 18)
  par(mfrow = c(1, 2))
  par(mar = c(2.1, 4.1, 2.1, 4.1))
  par(mgp = c(3, 1, 0))

  minSigPval <- min(-log10(pmax(sigFisherPvalsList[[x]], 1e-323)))
  maxSigPval <- max(-log10(pmax(sigFisherPvalsList[[x]], 1e-323)))
  minSigAdjPval <- min(-log10(pmax(sigFisherAdjPvalsList[[x]], 1e-323)))
  maxSigAdjPval <- max(-log10(pmax(sigFisherAdjPvalsList[[x]], 1e-323)))
  minCount <- min(c(sigWinCounts[[x]], sigRanCounts[[x]]))
  maxCount <- max(c(sigWinCounts[[x]], sigRanCounts[[x]]))

  # Define adjacent windows of width winSize nt
  winStarts <- seq(from = 1,
                   to = length(Cen180[[1]]),
                   by = winSize)
  minusLog10PvalPlot(xplot = winStarts,
                     word = names(sigFisherAdjPvalsList)[x],
                     Pvals = pmax(sigFisherPvalsList[[x]], 1e-323), PvalsCol = "dodgerblue3",
                     minPval = minSigPval, maxPval = maxSigPval,
                     AdjPvals = pmax(sigFisherAdjPvalsList[[x]], 1e-323), AdjPvalsCol = "red",
                     minAdjPval = minSigAdjPval, maxAdjPval = maxSigAdjPval,
                     RC = "")
  countPlot(xplot = winStarts,
            minCount = minCount, maxCount = maxCount,
            word = names(sigWinCounts)[x],
            winCounts = sigWinCounts[[x]], winCountsCol = "forestgreen",
            ranCounts = sigRanCounts[[x]], ranCountsCol = "grey50",
            RC = "")
  dev.off()
}, mc.cores = detectCores())

# And for reverse-complemented sequences in adjacent windows
# Plot chromosome profiles of P-values and counts for
# words that are over-represented in at least one genomic window
mclapply(seq_along(sigFisherAdjPvalsListRC), function(x) {
  pdf(paste0(plotDir,
             names(sigFisherAdjPvalsListRC)[x],
             "_inAdjacentWinsRC",
             "_wordSize", wordName,
             "_winSize", winName,
             "_FDR", FDRname, ".pdf"),
      height = 3, width = 18)
  par(mfrow = c(1, 2))
  par(mar = c(2.1, 4.1, 2.1, 4.1))
  par(mgp = c(3, 1, 0))

  minSigPval <- min(-log10(pmax(sigFisherPvalsListRC[[x]], 1e-323)))
  maxSigPval <- max(-log10(pmax(sigFisherPvalsListRC[[x]], 1e-323)))
  minSigAdjPval <- min(-log10(pmax(sigFisherAdjPvalsListRC[[x]], 1e-323)))
  maxSigAdjPval <- max(-log10(pmax(sigFisherAdjPvalsListRC[[x]], 1e-323)))
  minCount <- min(c(sigWinCountsRC[[x]], sigRanCountsRC[[x]]))
  maxCount <- max(c(sigWinCountsRC[[x]], sigRanCountsRC[[x]]))

  # Define adjacent windows of width winSize nt
  winStarts <- seq(from = 1,
                   to = length(Cen180[[1]]),
                   by = winSize)
  minusLog10PvalPlot(xplot = winStarts,
                     word = names(sigFisherAdjPvalsListRC)[x],
                     Pvals = pmax(sigFisherPvalsListRC[[x]], 1e-323), PvalsCol = "dodgerblue3",
                     minPval = minSigPval, maxPval = maxSigPval,
                     AdjPvals = pmax(sigFisherAdjPvalsListRC[[x]], 1e-323), AdjPvalsCol = "red",
                     minAdjPval = minSigAdjPval, maxAdjPval = maxSigAdjPval,
                     RC = "(RC)")
  countPlot(xplot = winStarts,
            minCount = minCount, maxCount = maxCount,
            word = names(sigWinCountsRC)[x],
            winCounts = sigWinCountsRC[[x]], winCountsCol = "forestgreen",
            ranCounts = sigRanCountsRC[[x]], ranCountsCol = "grey50",
            RC = "(RC)")
  dev.off()
}, mc.cores = detectCores())
