#!/applications/R/R-3.5.0/bin/Rscript

# Perform a Fisher's exact test on the counts of each word of wordSize nt
# in each genomic window and shuffled window of winSize nt

# Usage:
# wordsInAdjacentGenomicWins_fisherTest.R 5 100000 100kb 0.1

wordSize <- 5
wordName <- as.character(wordSize)
winSize <- 100000
winName <- "100kb"
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
#setwd("/home/ajt200/analysis/repeats/Piotr_kmers/")
dir1 <- paste0("./plots/wordSize", wordName)
dir2 <- paste0(dir1, "/winSize", winName)
plotDir <- paste0(dir2, "/FDR", FDRname, "/")
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
library(BSgenome)
library(BSgenome.Athaliana.TAIR.TAIR9)
ls("package:BSgenome.Athaliana.TAIR.TAIR9")
genome <- BSgenome.Athaliana.TAIR.TAIR9
print(genome)

# Genomic definitions
chrs <- c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5")
chrLens <- c(30427671, 19698289, 23459830, 18585056, 26975502)
centromeres <- c(15086045, 3607929, 13587786, 3956021, 11725024)
# Pericentromeric regions are as defined in Supplemental Table S26
# of Ziolkowski et al. (2017) Genes Dev. 31
pericenStart <- c(11330001, 990001, 10200001, 990001, 8890001)
pericenEnd <- c(18480000, 7540000, 16860000, 6850000, 15650000)

# Load word counts
winCountsList <- mclapply(seq_along(chrs), function(i) {
  load(paste0("./wordsInAdjacentGenomicWins_", chrs[i],
              "_wordSize", wordName,
              "_winSize", winName,
              ".RData"))
  winCounts
}, mc.cores = length(chrs))
ranCountsList <- mclapply(seq_along(chrs), function(i) { 
  load(paste0("./wordsInAdjacentShuffledWins_", chrs[i],
              "_wordSize", wordName,
              "_winSize", winName,
              ".RData"))
  ranCounts
}, mc.cores = length(chrs))

# Function to create a 2x2 contingency table of counts for a word of wordSize nt
# inside and outside a genomic window and a shuffled window of winSize nt
contingencyTable <- function(inWin, outWin, inRan, outRan) {
  conTab <- matrix(c(inWin, outWin, inRan, outRan), ncol = 2)
  rownames(conTab) <- c("in", "out")
  colnames(conTab) <- c("win", "ran")
  conTab
}

# Perform Fisher's exact test on each contingency table
fisherPvalsList <- lapply(seq_along(winCountsList), function(i) {
  mclapply(seq_along(winCountsList[[i]]), function(x) {
    sapply(seq_along(winCountsList[[i]][[x]]), function(y) {
      fisher.test(contingencyTable(inWin = winCountsList[[i]][[x]][y],
                                   outWin = sum(winCountsList[[i]][[x]])-winCountsList[[i]][[x]][y],
                                   inRan = ranCountsList[[i]][[x]][y],
                                   outRan = sum(ranCountsList[[i]][[x]])-ranCountsList[[i]][[x]][y]),
                  alternative = "greater")$p.value
    })
  }, mc.cores = detectCores())
})
# Assign word names to elements in list of list of P-values
for(i in 1:length(chrs)) {
  names(fisherPvalsList[[i]]) <- names(winCountsList[[i]])
}

# Adjust P-values by applying the Benjamini-Hochberg multiple testing correction
fisherAdjPvalsList <- lapply(seq_along(fisherPvalsList), function(i) {
  mclapply(seq_along(fisherPvalsList[[i]]), function(x) {
    p.adjust(p = fisherPvalsList[[i]][[x]], method = "BH")
  }, mc.cores = detectCores())
})
# Assign word names to elements in list of list of adjusted P-values
for(i in 1:length(chrs)) {
  names(fisherAdjPvalsList[[i]]) <- names(fisherPvalsList[[i]])
}

# Determine which words are over-represented in at least one window
overRepWordsBool <- mclapply(seq_along(fisherAdjPvalsList[[1]]), function(x) {
  minPvalvFDR <- NULL
  for(i in 1:length(chrs)) {
    minPvalChrvFDR <- min(fisherAdjPvalsList[[i]][[x]]) < FDR
    minPvalvFDR <- c(minPvalvFDR, minPvalChrvFDR)
  }
  sum(minPvalvFDR) > 0
}, mc.cores = detectCores())
names(overRepWordsBool) <- names(fisherAdjPvalsList[[1]])
overRepWords <- names(overRepWordsBool)[overRepWordsBool == TRUE]
# Confirm that the sum of "TRUE" elements in list overRepWordsBool is equal
# to the number of words in overRepWords
stopifnot(identical(sum(unlist(overRepWordsBool[names(overRepWordsBool) %in% overRepWords])),
                    length(overRepWords)))

# Retain words that are over-represented in at least one window
sigFisherAdjPvalsList <- lapply(seq_along(fisherAdjPvalsList), function(i) {
  fisherAdjPvalsList[[i]][names(fisherAdjPvalsList[[i]]) %in% overRepWords]
})
save(sigFisherAdjPvalsList,
     file = paste0("./wordsInAdjacentGenomicWins",
                   "_wordSize", wordName,
                   "_winSize", winName,
                   "_FDR", FDRname,
                   "_sigFisherAdjPvalsList.RData"))
sigFisherPvalsList <- lapply(seq_along(fisherPvalsList), function(i) {
  fisherPvalsList[[i]][names(fisherPvalsList[[i]]) %in% overRepWords]
})
save(sigFisherPvalsList,
     file = paste0("./wordsInAdjacentGenomicWins",
                   "_wordSize", wordName,
                   "_winSize", winName,
                   "_FDR", FDRname,
                   "_sigFisherPvalsList.RData"))

sigWinCountsList <- lapply(seq_along(winCountsList), function(i) {
  winCountsList[[i]][names(winCountsList[[i]]) %in% overRepWords]
})
save(sigWinCountsList,
     file = paste0("./wordsInAdjacentGenomicWins",
                   "_wordSize", wordName,
                   "_winSize", winName,
                   "_FDR", FDRname,
                   "_sigWinCountsList.RData"))

sigRanCountsList <- lapply(seq_along(ranCountsList), function(i) {
  ranCountsList[[i]][names(ranCountsList[[i]]) %in% overRepWords]
})
save(sigRanCountsList,
     file = paste0("./wordsInAdjacentShuffledWins",
                   "_wordSize", wordName,
                   "_winSize", winName,
                   "_FDR", FDRname,
                   "_sigRanCountsList.RData"))

# Sum up counts for all words that are over-represented in at least one genomic window
sigSumWinCounts <- mclapply(seq_along(sigWinCountsList), function(i) {
  sigWinCountsChr <- 0
  for(x in 1:length(sigWinCountsList[[i]])) {
    sigWinCountsChr <- sigWinCountsChr+sigWinCountsList[[i]][[x]]
  }
  sigWinCountsChr
}, mc.cores = length(chrs))

sigSumRanCounts <- mclapply(seq_along(sigRanCountsList), function(i) {
  sigRanCountsChr <- 0
  for(x in 1:length(sigRanCountsList[[i]])) {
    sigRanCountsChr <- sigRanCountsChr+sigRanCountsList[[i]][[x]]
  }
  sigRanCountsChr
}, mc.cores = length(chrs))

# Perform Fisher's exact test on contingency table of summed counts
fisherPvalsSumCounts <- mclapply(seq_along(sigSumWinCounts), function(i) {
  sapply(seq_along(sigSumWinCounts[[i]]), function(y) {
    fisher.test(contingencyTable(inWin = sigSumWinCounts[[i]][y],
                                 outWin = sum(sigSumWinCounts[[i]])-sigSumWinCounts[[i]][y],
                                 inRan = sigSumRanCounts[[i]][y],
                                 outRan = sum(sigSumRanCounts[[i]])-sigSumRanCounts[[i]][y]),
                alternative = "greater")$p.value
  })
}, mc.cores = length(chrs))
# Adjust P-values by applying the Benjamini-Hochberg multiple testing correction
fisherAdjPvalsSumCounts <- mclapply(seq_along(fisherPvalsSumCounts), function(i) {
  p.adjust(p = fisherPvalsSumCounts[[i]], method = "BH")
}, mc.cores = length(chrs))

# Function to plot -log10-transformed P-values derived from Fisher's exact tests
minusLog10PvalPlot <- function(xplot, word, chrom,
                               Pvals, PvalsCol, minPval, maxPval,
                               AdjPvals, AdjPvalsCol, minAdjPval, maxAdjPval) {
  plot(x = xplot, y = -log10(Pvals), col = PvalsCol, type = "l", lwd = 1.5,
       ylim = c(minPval,
                maxPval),
       xlim = c(0, max(chrLens)),
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
       xlim = c(0, max(chrLens)),
       xlab = "", ylab = "",
       xaxt = "n", yaxt = "n",
       main = paste0("Fisher's exact tests on ", word,
                     " (", wordName, " nt) counts in windows along chromosome ", chrom),
       cex.main = 1.5)
  p <- par('usr')
  text(p[2], mean(p[3:4]), cex = 1.5, adj = c(0.5, -1.8), xpd = NA, srt = -90, col = AdjPvalsCol,
       labels = bquote("-Log"[10]*"(BH-adjusted "*italic("P")*"-value) "*.(winName)^-1))
  axis(side = 4, cex.axis = 1.5, lwd.tick = 1.5)
  axis(side = 1, cex.axis = 1.5, lwd.tick = 1.5)
  abline(h = -log10(FDR), lty = 5, lwd = 1, col = AdjPvalsCol) 

  abline(v = c(1, xplot[length(xplot)]), lty = 1, lwd = 1, col = "black")
  abline(v = centromeres[chrom], lty = 5, lwd = 1, col = "black")
  abline(v = c(pericenStart[chrom], pericenEnd[chrom]), lty = 5, lwd = 1, col = "grey75")
  box(lwd = 1.5) 
}

# Function to plot word counts in windows containing genomic and shuffled sequence
countPlot <- function(xplot, word, chrom, minCount, maxCount,
                      winCounts, winCountsCol,
                      ranCounts, ranCountsCol) {
  plot(x = xplot, y = winCounts, col = winCountsCol, type = "l", lwd = 1.5,
       ylim = c(minCount, maxCount),
       xlim = c(0, max(chrLens)),
       xlab = "", ylab = "",
       xaxt = "n", yaxt = "n",
       main = paste0(word,
                     " (", wordName, " nt) counts in windows along chromosome ", chrom),
       cex.main = 1.5)
  lines(x = xplot, y = ranCounts, col = ranCountsCol, type = "l", lwd = 1.5)
  axis(side = 1, cex.axis = 1.5, lwd.tick = 1.5)
  axis(side = 2, cex.axis = 1.5, lwd.tick = 1.5)
  mtext(side = 2, line = 2.25, cex = 1, col = "black",
        text = bquote(.(word)*" counts "*.(winName)^-1))
  legend("topright",
         legend = c("Genomic", "Shuffled"),
         col = c(winCountsCol, ranCountsCol),
         text.col = c(winCountsCol, ranCountsCol),
         text.font = c(1, 1),
         ncol = 1, cex = 1.5, lwd = 1.5, bty = "n")

  abline(v = c(1, xplot[length(xplot)]), lty = 1, lwd = 1, col = "black")
  abline(v = centromeres[chrom], lty = 5, lwd = 1, col = "black")
  abline(v = c(pericenStart[chrom], pericenEnd[chrom]), lty = 5, lwd = 1, col = "grey75")
  box(lwd = 1.5) 
}

# Plot chromosome profiles of P-values and counts for
# words that are over-represented in at least one genomic window
mclapply(seq_along(sigFisherAdjPvalsList[[1]]), function(x) {
  pdf(paste0(plotDir,
             names(sigFisherAdjPvalsList[[1]])[x],
             "_inAdjacentWins",
             "_wordSize", wordName,
             "_winSize", winName,
             "_FDR", FDRname, ".pdf"),
      height = 15, width = 18)
  par(mfrow = c(5, 2))
  par(mar = c(2.1, 4.1, 2.1, 4.1))
  par(mgp = c(3, 1, 0))

  minSigPval <- NULL
  maxSigPval <- NULL
  minSigAdjPval <- NULL
  maxSigAdjPval <- NULL
  minCount <- NULL
  maxCount <- NULL
  for(i in 1:length(chrs)) {
    minSigPvalChr <- min(-log10(pmax(sigFisherPvalsList[[i]][[x]], 1e-323)))
    maxSigPvalChr <- max(-log10(pmax(sigFisherPvalsList[[i]][[x]], 1e-323)))
    minSigAdjPvalChr <- min(-log10(pmax(sigFisherAdjPvalsList[[i]][[x]], 1e-323)))
    maxSigAdjPvalChr <- max(-log10(pmax(sigFisherAdjPvalsList[[i]][[x]], 1e-323)))
    minCountChr <- min(c(sigWinCountsList[[i]][[x]], sigRanCountsList[[i]][[x]]))
    maxCountChr <- max(c(sigWinCountsList[[i]][[x]], sigRanCountsList[[i]][[x]]))
 
    minSigPval <- min(c(minSigPval, minSigPvalChr)) 
    maxSigPval <- max(c(maxSigPval, maxSigPvalChr))
    minSigAdjPval <- min(c(minSigAdjPval, minSigAdjPvalChr))
    maxSigAdjPval <- max(c(maxSigAdjPval, maxSigAdjPvalChr))
    minCount <- min(c(minCount, minCountChr))
    maxCount <- max(c(maxCount, maxCountChr))
  }

  for(i in 1:length(chrs)) {
    # Define adjacent windows of width winSize nt
    winStarts <- seq(from = 1,
                     to = length(genome[[i]]),
                     by = winSize)
    minusLog10PvalPlot(xplot = winStarts,
                       chrom = i,
                       word = names(sigFisherAdjPvalsList[[i]])[x],
                       Pvals = pmax(sigFisherPvalsList[[i]][[x]], 1e-323), PvalsCol = "dodgerblue3",
                       minPval = minSigPval, maxPval = maxSigPval,
                       AdjPvals = pmax(sigFisherAdjPvalsList[[i]][[x]], 1e-323), AdjPvalsCol = "red",
                       minAdjPval = minSigAdjPval, maxAdjPval = maxSigAdjPval)
    countPlot(xplot = winStarts,
              chrom = i,
              minCount = minCount, maxCount = maxCount,
              word = names(sigWinCountsList[[i]])[x],
              winCounts = sigWinCountsList[[i]][[x]], winCountsCol = "forestgreen",
              ranCounts = sigRanCountsList[[i]][[x]], ranCountsCol = "grey50")
  }
dev.off()
}, mc.cores = detectCores())

# Plot chromosome profiles of P-values corresponding to summed counts for
# words that are over-represented in at least one genomic window
pdf(paste0(plotDir,
           "overRepWordsSummedCountsInAdjacentWins",
           "_wordSize", wordName,
           "_winSize", winName,
           "_FDR", FDRname, ".pdf"),
    height = 15, width = 18)
par(mfrow = c(5, 2))
par(mar = c(2.1, 4.1, 2.1, 4.1))
par(mgp = c(3, 1, 0))

minSigPval <- NULL
maxSigPval <- NULL
minSigAdjPval <- NULL
maxSigAdjPval <- NULL
minCount <- NULL
maxCount <- NULL
for(i in 1:length(chrs)) {
  minSigPvalChr <- min(-log10(pmax(fisherPvalsSumCounts[[i]], 1e-323)))
  maxSigPvalChr <- max(-log10(pmax(fisherPvalsSumCounts[[i]], 1e-323)))
  minSigAdjPvalChr <- min(-log10(pmax(fisherAdjPvalsSumCounts[[i]], 1e-323)))
  maxSigAdjPvalChr <- max(-log10(pmax(fisherAdjPvalsSumCounts[[i]], 1e-323)))
  minCountChr <- min(c(sigSumWinCounts[[i]], sigSumRanCounts[[i]]))
  maxCountChr <- max(c(sigSumWinCounts[[i]], sigSumRanCounts[[i]]))

  minSigPval <- min(c(minSigPval, minSigPvalChr)) 
  maxSigPval <- max(c(maxSigPval, maxSigPvalChr))
  minSigAdjPval <- min(c(minSigAdjPval, minSigAdjPvalChr))
  maxSigAdjPval <- max(c(maxSigAdjPval, maxSigAdjPvalChr))
  minCount <- min(c(minCount, minCountChr))
  maxCount <- max(c(maxCount, maxCountChr))
}

for(i in 1:length(chrs)) {
  # Define adjacent windows of width winSize nt
  winStarts <- seq(from = 1,
                   to = length(genome[[i]]),
                   by = winSize)
  minusLog10PvalPlot(xplot = winStarts,
                     chrom = i,
                     word = "Summed word",
                     Pvals = pmax(fisherPvalsSumCounts[[i]], 1e-323), PvalsCol = "dodgerblue3",
                     minPval = minSigPval, maxPval = maxSigPval,
                     AdjPvals = pmax(fisherAdjPvalsSumCounts[[i]], 1e-323), AdjPvalsCol = "red",
                     minAdjPval = minSigAdjPval, maxAdjPval = maxSigAdjPval)
  countPlot(xplot = winStarts,
            chrom = i,
            minCount = minCount, maxCount = maxCount,
            word = "Summed word",
            winCounts = sigSumWinCounts[[i]], winCountsCol = "forestgreen",
            ranCounts = sigSumRanCounts[[i]], ranCountsCol = "grey50")
}
dev.off()
