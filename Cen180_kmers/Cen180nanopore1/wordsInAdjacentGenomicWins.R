#!/applications/R/R-3.5.0/bin/Rscript

# Count the occurrences of wordSize-nt sequences in
# each adjacent winSize-nt window

# Usage:
# wordsInAdjacentGenomicWins.R 10 1000 1kb

#wordSize <- 10
#wordName <- as.character(wordSize)
#winSize <- 1000
#winName <- "1kb"

args <- commandArgs(trailingOnly = T)
wordSize <- as.numeric(args[1])
wordName <- as.character(wordSize)
winSize <- as.numeric(args[2])
winName <- as.character(args[3])

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("BSgenome", version = "3.8")
#BiocManager::install("BSgenome.Athaliana.TAIR.TAIR9", version = "3.8")
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

chrs <- c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5")

# Create a DNAStringSet object containing all wordSize-nt
# words in Cen180
wordStarts <- seq(from = 1,
                  to = length(Cen180[[1]]),
                  by = 1)
wordStarts <- wordStarts[-wordStarts[(length(wordStarts)-wordSize+2):
                                     (length(wordStarts))]]
wordEnds <- seq(from = wordStarts[1]+wordSize-1,
                to = length(Cen180[[1]]),
                by = 1)
wordStringSet <- DNAStringSet(x = Cen180[[1]],
                              start = wordStarts,
                              end = wordEnds)
wordStringSetRC <- DNAStringSet(x = reverseComplement(Cen180[[1]]),
                                start = wordStarts,
                                end = wordEnds)
wordStringSet <- append(wordStringSet, wordStringSetRC)

print(paste0(wordName,
             "-nt words in Cen180, including duplicates:")) 
print(wordStringSet)

# Retain words that appear more than once, and keep just one occurrence
wordStringSet <- unique(wordStringSet[duplicated(wordStringSet)])
print(paste0(wordName,
             "-nt words that appear more than once in Cen180,
             excluding duplicates:"))
print(wordStringSet)
#  A DNAStringSet instance of length 217292

# Count number of occurrences of each word in Cen180
Cen180wordCounts <- unlist(mclapply(seq_along(wordStringSet), function(x) {
    countPattern(pattern = toString(wordStringSet[[x]]),
                 subject = append(Cen180[[1]], reverseComplement(Cen180[[1]])),
                 max.mismatch = 0)
}, mc.cores = detectCores()))

# Retain words that occur >= 10 times
wordStringSet <- wordStringSet[Cen180wordCounts >= 10]
#  A DNAStringSet instance of length 14628

# Define adjacent windows of width winSize nt in Cen180
winStarts <- seq(from = 1,
                 to = length(Cen180[[1]]),
                 by = winSize)
winEnds <- c(seq(from = winStarts[1]+winSize-1,
                 to = length(Cen180[[1]]),
                 by = winSize),
             length(Cen180[[1]]))
winStringSet <- DNAStringSet(x = Cen180[[1]],
                             start = winStarts,
                             end = winEnds)
winStringSetRC <- reverseComplement(winStringSet)

# Generate shuffled sequence strings in adjacent windows
ranStringSet <- endoapply(winStringSet, sample)
ranStringSetRC <- endoapply(winStringSetRC, sample)
# Check that base frequency is the same for real and shuffled sequences
stopifnot(
  sum(unlist(mclapply(seq_along(winStringSet), function(x) {
    identical(
      alphabetFrequency(winStringSet[x]),
      alphabetFrequency(ranStringSet[x])
    )
  }, mc.cores = detectCores())))
== length(winStringSet))
stopifnot(
  sum(unlist(mclapply(seq_along(winStringSetRC), function(x) {
    identical(
      alphabetFrequency(winStringSetRC[x]),
      alphabetFrequency(ranStringSetRC[x])
    )
  }, mc.cores = detectCores())))
== length(winStringSetRC))
  

# Count the occurrences of wordSize-nt sequences in
# each winSize-nt window
# Windows containing real sequences
winCounts <- mclapply(seq_along(wordStringSet), function(x) {
  vcountPattern(pattern = toString(wordStringSet[[x]]),
                subject = winStringSet,
                max.mismatch = 0)
}, mc.cores = detectCores())

names(winCounts) <- mclapply(seq_along(winCounts), function(x) {
  toString(wordStringSet[x])
}, mc.cores = detectCores())

save(winCounts,
     file = paste0("./wordsInAdjacentCen180nanopore1Wins",
                   "_wordSize", wordName,
                   "_winSize", winName,
                   ".RData"))

# Windows containing shuffled sequences
ranCounts <- mclapply(seq_along(wordStringSet), function(x) {
  vcountPattern(pattern = toString(wordStringSet[[x]]),
                subject = ranStringSet,
                max.mismatch = 0)
}, mc.cores = detectCores())

names(ranCounts) <- mclapply(seq_along(ranCounts), function(x) {
  toString(wordStringSet[x])
}, mc.cores = detectCores())

save(ranCounts,
     file = paste0("./wordsInAdjacentCen180nanopore1ShuffledWins",
                   "_wordSize", wordName,
                   "_winSize", winName,
                   ".RData"))

## Reverse-complemented sequences in adjacent windows
# Count the occurrences of wordSize-nt sequences in
# each winSize-nt window
# Windows containing real sequences
winCountsRC <- mclapply(seq_along(wordStringSet), function(x) {
  vcountPattern(pattern = toString(wordStringSet[[x]]),
                subject = winStringSetRC,
                max.mismatch = 0)
}, mc.cores = detectCores())

names(winCountsRC) <- mclapply(seq_along(winCountsRC), function(x) {
  toString(wordStringSet[x])
}, mc.cores = detectCores())

save(winCountsRC,
     file = paste0("./wordsInAdjacentCen180nanopore1WinsRC",
                   "_wordSize", wordName,
                   "_winSize", winName,
                   ".RData"))

# Windows containing shuffled sequences
ranCountsRC <- mclapply(seq_along(wordStringSet), function(x) {
  vcountPattern(pattern = toString(wordStringSet[[x]]),
                subject = ranStringSetRC,
                max.mismatch = 0)
}, mc.cores = detectCores())

names(ranCountsRC) <- mclapply(seq_along(ranCountsRC), function(x) {
  toString(wordStringSet[x])
}, mc.cores = detectCores())

save(ranCountsRC,
     file = paste0("./wordsInAdjacentCen180nanopore1ShuffledWinsRC",
                   "_wordSize", wordName,
                   "_winSize", winName,
                   ".RData"))
