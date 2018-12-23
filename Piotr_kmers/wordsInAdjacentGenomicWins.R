#!/applications/R/R-3.5.0/bin/Rscript

# Count the occurrences of wordSize-nt sequences in
# each adjacent winSize-nt window

# Usage:
# wordsInAdjacentGenomicWins.R 3 1000000 1Mb

#wordSize <- 3
#wordName <- as.character(wordSize)
#winSize <- 1000000
#winName <- "1Mb"

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
library(BSgenome)
library(BSgenome.Athaliana.TAIR.TAIR9)
ls("package:BSgenome.Athaliana.TAIR.TAIR9")
genome <- BSgenome.Athaliana.TAIR.TAIR9
print(genome)

chrs <- c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5")

# Create a DNAStringSet object containing all wordSize-nt
# words in the genome, excluding duplicates
wordStringSet <- NULL
for(i in 1:5) { 
  wordStarts <- seq(from = 1,
                    to = length(genome[[i]]),
                    by = 1)
  wordStarts <- wordStarts[-wordStarts[(length(wordStarts)-wordSize+2):(length(wordStarts))]]
  wordEnds <- seq(from = wordStarts[1]+wordSize-1,
                  to = length(genome[[i]]),
                  by = 1)
  wordStringSetChr <- DNAStringSet(x = genome[[i]],
                                   start = wordStarts,
                                   end = wordEnds)
  wordStringSet <- append(wordStringSet, wordStringSetChr)
}
print(paste0(wordName,
             "-nt words in the genome, including duplicates:")) 
print(wordStringSet)
wordStringSet <- unique(wordStringSet)
print(paste0(wordName,
             "-nt words in the genome, excluding duplicates:")) 
print(wordStringSet)

foreach(i = 1:5) %dopar% {
  # Define adjacent windows of width winSize nt
  winStarts <- seq(from = 1,
                   to = length(genome[[i]]),
                   by = winSize)
  winEnds <- c(seq(from = winStarts[1]+winSize-1,
                   to = length(genome[[i]]),
                   by = winSize),
               length(genome[[i]]))
  winStringSet <- DNAStringSet(x = genome[[i]],
                               start = winStarts,
                               end = winEnds)

  
  # Generate shuffled sequence strings in adjacent windows
  ranStringSet <- endoapply(winStringSet, sample)

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
       file = paste0("./wordsInAdjacentGenomicWins_", chrs[i],
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
       file = paste0("./wordsInAdjacentShuffledWins_", chrs[i],
                     "_wordSize", wordName,
                     "_winSize", winName,
                     ".RData"))
}
