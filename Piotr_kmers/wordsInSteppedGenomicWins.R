#!/applications/R/R-3.5.0/bin/Rscript

# Count the occurrences of wordSize-nt sequences in each winSize-nt window with
# a step size of stepSize nt along each chromosome

# Usage:
# wordsInSteppedGenomicWins.R 3 1000000 1Mb 100000 100kb

#wordSize <- 3
#wordName <- as.character(wordSize)
#winSize <- 1000000
#winName <- "1Mb"
#stepSize <- 100000
#stepName <- "100kb"

args <- commandArgs(trailingOnly = T)
wordSize <- as.numeric(args[1])
wordName <- as.character(wordSize)
winSize <- as.numeric(args[2])
winName <- as.character(args[3])
stepSize <- as.numeric(args[4])
stepName <- as.character(args[5])

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

#words <- NULL
#for(j in 1:(length(genome[[i]])-wordSize+1)) {
#  word <- toString(genome[[i]][j:(j+wordSize-1)])
#  words <- c(words, word)
#}

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
  # Define sliding windows of width winSize nt,
  # with a step of stepSize nt
  winStarts <- seq(from = 1,
                   to = length(genome[[i]])-winSize,
                   by = stepSize)
  winEnds <- c(seq(from = winStarts[1]+winSize-1,
                   to = length(genome[[i]]),
                   by = stepSize),
               length(genome[[i]]))
  if(length(genome[[i]])-winStarts[length(winStarts)] >= winSize) {
    winStarts <- c(winStarts,
                   winStarts[length(winStarts)]+stepSize)
  }
  winStringSet <- DNAStringSet(x = genome[[i]],
                               start = winStarts,
                               end = winEnds)

  # Generate shuffled sequence strings in stepped windows
  nonOverlap1 <- winStringSet[[1]][1:stepSize]
  nonOverlap1Shuffle <- nonOverlap1[sample(length(nonOverlap1))]
  stepOverlap1 <- winStringSet[[1]][(stepSize+1):(length(winStringSet[[1]]))]
  stepOverlap1Shuffle <- stepOverlap1[sample(length(stepOverlap1))]
  ranString1 <- DNAStringSet(paste0(nonOverlap1Shuffle,
                                    stepOverlap1Shuffle))
  #ranString1 <- winStringSet[[1]][sample(length(winStringSet[[1]]))]
  ranStringSet <- DNAStringSet(ranString1)
  for(x in 2:length(winStringSet)) {
    # Obtain the last string in ranStringSet
    lastString <- ranStringSet[[length(ranStringSet)]]

    # Obtain and shuffle the sequence in winStringSet[[x]]
    # that does not overlap the sequence in winStringSet[[x-1]]
    nonOverlap <- winStringSet[[x]][(winSize-stepSize+1):(length(winStringSet[[x]]))]
    nonOverlapShuffle <- nonOverlap[sample(length(nonOverlap))]

    # Paste the last string in ranStringSet and
    # the shuffled non-overlapping sequence together
    ranString <- DNAStringSet(paste0(lastString[(stepSize+1):(length(lastString))],
                                     nonOverlapShuffle))
    ranStringSet <- append(ranStringSet, ranString)
  }

  # Count the occurrences of wordSize-nt sequences in
  # each winSize-nt window with a step size of stepSize nt
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
       file = paste0("./wordsInSteppedGenomicWins_Chr", i,
                     "_wordSize", wordName,
                     "_winSize", winName,
                     "_stepSize", stepName,
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
       file = paste0("./wordsInSteppedShuffledWins_Chr", i,
                     "_wordSize", wordName,
                     "_winSize", winName,
                     "_stepSize", stepName,
                     ".RData"))
}
