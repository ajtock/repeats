#!/applications/R/R-3.5.0/bin/Rscript

# Forge a BSgenome data package for concatenated Cen180 sequences in
# nanopore contigs (concatenated by Ian Henderson)

# Usage:
# ./BSgenome_pkg_forge.R Cen180nanopore1

args <- commandArgs(trailingOnly = T)
packageNamePt4 <- args[1]

library(BSgenome)

seed_files <- system.file("extdata", "HendersonLab", package = "BSgenome") 
tail(list.files(seed_files, pattern = "-seed$"))

Cen180cat_seed <- list.files(seed_files,
                             pattern = paste0("\\.", packageNamePt4, "-seed$"),
                             full.names = TRUE)
cat(readLines(Cen180cat_seed), sep = "\n")

forgeBSgenomeDataPkg(Cen180cat_seed)

