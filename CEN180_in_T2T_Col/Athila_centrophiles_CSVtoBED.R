#!/applications/R/R-3.5.0/bin/Rscript

# Convert centromeric Athila LTR element coordinates identified in T2T_Col
# from CSV (1-based start coordinates) to BED format (0-based start coordinates)

# Usage:
# ./Athila_centrophiles_CSVtoBED.R 'Chr1,Chr2,Chr3,Chr4,Chr5'

#chrName <- unlist(strsplit("Chr1,Chr2,Chr3,Chr4,Chr5",
#                           split = ","))

args <- commandArgs(trailingOnly = T)
chrName <- unlist(strsplit(args[1],
                           split = ","))

options(stringsAsFactors = F)
library(GenomicRanges)

# Genomic definitions
fai <- read.table("T2T_Col.fa.fai", header = F)
chrs <- fai$V1[1:5]
chrLens <- fai$V2[1:5]
regionGR <- GRanges(seqnames = chrs,
                    ranges = IRanges(start = 1,
                                     end = chrLens),
                    strand = "*")

# Load table of Athila sequence coordinates
# with weighted SNVs vs consensus
tab <- read.csv("athila.rug.csv", header = T)
tab$chr <- gsub(pattern = "^", replacement = "Chr",
                x = tab$chr)
tab <- tab[!is.na(tab$ltr1_start.1),]
AthilaGR <- GRanges(seqnames = tab$chr,
                    ranges = IRanges(start = tab$ltr1_start.1,
                                     end = tab$ltr2_end.1),
                    strand = "*",
                    name = tab$name)
# Sort GRanges object (by chromosome, strand, start coordinate and, finally, end coordinate)
AthilaGR <- sort(AthilaGR)
AthilaGR <- AthilaGR[seqnames(AthilaGR) %in% chrName]
Athila_bed <- data.frame(chr = as.character(seqnames(AthilaGR)),
                         start = as.integer(start(AthilaGR)-1),
                         end = as.integer(end(AthilaGR)),
                         name = as.character(AthilaGR$name),
                         score = "NA",
                         strand = as.character(strand(AthilaGR)))
write.table(Athila_bed,
            file = paste0("CENAthila_in_T2T_Col_",
                          paste0(chrName, collapse = "_"), ".bed"),
            quote = F, sep = "\t", row.names = F, col.names = F)
