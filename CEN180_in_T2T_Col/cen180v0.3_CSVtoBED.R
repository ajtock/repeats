#!/applications/R/R-3.5.0/bin/Rscript

# Convert CEN180 sequence coordinates identified in T2T_Col
# from CSV (1-based start coordinates) into BED format (0-based start coordinates)
# Convert coordinates corresponding to the intervening sequences between CEN180
# sequences from CSV into BED format

# Usage:
# ./cen180v0.3_CSVtoBED.R 'Chr1,Chr2,Chr3,Chr4,Chr5'

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
CENstart <- c(14840750, 3724530, 13597090, 4203495, 11783990)
CENend <- c(17558182, 5946091, 15733029, 6977107, 14551874)
CENGR <- GRanges(seqnames = chrs,
                 ranges = IRanges(start = CENstart,
                                  end = CENend),
                 strand = "*")

# Load table of CEN180 sequence coordinates
# with weighted SNVs vs consensus
tab <- read.csv("cen180v0.3.csv", header = T)
tab$chromosome <- gsub(pattern = "^", replacement = "Chr",
                       x = tab$chromosome)
CEN180GR <- GRanges(seqnames = tab$chromosome,
                    ranges = IRanges(start = tab$start,
                                     end = tab$end),
                    strand = tab$strand,
                    index = tab$index,
                    weightedSNV = tab$weightedSNV,
                    HORlengthsSum = tab$HORlengthsSum,
                    HORcount = tab$HORcount)
# Sort GRanges object (by chromosome, strand, start coordinate and, finally, end coordinate)
# Sorting data.frame by multiple columns (including that correspodning to strand)
# would be more complicated
# Necessary to include sorting by strand because intervening CEN180 sequences on the
# opposite strand would otherwise result in non-detection of tandem repeats on the same strand
CEN180GR <- sort(CEN180GR)
CEN180GR <- CEN180GR[seqnames(CEN180GR) %in% chrName]
CEN180_bed <- data.frame(chr = as.character(seqnames(CEN180GR)),
                         start = start(CEN180GR)-1,
                         end = end(CEN180GR),
                         name = CEN180GR$index,
                         score = CEN180GR$weightedSNV,
                         strand = as.character(strand(CEN180GR)),
                         HORlengthsSum = CEN180GR$HORlengthsSum,
                         HORcount = CEN180GR$HORcount)
write.table(CEN180_bed,
            file = paste0("CEN180_in_T2T_Col_",
                          paste0(chrName, collapse = "_"), ".bed"),
            quote = F, sep = "\t", row.names = F, col.names = F)

# Define function to select randomly positioned loci of the same
# width distribution as CEN180_bed
ranLocStartSelect <- function(coordinates, n) {
  sample(x = coordinates,
         size = n,
         replace = FALSE)
}

# Disable scientific notation (e.g., 59000000 rather than 5.9e+07)
options(scipen = 100)

# Define seed so that random selections are reproducible
set.seed(76492749)

# Apply ranLocStartSelect() on a per-chromosome basis so that
# ranLocGR contains the same number of loci per chromosome as CEN180GR
chrs <- chrs[chrs %in% chrName]
ranLocGR <- GRanges()
for(i in 1:length(chrs)) {
  CEN180ChrGR <- CEN180GR[seqnames(CEN180GR) == chrs[i]]
  regionChrGR <- regionGR[seqnames(regionGR) == chrs[i]]
  # Contract regionChrGR so that random loci and 2-kb flanking regions
  # do not extend beyond chromosome ends
  end(regionChrGR) <- end(regionChrGR)-max(width(CEN180ChrGR))-2000
  start(regionChrGR) <- start(regionChrGR)+2000
  ranLocChrStart <- ranLocStartSelect(coordinates = unlist(lapply(seq_along(regionChrGR), function(x) {           
                                                             start(regionChrGR[x]) : end(regionChrGR[x])          
                                                           })),
                                      n = length(CEN180ChrGR))
  ranLocChrGR <- GRanges(seqnames = chrs[i],
                         ranges = IRanges(start = ranLocChrStart,
                                          width = width(CEN180ChrGR)),
                         strand = strand(CEN180ChrGR))
  ranLocGR <- append(ranLocGR, ranLocChrGR)
}
stopifnot(identical(width(ranLocGR), width(CEN180GR)))
stopifnot(identical(as.character(seqnames(ranLocGR)), as.character(seqnames(CEN180GR))))
stopifnot(identical(strand(ranLocGR), strand(CEN180GR)))
ranLoc_bed <- data.frame(chr = as.character(seqnames(ranLocGR)),
                         start = start(ranLocGR)-1,
                         end = end(ranLocGR),
                         name = 1:length(ranLocGR),
                         score = "NA",
                         strand = strand(ranLocGR))
write.table(ranLoc_bed,
            file = paste0("CEN180_in_T2T_Col_",
                          paste0(chrName, collapse = "_"), "_randomLoci.bed"),
            quote = F, sep = "\t", row.names = F, col.names = F)

# Define adjacent centromeric loci each with widths equal to the median width of
# CEN180 sequences in (a) specified chromosome(s)
# Extract those that do not overlap CEN180 sequences
CENgapGR <- GRanges()
for(i in 1:length(chrs)) {
  CEN180ChrGR <- CEN180GR[seqnames(CEN180GR) == chrs[i]]
  CENChrGR <- CENGR[seqnames(CENGR) == chrs[i]]
  CENgapChrStart <- seq(from = start(CENChrGR),
                        to = end(CENChrGR),
                        by = median(width(CEN180ChrGR)))  
  CENgapChrGR <- GRanges(seqnames = chrs[i],
                         ranges = IRanges(start = CENgapChrStart,
                                          width = median(width(CEN180ChrGR))),
                         strand = "*")
  CENgapGR <- append(CENgapGR, CENgapChrGR)
}
CENgapGR_CEN180GR_ol <- findOverlaps(query = CEN180GR,
                                     subject = CENgapGR,
                                     type = "any",
                                     select = "all",
                                     ignore.strand = TRUE)
CENgapGR <- CENgapGR[-subjectHits(CENgapGR_CEN180GR_ol)]

# Sort GRanges object (by chromosome, strand, start coordinate and, finally, end coordinate)
CENgapGR <- sort(CENgapGR)
CENgapGR <- CENgapGR[seqnames(CENgapGR) %in% chrName]
CENgap_bed <- data.frame(chr = as.character(seqnames(CENgapGR)),
                         start = start(CENgapGR)-1,
                         end = end(CENgapGR),
                         name = 1:length(CENgapGR),
                         score = "NA",
                         strand = strand(CENgapGR))
write.table(CENgap_bed,
            file = paste0("CENgap_in_T2T_Col_",
                          paste0(chrName, collapse = "_"), ".bed"),
            quote = F, sep = "\t", row.names = F, col.names = F)
