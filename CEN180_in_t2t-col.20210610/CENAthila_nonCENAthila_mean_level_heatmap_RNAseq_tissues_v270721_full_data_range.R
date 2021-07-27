#!/applications/R/R-4.0.0/bin/Rscript

# author: Andy Tock
# contact: ajt200@cam.ac.uk
# date: 14.07.2021

# Calculate mean levels between each CENAthila and nonCENAthila start and end coordinate,
# and plot as heatmaps

# Usage:
# /applications/R/R-4.0.0/bin/Rscript CENAthila_nonCENAthila_mean_level_heatmap_RNAseq_tissues_v270721_full_data_range.R 'Chr1,Chr2,Chr3,Chr4,Chr5' both 2000 2000 '2kb' 10

#chrName <- unlist(strsplit("Chr1,Chr2,Chr3,Chr4,Chr5",
#                           split = ","))
#align <- "both"
#Athila_bodyLength <- 2000
#upstream <- 2000
#downstream <- 2000
#flankName <- "2kb"
#Athila_binSize <- 10

args <- commandArgs(trailingOnly = T)
chrName <- unlist(strsplit(args[1],
                           split = ","))
align <- args[2]
Athila_bodyLength <- as.numeric(args[3])
upstream <- as.numeric(args[4])
downstream <- as.numeric(args[4])
flankName <- args[5]
Athila_binSize <- as.numeric(args[6])

options(stringsAsFactors = F)
library(parallel)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(grid)
library(gridExtra)
library(extrafont)
library(ComplexHeatmap)
library(dplyr)
library(RColorBrewer)
library(viridis)
library(scales)
library(circlize)

outDir <- paste0(paste0(chrName, collapse = "_"), "/")
plotDir <- paste0(outDir, "plots/")
system(paste0("[ -d ", outDir, " ] || mkdir -p ", outDir))
system(paste0("[ -d ", plotDir, " ] || mkdir -p ", plotDir))

# Load feature matrices for each data set
# RNA
RNANamesDir <- "RNAseq_tissues_Klepikova_Penin_2016_PlantJ/snakemake_RNAseq_STAR_t2t-col.20210610"
RNANames <- system(paste0("ls /home/ajt200/analysis/", RNANamesDir, "/mapped/both/bw"), intern = T)
RNANames <- gsub(pattern = "_MappedOn_.+", replacement = "", x = RNANames)
RNANamesPlot <- gsub(pattern = "Col_0_", replacement = "", x = RNANames)
#RNANamesPlot <- gsub(pattern = "_RNAseq_.+", replacement = "", x = RNANamesPlot)
#RNANamesPlot <- gsub(pattern = "_", replacement = " ", x = RNANamesPlot)
substr(RNANamesPlot, start = 1, stop = 1) <- toupper(substr(RNANamesPlot, start = 1, stop = 1))

RNADirs <- sapply(seq_along(RNANames), function(x) {
  paste0("/home/ajt200/analysis/",
         RNANamesDir,
         "/mapped/")
})

RNA_AthilaMats <- mclapply(seq_along(RNANames), function(x) {
    as.matrix(rbind(read.table(paste0(RNADirs[x], "CENAthilaProfiles/matrices/",
                                      RNANames[x],
                                      "_MappedOn_t2t-col.20210610_", align, "_sort_norm_CENAthila_in_",
                                      paste0(chrName, collapse = "_"), "_matrix_bin", Athila_binSize, "bp_flank", flankName, ".tab"),
                               header = F, skip = 3),
                    read.table(paste0(RNADirs[x], "CENAthilaProfiles/matrices/",
                                      RNANames[x],
                                      "_MappedOn_t2t-col.20210610_", align, "_sort_norm_nonCENAthila_in_",
                                      paste0(chrName, collapse = "_"), "_matrix_bin", Athila_binSize, "bp_flank", flankName, ".tab"),
                               header = F, skip = 3)))
}, mc.cores = length(RNANames))

# Add column names
for(x in seq_along(RNA_AthilaMats)) {
  colnames(RNA_AthilaMats[[x]]) <- c(paste0("u", 1:(upstream/Athila_binSize)),
                                     paste0("t", ((upstream/Athila_binSize)+1):((upstream+Athila_bodyLength)/Athila_binSize)),
                                     paste0("d", (((upstream+Athila_bodyLength)/Athila_binSize)+1):(((upstream+Athila_bodyLength)/Athila_binSize)+(downstream/Athila_binSize))))
}

RNA_AthilaMats_bodies <- lapply(seq_along(RNA_AthilaMats), function(x) {
  RNA_AthilaMats[[x]][,((upstream/Athila_binSize)+1):((upstream+Athila_bodyLength)/Athila_binSize)]
})
RNA_AthilaMats_bodiesRowMeans <- lapply(seq_along(RNA_AthilaMats_bodies), function(x) {
  rowMeans(RNA_AthilaMats_bodies[[x]], na.rm = T)
})
RNA_AthilaMats_bodiesRowSums <- lapply(seq_along(RNA_AthilaMats_bodies), function(x) {
  rowSums(RNA_AthilaMats_bodies[[x]], na.rm = T)
})
RNA_AthilaMats_bodiesRowMeans_mat <- do.call(cbind, RNA_AthilaMats_bodiesRowMeans)
colnames(RNA_AthilaMats_bodiesRowMeans_mat) <- paste0(RNANamesPlot, "_mean")
RNA_AthilaMats_bodiesRowSums_mat <- do.call(cbind, RNA_AthilaMats_bodiesRowSums)
colnames(RNA_AthilaMats_bodiesRowSums_mat) <- paste0(RNANamesPlot, "_sum")


## Load table of centromeric gaps (Athila-containing discontinuities in CEN180 sequences)
# and Athila sequence coordinates
tab <- read.table("/home/ajt200/analysis/nanopore/t2t-col.20210610/annotation/t2t-col.20210610_ATHILA.inout.FL.soloLTR.fragments.tsv",
                  header = T)
tab$phylo <- toupper(tab$phylo)
tab$centr_pos <- paste0(tab$centr_pos, "_centr")
tab$phylo <- gsub(pattern = "ATHILA4C", replacement = "ATHILA4", x = tab$phylo)
tab$phylo <- gsub(pattern = "ATHILA4-4C", replacement = "ATHILA4", x = tab$phylo)

CENATHILA <- tab[tab$centr_pos == "in_centr" & tab$quality == "intact",]
nonCENATHILA <- tab[tab$centr_pos == "out_centr" & tab$quality == "intact",]
CENATHILA$phylo <- gsub("^", "CEN", CENATHILA$phylo)
nonCENATHILA$phylo <- gsub("^", "nonCEN", nonCENATHILA$phylo)
CENATHILA <- CENATHILA[ with(CENATHILA, order(chr, genome_left_coord_FL, genome_right_coord_FL, decreasing = F)), ]
nonCENATHILA <- nonCENATHILA[ with(nonCENATHILA, order(chr, genome_left_coord_FL, genome_right_coord_FL, decreasing = F)), ]

stopifnot(identical(colnames(CENATHILA), colnames(nonCENATHILA)))

CENATHILA <- data.frame(chr = CENATHILA$chr,
#                        gap_start = CENATHILA$gap_start,
#                        gap_end = CENATHILA$gap_end,
#                        gap_width = CENATHILA$gap_width,
#                        gap_name = CENATHILA$gap_name,
                        phylo = CENATHILA$phylo,
                        start = CENATHILA$genome_left_coord_FL,
                        end = CENATHILA$genome_right_coord_FL,
#                        length_FL = CENATHILA$length_FL,
                        strand = CENATHILA$direction,
                        name = CENATHILA$TE_ID,
#                        left_TSD = CENATHILA$left_TSD,
#                        right_TSD = CENATHILA$right_TSD,
#                        TSD = CENATHILA$TSD,
#                        quality = CENATHILA$quality,
#                        genome_5LTR_start = CENATHILA$genome_5LTR_start,
#                        genome_5LTR_end = CENATHILA$genome_5LTR_end,
#                        genome_3LTR_start = CENATHILA$genome_3LTR_start,
#                        genome_3LTR_end = CENATHILA$genome_3LTR_end,
#                        internal.length = abs(CENATHILA$genome_3LTR_start-CENATHILA$genome_5LTR_end+1),
                        position = CENATHILA$centr_pos)

nonCENATHILA <- data.frame(chr = nonCENATHILA$chr,
#                           gap_start = nonCENATHILA$gap_start,
#                           gap_end = nonCENATHILA$gap_end,
#                           gap_width = nonCENATHILA$gap_width,
#                           gap_name = nonCENATHILA$gap_name,
                           phylo = nonCENATHILA$phylo,
                           start = nonCENATHILA$genome_left_coord_FL,
                           end = nonCENATHILA$genome_right_coord_FL,
#                           length_FL = nonCENATHILA$length_FL,
                           strand = nonCENATHILA$direction,
                           name = nonCENATHILA$TE_ID,
#                           left_TSD = nonCENATHILA$left_TSD,
#                           right_TSD = nonCENATHILA$right_TSD,
#                           TSD = nonCENATHILA$TSD,
#                           quality = nonCENATHILA$quality,
#                           genome_5LTR_start = nonCENATHILA$genome_5LTR_start,
#                           genome_5LTR_end = nonCENATHILA$genome_5LTR_end,
#                           genome_3LTR_start = nonCENATHILA$genome_3LTR_start,
#                           genome_3LTR_end = nonCENATHILA$genome_3LTR_end,
#                           internal.length = abs(nonCENATHILA$genome_3LTR_start-nonCENATHILA$genome_5LTR_end+1),
                           position = nonCENATHILA$centr_pos)

ATHILA <- rbind(CENATHILA, nonCENATHILA)
tab_extend <- data.frame(ATHILA,
                         RNA_AthilaMats_bodiesRowMeans_mat)

# Define heatmap colours
rich12 <- function() {manual_pal(values = c("#000040","#000093","#0020E9","#0076FF","#00B8C2","#04E466","#49FB25","#E7FD09","#FEEA02","#FFC200","#FF8500","#FF3300"))}
rich10 <- function() {manual_pal(values = c("#000041","#0000A9","#0049FF","#00A4DE","#03E070","#5DFC21","#F6F905","#FFD701","#FF9500","#FF3300"))}
rich8 <- function() {manual_pal(values = c("#000041","#0000CB","#0081FF","#02DA81","#80FE1A","#FDEE02","#FFAB00","#FF3300"))}
rich6 <- function() {manual_pal(values = c("#000043","#0033FF","#01CCA4","#BAFF12","#FFCC00","#FF3300"))}
rich8to6equal <- c("#0000CB", "#0081FF", "#87CEFA", "#FDEE02", "#FFAB00", "#FF3300")
revSpectralScale11 <- rev(brewer.pal(11, "Spectral"))
viridisScale6 <- viridis_pal()(6)

# Heatmap plotting function
featureHeatmap <- function(mat,
                           colFun,
                           datName,
                           rowOrder) {
  Heatmap(matrix = mat,
          col = colFun,
          row_order = rowOrder,
#          column_title = sub(" \\(.+\\)", "", datName),
          column_title = datName,
          column_title_rot = 0,
          column_title_gp = gpar(fontsize = 13),
          column_labels = sub("_RNAseq_.+", "", colnames(mat)),
          column_names_rot = 90,
          column_names_centered = FALSE,
          cluster_columns = FALSE,
          cluster_column_slices = FALSE,
          cluster_rows = FALSE,
          cluster_row_slices = FALSE,
          heatmap_legend_param = list(title = datName,
                                      title_position = "topcenter",
                                      title_gp = gpar(font = 2, fontsize = 12),
                                      legend_direction = "horizontal",
                                      labels_gp = gpar(fontsize = 10)),
          heatmap_width = unit(2, "npc"),
          heatmap_height = unit(4, "npc"),
          column_gap = unit(0, "mm"),
          row_gap = unit(1.0, "mm"),
          row_split = factor(tab_extend$phylo,
                             levels = sort(unique(as.character(tab_extend$phylo)))),
          row_title = NULL,
          show_row_names = FALSE,
          border = FALSE,
          # If converting into png with pdfTotiffTopng.sh,
          # set use_raster to FALSE
          use_raster = FALSE)
          #use_raster = TRUE, raster_device = "png", raster_quality = 4)
}

RNA_mat <- as.matrix(tab_extend[,which(grepl("_RNA", colnames(tab_extend)))])
#sRNA_21nt_mat <- as.matrix(tab_extend[,which(grepl("_sRNA_21nt", colnames(tab_extend)))])
#sRNA_22nt_mat <- as.matrix(tab_extend[,which(grepl("_sRNA_22nt", colnames(tab_extend)))])
#sRNA_24nt_mat <- as.matrix(tab_extend[,which(grepl("_sRNA_24nt", colnames(tab_extend)))])
#CpG_mat <- as.matrix(tab_extend[,which(grepl("_CpG_", colnames(tab_extend)))])
#CHG_mat <- as.matrix(tab_extend[,which(grepl("_CHG_", colnames(tab_extend)))])
#CHH_mat <- as.matrix(tab_extend[,which(grepl("_CHH_", colnames(tab_extend)))])

#RNA_colFun <- colorRamp2(quantile(RNA_mat,
#                                  c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95),
#                                  na.rm = T),
#                         heat.colors(6))
#sRNA_21nt_colFun <- colorRamp2(quantile(sRNA_21nt_mat,
#                                        c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95),
#                                        na.rm = T),
#                               plasma(6))
#sRNA_22nt_colFun <- colorRamp2(quantile(sRNA_22nt_mat,
#                                        c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95),
#                                        na.rm = T),
#                               plasma(6))
#sRNA_24nt_colFun <- colorRamp2(quantile(sRNA_24nt_mat,
#                                        c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95),
#                                        na.rm = T),
#                               plasma(6))
#CpG_colFun <- colorRamp2(quantile(CpG_mat,
#                                  c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95),
#                                  na.rm = T),
#                         viridis(6))
#CHG_colFun <- colorRamp2(quantile(CHG_mat,
#                                  c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95),
#                                  na.rm = T),
#                         viridis(6))
#CHH_colFun <- colorRamp2(quantile(CHH_mat,
#                                  c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95),
#                                  na.rm = T),
#                         viridis(6))

RNA_colFun <- colorRamp2(quantile(RNA_mat,
                                  c(0.01, 0.2, 0.4, 0.6, 0.8, 0.99),
                                  na.rm = T),
                         plasma(6))
#                         heat.colors(6))
log2RNA_colFun <- colorRamp2(quantile(log2(RNA_mat+1),
                                      c(0.01, 0.2, 0.4, 0.6, 0.8, 0.99),
                                      na.rm = T),
                             plasma(6))
#                             heat.colors(6))
#sRNA_21nt_colFun <- colorRamp2(quantile(sRNA_21nt_mat,
#                                        c(0.01, 0.2, 0.4, 0.6, 0.8, 0.99),
#                                        na.rm = T),
#                               plasma(6))
#sRNA_22nt_colFun <- colorRamp2(quantile(sRNA_22nt_mat,
#                                        c(0.01, 0.2, 0.4, 0.6, 0.8, 0.99),
#                                        na.rm = T),
#                               plasma(6))
#sRNA_24nt_colFun <- colorRamp2(quantile(sRNA_24nt_mat,
#                                        c(0.01, 0.2, 0.4, 0.6, 0.8, 0.99),
#                                        na.rm = T),
#                               plasma(6))
#CpG_colFun <- colorRamp2(quantile(CpG_mat,
#                                  c(0.01, 0.2, 0.4, 0.6, 0.8, 0.99),
#                                  na.rm = T),
#                         viridis(6))
#CHG_colFun <- colorRamp2(quantile(CHG_mat,
#                                  c(0.01, 0.2, 0.4, 0.6, 0.8, 0.99),
#                                  na.rm = T),
#                         viridis(6))
#CHH_colFun <- colorRamp2(quantile(CHH_mat,
#                                  c(0.01, 0.2, 0.4, 0.6, 0.8, 0.99),
#                                  na.rm = T),
#                         viridis(6))

RNA_htmp <- featureHeatmap(mat = RNA_mat,
                           colFun = RNA_colFun,
                           datName = "RNA-seq (TPM)",
                           rowOrder = c(1:nrow(RNA_mat)))
log2RNA_htmp <- featureHeatmap(mat = log2(RNA_mat+1),
                               colFun = log2RNA_colFun,
                               datName = "Log2(RNA-seq (TPM)+1)",
                               rowOrder = c(1:nrow(RNA_mat)))
#sRNA_21nt_htmp <- featureHeatmap(mat = sRNA_21nt_mat,
#                                 colFun = sRNA_21nt_colFun,
#                                 datName = "21-nt (TPM)",
#                                 rowOrder = c(1:nrow(sRNA_21nt_mat)))
#sRNA_22nt_htmp <- featureHeatmap(mat = sRNA_22nt_mat,
#                                 colFun = sRNA_22nt_colFun,
#                                 datName = "22-nt (TPM)",
#                                 rowOrder = c(1:nrow(sRNA_22nt_mat)))
#sRNA_24nt_htmp <- featureHeatmap(mat = sRNA_24nt_mat,
#                                 colFun = sRNA_24nt_colFun,
#                                 datName = "24-nt (TPM)",
#                                 rowOrder = c(1:nrow(sRNA_24nt_mat)))
#CpG_htmp <- featureHeatmap(mat = CpG_mat,
#                           colFun = CpG_colFun,
#                           datName = "mCpG (%)",
#                           rowOrder = c(1:nrow(CpG_mat)))
#CHG_htmp <- featureHeatmap(mat = CHG_mat,
#                           colFun = CHG_colFun,
#                           datName = "mCHG (%)",
#                           rowOrder = c(1:nrow(CHG_mat)))
#CHH_htmp <- featureHeatmap(mat = CHH_mat,
#                           colFun = CHH_colFun,
#                           datName = "mCHH (%)",
#                           rowOrder = c(1:nrow(CHH_mat)))

fam_htmp <- Heatmap(mat = matrix(tab_extend$phylo),
                    col = c("CENATHILA1" = rich10()(10)[9], "CENATHILA2" = rich10()(10)[8], "CENATHILA5" = rich10()(10)[4],
                            "CENATHILA6A" = rich10()(10)[3], "CENATHILA6B" = rich10()(10)[2],
                            "nonCENATHILA0" = rich10()(10)[10], "nonCENATHILA1" = rich10()(10)[9], "nonCENATHILA2" = rich10()(10)[8],
                            "nonCENATHILA3" = rich10()(10)[7], "nonCENATHILA4" = rich10()(10)[6],
                            #"nonCENATHILA4C" = rich10()(10)[5],
                            "nonCENATHILA5" = rich10()(10)[4], "nonCENATHILA6A" = rich10()(10)[3], "nonCENATHILA6B" = rich10()(10)[2],
                            "nonCENATHILA7A" = rich10()(10)[1]),
                    row_order = c(1:nrow(matrix(tab_extend$phylo))),
                    column_title = "Fam.",
                    column_title_rot = 0,
                    column_title_gp = gpar(fontsize = 13),
                    column_labels = NULL,
                    cluster_rows = FALSE,
                    cluster_row_slices = FALSE,
                    heatmap_legend_param = list(title = "Family",
                                                title_position = "topcenter",
                                                title_gp = gpar(font = 2, fontsize = 12),
                                                legend_direction = "vertical",
                                                nrow = length(unique(as.character(tab_extend$phylo))),
                                                labels_gp = gpar(fontsize = 10)),
                    heatmap_width = unit(10, "mm"),
                    heatmap_height = unit(4, "npc"),
                    column_gap = unit(0, "mm"),
                    row_gap = unit(1.0, "mm"),
                    row_split = factor(tab_extend$phylo,
                                       levels = sort(unique(as.character(tab_extend$phylo)))),
                    row_title = NULL,
                    row_title_side = "right",
                    row_title_rot = 0,
                    border = FALSE,
                    # If converting into png with pdfTotiffTopng.sh,
                    # set use_raster to FALSE
                    use_raster = FALSE)
                    #use_raster = TRUE, raster_device = "png", raster_quality = 4)

reg_htmp <- Heatmap(mat = matrix(tab_extend$phylo),
                    col = c("CENATHILA1" = "#008ebe", "CENATHILA2" = "#008ebe", "CENATHILA5" = "#008ebe",
                            "CENATHILA6A" = "#008ebe", "CENATHILA6B" = "#008ebe",
                            "nonCENATHILA0" = "#be3000", "nonCENATHILA1" = "#be3000", "nonCENATHILA2" = "#be3000",
                            "nonCENATHILA3" = "#be3000", "nonCENATHILA4" = "#be3000",
                            #"nonCENATHILA4C" = "#be3000",
                            "nonCENATHILA5" = "#be3000", "nonCENATHILA6A" = "#be3000", "nonCENATHILA6B" = "#be3000",
                            "nonCENATHILA7A" = "#be3000"),
                    row_order = c(1:nrow(matrix(tab_extend$phylo))),
                    column_title = "Reg.",
                    column_title_rot = 0,
                    column_title_gp = gpar(fontsize = 13),
                    column_labels = NULL,
                    cluster_rows = FALSE,
                    cluster_row_slices = FALSE,
                    heatmap_legend_param = list(title = "Region",
                                                title_position = "topcenter",
                                                title_gp = gpar(font = 2, fontsize = 12),
                                                legend_direction = "vertical",
                                                nrow = length(unique(as.character(tab_extend$phylo))),
                                                labels_gp = gpar(fontsize = 10)),
                    heatmap_width = unit(10, "mm"),
                    heatmap_height = unit(4, "npc"),
                    column_gap = unit(0, "mm"),
                    row_gap = unit(1.0, "mm"),
                    row_split = factor(tab_extend$phylo,
                                       levels = sort(unique(as.character(tab_extend$phylo)))),
                    row_title = NULL,
                    row_title_side = "right",
                    row_title_rot = 0,
                    border = FALSE,
                    # If converting into png with pdfTotiffTopng.sh,
                    # set use_raster to FALSE
                    use_raster = FALSE)
                    #use_raster = TRUE, raster_device = "png", raster_quality = 4)

#htmps <- reg_htmp + fam_htmp + RNA_htmp + sRNA_21nt_htmp + sRNA_22nt_htmp + sRNA_24nt_htmp + CpG_htmp + CHG_htmp + CHH_htmp

htmps <- reg_htmp + fam_htmp + RNA_htmp

legendGap <- unit(15, "mm")

pdf(paste0(plotDir,
           "CENATHILA_nonCENATHILA_in_t2t-col.20210610_",
           paste0(chrName, collapse = "_"), "_mean_RNAseq_tissues_Klepikova_2016_PlantJ_heatmap_colourQuantiles0.01to0.99_v270721.pdf"),
    width = 5.5*length(htmps), height = 10)
draw(htmps,
     gap = unit(1, "mm"),
     column_title = "ATHILA",
     column_title_gp = gpar(font = 2, fontsize = 16),
     heatmap_legend_side = "bottom",
     legend_gap = legendGap) 
dev.off()

htmps <- reg_htmp + fam_htmp + log2RNA_htmp

legendGap <- unit(15, "mm")

pdf(paste0(plotDir,
           "CENATHILA_nonCENATHILA_in_t2t-col.20210610_",
           paste0(chrName, collapse = "_"), "_log2_mean_RNAseq_tissues_Klepikova_2016_PlantJ_heatmap_colourQuantiles0.01to0.99_v270721.pdf"),
    width = 5.5*length(htmps), height = 10)
draw(htmps,
     gap = unit(1, "mm"),
     column_title = "ATHILA",
     column_title_gp = gpar(font = 2, fontsize = 16),
     heatmap_legend_side = "bottom",
     legend_gap = legendGap) 
dev.off()

write.table(tab_extend,
            file = paste0("/home/ajt200/analysis/nanopore/t2t-col.20210610/annotation/CENAthila/",
                          "t2t_ATHILA_with_mean_RNA_tissues_Klepikova_2016_PlantJ_v270721.tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)
