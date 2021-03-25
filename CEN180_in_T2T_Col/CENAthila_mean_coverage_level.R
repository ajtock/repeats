#!/applications/R/R-4.0.0/bin/Rscript

# author: Andy Tock
# contact: ajt200@cam.ac.uk
# date: 24.02.2021

# Calculate mean levels between each CENAthila start and end coordinate

# Usage:
# /applications/R/R-4.0.0/bin/Rscript CENAthila_mean_coverage_level.R 'Chr1,Chr2,Chr3,Chr4,Chr5' both 180 2000 2000 2000 '2kb' 10 10 10bp 10bp

#chrName <- unlist(strsplit("Chr1,Chr2,Chr3,Chr4,Chr5",
#                           split = ","))
#align <- "both"
#bodyLength <- 180
#Athila_bodyLength <- 2000
#TEsf_bodyLength <- 2000
#upstream <- 2000
#downstream <- 2000
#flankName <- "2kb"
#binSize <- 10
#Athila_binSize <- 10
#binName <- "10bp"
#Athila_binName <- "10bp"

args <- commandArgs(trailingOnly = T)
chrName <- unlist(strsplit(args[1],
                           split = ","))
align <- args[2]
bodyLength <- as.numeric(args[3])
Athila_bodyLength <- as.numeric(args[4])
TEsf_bodyLength <- as.numeric(args[5])
upstream <- as.numeric(args[6])
downstream <- as.numeric(args[6])
flankName <- args[7]
binSize <- as.numeric(args[8])
Athila_binSize <- as.numeric(args[9])
binName <- args[10]
Athila_binName <- args[11]

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
RNANames <- unlist(strsplit("WT_RNAseq_Rep1,met1_RNAseq_Rep1",
                             split = ","))
RNANamesDir <- unlist(strsplit("160601_Kyuha_RNAseq/snakemake_RNAseq_STAR_T2T_Col,160601_Kyuha_RNAseq/snakemake_RNAseq_STAR_T2T_Col",
                                split = ","))
RNANamesPlot <- unlist(strsplit("wt_RNA,met1_RNA",
                                 split = ","))

RNADirs <- sapply(seq_along(RNANames), function(x) {
  paste0("/home/ajt200/analysis/",
         RNANamesDir[x],
         "/mapped/")
})

RNA_AthilaMats <- mclapply(seq_along(RNANames), function(x) {
    as.matrix(read.table(paste0(RNADirs[x], "CEN180profiles/matrices/",
                                RNANames[x],
                                "_MappedOn_T2T_Col_", align, "_sort_norm_CENAthila_in_",
                                paste0(chrName, collapse = "_"), "_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                         header = F, skip = 3))
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

# sRNA
sRNANames <- unlist(strsplit("Col_0_sRNA_ERR966148,met1_3_sRNA_ERR966149",
                             split = ","))
sRNANamesDir <- unlist(strsplit("sRNAseq_leaf_Rigal_Mathieu_2016_PNAS/snakemake_sRNAseq_T2T_Col,sRNAseq_leaf_Rigal_Mathieu_2016_PNAS/snakemake_sRNAseq_T2T_Col",
                                split = ","))
sRNANamesPlot <- unlist(strsplit("wt,met1",
                                 split = ","))
sRNAsize <- unlist(strsplit("21,22,24",
                            split = ","))
sRNANames <- rep(sRNANames, 3)
sRNANamesDir <- rep(sRNANamesDir, 3)
sRNAsize <- as.character(sapply(sRNAsize, function(x) rep(x, 2)))
sRNANamesPlot <- paste0(rep(sRNANamesPlot, 3), "_sRNA_", sRNAsize, "nt")

sRNADirs <- sapply(seq_along(sRNANames), function(x) {
  paste0("/home/ajt200/analysis/",
         sRNANamesDir[x],
         "/mapped/")
})

sRNA_AthilaMats <- mclapply(seq_along(sRNANames), function(x) {
    as.matrix(read.table(paste0(sRNADirs[x], "CEN180profiles/matrices/",
                                sRNANames[x],
                                "_MappedOn_T2T_Col_", align, "_", sRNAsize[x], "nt_sort_norm_CENAthila_in_",
                                paste0(chrName, collapse = "_"), "_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                         header = F, skip = 3))
}, mc.cores = length(sRNANames))

# Add column names
for(x in seq_along(sRNA_AthilaMats)) {
  colnames(sRNA_AthilaMats[[x]]) <- c(paste0("u", 1:(upstream/Athila_binSize)),
                                      paste0("t", ((upstream/Athila_binSize)+1):((upstream+Athila_bodyLength)/Athila_binSize)),
                                      paste0("d", (((upstream+Athila_bodyLength)/Athila_binSize)+1):(((upstream+Athila_bodyLength)/Athila_binSize)+(downstream/Athila_binSize))))
}

sRNA_AthilaMats_bodies <- lapply(seq_along(sRNA_AthilaMats), function(x) {
  sRNA_AthilaMats[[x]][,((upstream/Athila_binSize)+1):((upstream+Athila_bodyLength)/Athila_binSize)]
})
sRNA_AthilaMats_bodiesRowMeans <- lapply(seq_along(sRNA_AthilaMats_bodies), function(x) {
  rowMeans(sRNA_AthilaMats_bodies[[x]], na.rm = T)
})
sRNA_AthilaMats_bodiesRowSums <- lapply(seq_along(sRNA_AthilaMats_bodies), function(x) {
  rowSums(sRNA_AthilaMats_bodies[[x]], na.rm = T)
})
sRNA_AthilaMats_bodiesRowMeans_mat <- do.call(cbind, sRNA_AthilaMats_bodiesRowMeans)
colnames(sRNA_AthilaMats_bodiesRowMeans_mat) <- paste0(sRNANamesPlot, "_mean")
sRNA_AthilaMats_bodiesRowSums_mat <- do.call(cbind, sRNA_AthilaMats_bodiesRowSums)
colnames(sRNA_AthilaMats_bodiesRowSums_mat) <- paste0(sRNANamesPlot, "_sum")

# DNAmeth
DNAmethNames <- unlist(strsplit("Col_0_BSseq_Rep1_ERR965674,Col_0_BSseq_Rep2_ERR965675,met1_3_BSseq_Rep1_ERR965676,met1_3_BSseq_Rep2_ERR965677",
                             split = ","))
DNAmethNamesDir <- unlist(strsplit("BSseq_leaf_Rigal_Mathieu_2016_PNAS/snakemake_BSseq_T2T_Col,BSseq_leaf_Rigal_Mathieu_2016_PNAS/snakemake_BSseq_T2T_Col,BSseq_leaf_Rigal_Mathieu_2016_PNAS/snakemake_BSseq_T2T_Col,BSseq_leaf_Rigal_Mathieu_2016_PNAS/snakemake_BSseq_T2T_Col",
                                split = ","))
DNAmethNamesPlot <- unlist(strsplit("wt_Rep1,wt_Rep2,met1_Rep1,met1_Rep2",
                                 split = ","))
context <- c("CpG", "CHG", "CHH")
DNAmethNames <- rep(DNAmethNames, 3)
DNAmethNamesDir <- rep(DNAmethNamesDir, 3)
context <- as.character(sapply(context, function(x) rep(x, 4)))
DNAmethNamesPlot <- paste0(rep(DNAmethNamesPlot, 3), "_", context)

DNAmethDirs <- sapply(seq_along(DNAmethNamesDir), function(x) {
  if(grepl("nanopolish", DNAmethNamesDir[x])) {
    paste0("/home/ajt200/analysis/",
           DNAmethNamesDir[x],
           "/")
  } else {
    paste0("/home/ajt200/analysis/",
           DNAmethNamesDir[x],
           "/coverage/")
  }
})

DNAmeth_AthilaMats <- mclapply(seq_along(DNAmethNames), function(x) {
    as.matrix(read.table(paste0(DNAmethDirs[x], "CEN180profiles/matrices/",
                                DNAmethNames[x],
                                "_MappedOn_T2T_Col_", context[x], "_CENAthila_in_",
                                paste0(chrName, collapse = "_"), "_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                         header = F, skip = 3))
}, mc.cores = length(DNAmethNames))

# Add column names
for(x in seq_along(DNAmeth_AthilaMats)) {
  colnames(DNAmeth_AthilaMats[[x]]) <- c(paste0("u", 1:(upstream/Athila_binSize)),
                                         paste0("t", ((upstream/Athila_binSize)+1):((upstream+Athila_bodyLength)/Athila_binSize)),
                                         paste0("d", (((upstream+Athila_bodyLength)/Athila_binSize)+1):(((upstream+Athila_bodyLength)/Athila_binSize)+(downstream/Athila_binSize))))
}

DNAmeth_AthilaMats_bodies <- lapply(seq_along(DNAmeth_AthilaMats), function(x) {
  DNAmeth_AthilaMats[[x]][,((upstream/Athila_binSize)+1):((upstream+Athila_bodyLength)/Athila_binSize)]
})
DNAmeth_AthilaMats_bodiesRowMeans <- lapply(seq_along(DNAmeth_AthilaMats_bodies), function(x) {
  rowMeans(DNAmeth_AthilaMats_bodies[[x]], na.rm = T)
})
DNAmeth_AthilaMats_bodiesRowMeans_mat <- do.call(cbind, DNAmeth_AthilaMats_bodiesRowMeans)
colnames(DNAmeth_AthilaMats_bodiesRowMeans_mat) <- paste0(DNAmethNamesPlot, "_mean")


# Load table of centromeric Athila element coordinates for addition of coverage rowMeans and rowSums columns
tab <- read.table(paste0("/home/ajt200/analysis/nanopore/T2T_Col/annotation/TEs/",
                         "t2t_Athila_master.txt_coord.corrected_50Athila_strand.bed.clean.matrix"),
                  header = T, na.strings = "na")
colnames(tab)[1] <- "chr"
tab$chr <- gsub(pattern = "^", replacement = "Chr",
                x = tab$chr)
tab_extend <- data.frame(tab,
                         RNA_AthilaMats_bodiesRowMeans_mat,
                         sRNA_AthilaMats_bodiesRowMeans_mat,
                         DNAmeth_AthilaMats_bodiesRowMeans_mat)
#                         RNA_AthilaMats_bodiesRowSums_mat, 
#                         sRNA_AthilaMats_bodiesRowSums_mat) 


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
          column_title = sub(" \\(.+\\)", "", datName),
          column_title_rot = 0,
          column_title_gp = gpar(fontsize = 13),
          column_labels = sub("_.+", "", colnames(mat)),
          column_names_rot = 0,
          column_names_centered = TRUE,
          cluster_columns = FALSE,
          cluster_column_slices = FALSE,
          cluster_rows = FALSE,
          cluster_row_slices = FALSE,
          heatmap_legend_param = list(title = sub(" \\(.+\\)", "", datName),
                                      title_position = "topcenter",
                                      title_gp = gpar(font = 2, fontsize = 12),
                                      legend_direction = "horizontal",
                                      labels_gp = gpar(fontsize = 10)),
          heatmap_width = unit(2, "npc"),
          heatmap_height = unit(4, "npc"),
          column_gap = unit(0, "mm"),
          row_gap = unit(1.0, "mm"),
          row_split = factor(tab_extend$class1,
                             levels = unique(as.character(tab_extend$class1))),
          row_title = NULL,
          show_row_names = FALSE,
          border = FALSE,
          # If converting into png with pdfTotiffTopng.sh,
          # set use_raster to FALSE
          use_raster = FALSE)
          #use_raster = TRUE, raster_device = "png", raster_quality = 4)
}

RNA_mat <- as.matrix(tab_extend[,which(grepl("_RNA", colnames(tab_extend)))])
sRNA_21nt_mat <- as.matrix(tab_extend[,which(grepl("_sRNA_21nt", colnames(tab_extend)))])
sRNA_22nt_mat <- as.matrix(tab_extend[,which(grepl("_sRNA_22nt", colnames(tab_extend)))])
sRNA_24nt_mat <- as.matrix(tab_extend[,which(grepl("_sRNA_24nt", colnames(tab_extend)))])
CpG_mat <- as.matrix(tab_extend[,which(grepl("_CpG_", colnames(tab_extend)))])
CHG_mat <- as.matrix(tab_extend[,which(grepl("_CHG_", colnames(tab_extend)))])
CHH_mat <- as.matrix(tab_extend[,which(grepl("_CHH_", colnames(tab_extend)))])

RNA_colFun <- colorRamp2(quantile(RNA_mat,
                                  c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95),
                                  na.rm = T),
                         heat.colors(6))
sRNA_21nt_colFun <- colorRamp2(quantile(sRNA_21nt_mat,
                                        c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95),
                                        na.rm = T),
                               plasma(6))
sRNA_22nt_colFun <- colorRamp2(quantile(sRNA_22nt_mat,
                                        c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95),
                                        na.rm = T),
                               plasma(6))
sRNA_24nt_colFun <- colorRamp2(quantile(sRNA_24nt_mat,
                                        c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95),
                                        na.rm = T),
                               plasma(6))
CpG_colFun <- colorRamp2(quantile(CpG_mat,
                                  c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95),
                                  na.rm = T),
                         viridis(6))
CHG_colFun <- colorRamp2(quantile(CHG_mat,
                                  c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95),
                                  na.rm = T),
                         viridis(6))
CHH_colFun <- colorRamp2(quantile(CHH_mat,
                                  c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95),
                                  na.rm = T),
                         viridis(6))

RNA_htmp <- featureHeatmap(mat = RNA_mat,
                           colFun = RNA_colFun,
                           datName = "RNA-seq (TPM)",
                           rowOrder = c(1:nrow(RNA_mat)))
sRNA_21nt_htmp <- featureHeatmap(mat = sRNA_21nt_mat,
                                 colFun = sRNA_21nt_colFun,
                                 datName = "21-nt sRNAs (TPM)",
                                 rowOrder = c(1:nrow(sRNA_21nt_mat)))
sRNA_22nt_htmp <- featureHeatmap(mat = sRNA_22nt_mat,
                                 colFun = sRNA_22nt_colFun,
                                 datName = "22-nt sRNAs (TPM)",
                                 rowOrder = c(1:nrow(sRNA_22nt_mat)))
sRNA_24nt_htmp <- featureHeatmap(mat = sRNA_24nt_mat,
                                 colFun = sRNA_24nt_colFun,
                                 datName = "24-nt sRNAs (TPM)",
                                 rowOrder = c(1:nrow(sRNA_24nt_mat)))
CpG_htmp <- featureHeatmap(mat = CpG_mat,
                           colFun = CpG_colFun,
                           datName = "CpG (%)",
                           rowOrder = c(1:nrow(CpG_mat)))
CHG_htmp <- featureHeatmap(mat = CHG_mat,
                           colFun = CHG_colFun,
                           datName = "CHG (%)",
                           rowOrder = c(1:nrow(CHG_mat)))
CHH_htmp <- featureHeatmap(mat = CHH_mat,
                           colFun = CHH_colFun,
                           datName = "CHH (%)",
                           rowOrder = c(1:nrow(CHH_mat)))

fam_htmp <- Heatmap(mat = matrix(tab_extend$class1),
                    col = c("Athila" = rich6()(6)[6], "Athila2" = rich6()(6)[5], "Athila5" = rich6()(6)[4],
                            "Athila6A" = rich6()(6)[3], "Athila6B" = rich6()(6)[2], "ATGP2_AT2TE00075" = rich6()(6)[1]),
                    row_order = c(1:nrow(matrix(tab_extend$class1))),
                    column_title = "Family",
                    column_title_rot = 0,
                    column_title_gp = gpar(fontsize = 13),
                    column_labels = NULL,
                    cluster_rows = FALSE,
                    cluster_row_slices = FALSE,
                    heatmap_legend_param = list(title = "Family",
                                                title_position = "topcenter",
                                                title_gp = gpar(font = 2, fontsize = 12),
                                                legend_direction = "horizontal",
                                                nrow = 1,
                                                labels_gp = gpar(fontsize = 10)),
                    heatmap_width = unit(10, "mm"),
                    heatmap_height = unit(4, "npc"),
                    column_gap = unit(0, "mm"),
                    row_gap = unit(1.0, "mm"),
                    row_split = factor(tab_extend$class1,
                                       levels = unique(as.character(tab_extend$class1))),
                    row_title_side = "right",
                    row_title_rot = 0,
                    border = FALSE,
                    # If converting into png with pdfTotiffTopng.sh,
                    # set use_raster to FALSE
                    use_raster = FALSE)
                    #use_raster = TRUE, raster_device = "png", raster_quality = 4)

htmps <- RNA_htmp + sRNA_21nt_htmp + sRNA_22nt_htmp + sRNA_24nt_htmp + CpG_htmp + CHG_htmp + CHH_htmp + fam_htmp 

pdf(paste0(plotDir,
           "CENAthila_in_T2T_Col_",
           paste0(chrName, collapse = "_"), "_mean_RNAseq_sRNAseq_DNAmeth_heatmap.pdf"),
    width = 1.5*length(htmps), height = 8)
draw(htmps,
     heatmap_legend_side = "bottom",
     gap = unit(c(1), "mm"))
dev.off()

write.table(tab_extend,
            file = paste0("/home/ajt200/analysis/nanopore/T2T_Col/annotation/TEs/",
                          "t2t_Athila_master.txt_coord.corrected_50Athila_strand_with_mean_RNA_sRNA_DNAmeth.bed.clean.matrix"),
            quote = F, sep = "\t", row.names = F, col.names = T)
