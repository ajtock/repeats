#!/bin/bash

# Convert file containing feature coordinates in BED format
# to bedGraph format, and convert bedGraph to bigWig for use with
# hicPlotMatrix from HiCExplorer

# Usage:
# ./bed2bg2bw.sh CENAthila T2T_Col Chr1_Chr2_Chr3_Chr4_Chr5

prefix=$1
refbase=$2
chr=$3

source activate ChIPseq_mapping

awk 'BEGIN {OFS="\t"}; {print $1, $2, $3, 1}' ${prefix}_in_${refbase}_${chr}.bed > ${prefix}_in_${refbase}_${chr}.bedgraph
sed -i '1i track type=bedGraph' ${prefix}_in_${refbase}_${chr}.bedgraph
# USAGE: bedGraphToBigWig in.bedGraph chrom.sizes out.bw
# where in.bedGraph is a four-column file in the format:
#       <chrom> <start> <end> <value>
# and chrom.sizes is a two-column file/URL: <chromosome name> <size in bases>
# and out.bw is the output indexed big wig file.
# The input bedGraph file must be sorted, use the unix sort command:
sed '1d' ${prefix}_in_${refbase}_${chr}.bedgraph | LC_COLLATE=C sort -k1,1 -k2,2n > ${prefix}_in_${refbase}_${chr}.bedgraph_sorted
bedGraphToBigWig ${prefix}_in_${refbase}_${chr}.bedgraph_sorted /home/ajt200/analysis/nanopore/${refbase}/${refbase}.fa.sizes ${prefix}_in_${refbase}_${chr}.bw
rm ${prefix}_in_${refbase}_${chr}.bedgraph_sorted

conda deactivate
