#!/bin/bash

# wSNV 4 quantiles
/applications/R/R-3.5.0/bin/Rscript CEN180_quantile_stat_density_mean_95CI_plot.R wSNV SNVs 'Chr1,Chr2,Chr3,Chr4,Chr5' wSNV 4 1.00 0.2 '%2.0f' '%3.1f' '%2.0f' 0.65
/applications/R/R-3.5.0/bin/Rscript CEN180_quantile_stat_density_mean_95CI_plot.R map_K150_E4_in_bodies 'Mappability (k=150 e=4)' 'Chr1,Chr2,Chr3,Chr4,Chr5' wSNV 4 1.00 0.2 '%1.2f' '%3.1f' '%1.2f' 0.65
/applications/R/R-3.5.0/bin/Rscript CEN180_quantile_stat_density_mean_95CI_plot.R CENH3_in_bodies 'CENH3' 'Chr1,Chr2,Chr3,Chr4,Chr5' wSNV 4 1.00 0.2 '%1.1f' '%3.1f' '%1.1f' 0.65
/applications/R/R-3.5.0/bin/Rscript CEN180_quantile_stat_density_mean_95CI_plot.R H3K9me2_in_bodies 'H3K9me2' 'Chr1,Chr2,Chr3,Chr4,Chr5' wSNV 4 1.00 0.2 '%1.1f' '%3.1f' '%1.1f' 0.65
/applications/R/R-3.5.0/bin/Rscript CEN180_quantile_stat_density_mean_95CI_plot.R H3K27me1_in_bodies 'H3K27me1' 'Chr1,Chr2,Chr3,Chr4,Chr5' wSNV 4 1.00 0.2 '%1.1f' '%3.1f' '%1.2f' 0.65
/applications/R/R-3.5.0/bin/Rscript CEN180_quantile_stat_density_mean_95CI_plot.R HORlengthsSum 'Activity' 'Chr1,Chr2,Chr3,Chr4,Chr5' wSNV 4 1.00 0.2 '%1.0f' '%1.2f' '%1.0f' 0.65
/applications/R/R-3.5.0/bin/Rscript CEN180_quantile_stat_density_mean_95CI_plot.R array_size 'Array size' 'Chr1,Chr2,Chr3,Chr4,Chr5' wSNV 4 1.00 0.2 '%1.0f' '%1.2f' '%1.0f' 0.65

for i in $(seq 1 5)
do
  /applications/R/R-3.5.0/bin/Rscript CEN180_quantile_stat_density_mean_95CI_plot.R wSNV SNVs Chr${i} wSNV 4 1.00 0.2 '%2.0f' '%3.1f' '%2.0f' 0.65
  /applications/R/R-3.5.0/bin/Rscript CEN180_quantile_stat_density_mean_95CI_plot.R map_K150_E4_in_bodies 'Mappability (k=150 e=4)' Chr${i} wSNV 4 1.00 0.2 '%1.2f' '%3.1f' '%1.2f' 0.65
  /applications/R/R-3.5.0/bin/Rscript CEN180_quantile_stat_density_mean_95CI_plot.R CENH3_in_bodies 'CENH3' Chr${i} wSNV 4 1.00 0.2 '%1.1f' '%3.1f' '%1.1f' 0.65
  /applications/R/R-3.5.0/bin/Rscript CEN180_quantile_stat_density_mean_95CI_plot.R H3K9me2_in_bodies 'H3K9me2' Chr${i} wSNV 4 1.00 0.2 '%1.1f' '%3.1f' '%1.1f' 0.65
  /applications/R/R-3.5.0/bin/Rscript CEN180_quantile_stat_density_mean_95CI_plot.R H3K27me1_in_bodies 'H3K27me1' Chr${i} wSNV 4 1.00 0.2 '%1.1f' '%3.1f' '%1.2f' 0.65
  /applications/R/R-3.5.0/bin/Rscript CEN180_quantile_stat_density_mean_95CI_plot.R HORlengthsSum 'Activity' Chr${i} wSNV 4 1.00 0.2 '%1.0f' '%1.2f' '%1.0f' 0.65
  /applications/R/R-3.5.0/bin/Rscript CEN180_quantile_stat_density_mean_95CI_plot.R array_size 'Array size' Chr${i} wSNV 4 1.00 0.2 '%1.0f' '%1.2f' '%1.0f' 0.65
done


# HORlengthsSum 4 quantiles
/applications/R/R-3.5.0/bin/Rscript CEN180_quantile_stat_density_mean_95CI_plot.R wSNV SNVs 'Chr1,Chr2,Chr3,Chr4,Chr5' HORlengthsSum 4 1.00 0.2 '%2.0f' '%3.1f' '%2.0f' 0.65
/applications/R/R-3.5.0/bin/Rscript CEN180_quantile_stat_density_mean_95CI_plot.R map_K150_E4_in_bodies 'Mappability (k=150 e=4)' 'Chr1,Chr2,Chr3,Chr4,Chr5' HORlengthsSum 4 1.00 0.2 '%1.2f' '%3.1f' '%1.2f' 0.65
/applications/R/R-3.5.0/bin/Rscript CEN180_quantile_stat_density_mean_95CI_plot.R CENH3_in_bodies 'CENH3' 'Chr1,Chr2,Chr3,Chr4,Chr5' HORlengthsSum 4 1.00 0.2 '%1.1f' '%3.1f' '%1.1f' 0.65
/applications/R/R-3.5.0/bin/Rscript CEN180_quantile_stat_density_mean_95CI_plot.R H3K9me2_in_bodies 'H3K9me2' 'Chr1,Chr2,Chr3,Chr4,Chr5' HORlengthsSum 4 1.00 0.2 '%1.1f' '%3.1f' '%1.1f' 0.65
/applications/R/R-3.5.0/bin/Rscript CEN180_quantile_stat_density_mean_95CI_plot.R H3K27me1_in_bodies 'H3K27me1' 'Chr1,Chr2,Chr3,Chr4,Chr5' HORlengthsSum 4 1.00 0.2 '%1.1f' '%3.1f' '%1.2f' 0.65
/applications/R/R-3.5.0/bin/Rscript CEN180_quantile_stat_density_mean_95CI_plot.R HORlengthsSum 'Activity' 'Chr1,Chr2,Chr3,Chr4,Chr5' HORlengthsSum 4 1.00 0.2 '%1.0f' '%1.2f' '%1.0f' 0.65
/applications/R/R-3.5.0/bin/Rscript CEN180_quantile_stat_density_mean_95CI_plot.R array_size 'Array size' 'Chr1,Chr2,Chr3,Chr4,Chr5' HORlengthsSum 4 1.00 0.2 '%1.0f' '%1.2f' '%1.0f' 0.65

for i in $(seq 1 5)
do
  /applications/R/R-3.5.0/bin/Rscript CEN180_quantile_stat_density_mean_95CI_plot.R wSNV SNVs Chr${i} HORlengthsSum 4 1.00 0.2 '%2.0f' '%3.1f' '%2.0f' 0.65
  /applications/R/R-3.5.0/bin/Rscript CEN180_quantile_stat_density_mean_95CI_plot.R map_K150_E4_in_bodies 'Mappability (k=150 e=4)' Chr${i} HORlengthsSum 4 1.00 0.2 '%1.2f' '%3.1f' '%1.2f' 0.65
  /applications/R/R-3.5.0/bin/Rscript CEN180_quantile_stat_density_mean_95CI_plot.R CENH3_in_bodies 'CENH3' Chr${i} HORlengthsSum 4 1.00 0.2 '%1.1f' '%3.1f' '%1.1f' 0.65
  /applications/R/R-3.5.0/bin/Rscript CEN180_quantile_stat_density_mean_95CI_plot.R H3K9me2_in_bodies 'H3K9me2' Chr${i} HORlengthsSum 4 1.00 0.2 '%1.1f' '%3.1f' '%1.1f' 0.65
  /applications/R/R-3.5.0/bin/Rscript CEN180_quantile_stat_density_mean_95CI_plot.R H3K27me1_in_bodies 'H3K27me1' Chr${i} HORlengthsSum 4 1.00 0.2 '%1.1f' '%3.1f' '%1.2f' 0.65
  /applications/R/R-3.5.0/bin/Rscript CEN180_quantile_stat_density_mean_95CI_plot.R HORlengthsSum 'Activity' Chr${i} HORlengthsSum 4 1.00 0.2 '%1.0f' '%1.2f' '%1.0f' 0.65
  /applications/R/R-3.5.0/bin/Rscript CEN180_quantile_stat_density_mean_95CI_plot.R array_size 'Array size' Chr${i} HORlengthsSum 4 1.00 0.2 '%1.0f' '%1.2f' '%1.0f' 0.65
done


# CENH3_in_bodies 4 quantiles
/applications/R/R-3.5.0/bin/Rscript CEN180_quantile_stat_density_mean_95CI_plot.R wSNV SNVs 'Chr1,Chr2,Chr3,Chr4,Chr5' CENH3_in_bodies 4 1.00 0.2 '%2.0f' '%3.1f' '%2.0f' 0.65
/applications/R/R-3.5.0/bin/Rscript CEN180_quantile_stat_density_mean_95CI_plot.R map_K150_E4_in_bodies 'Mappability (k=150 e=4)' 'Chr1,Chr2,Chr3,Chr4,Chr5' CENH3_in_bodies 4 1.00 0.2 '%1.2f' '%3.1f' '%1.2f' 0.65
/applications/R/R-3.5.0/bin/Rscript CEN180_quantile_stat_density_mean_95CI_plot.R CENH3_in_bodies 'CENH3' 'Chr1,Chr2,Chr3,Chr4,Chr5' CENH3_in_bodies 4 1.00 0.2 '%1.1f' '%3.1f' '%1.1f' 0.65
/applications/R/R-3.5.0/bin/Rscript CEN180_quantile_stat_density_mean_95CI_plot.R H3K9me2_in_bodies 'H3K9me2' 'Chr1,Chr2,Chr3,Chr4,Chr5' CENH3_in_bodies 4 1.00 0.2 '%1.1f' '%3.1f' '%1.1f' 0.65
/applications/R/R-3.5.0/bin/Rscript CEN180_quantile_stat_density_mean_95CI_plot.R H3K27me1_in_bodies 'H3K27me1' 'Chr1,Chr2,Chr3,Chr4,Chr5' CENH3_in_bodies 4 1.00 0.2 '%1.1f' '%3.1f' '%1.2f' 0.65
/applications/R/R-3.5.0/bin/Rscript CEN180_quantile_stat_density_mean_95CI_plot.R HORlengthsSum 'Activity' 'Chr1,Chr2,Chr3,Chr4,Chr5' CENH3_in_bodies 4 1.00 0.2 '%1.0f' '%1.2f' '%1.0f' 0.65
/applications/R/R-3.5.0/bin/Rscript CEN180_quantile_stat_density_mean_95CI_plot.R array_size 'Array size' 'Chr1,Chr2,Chr3,Chr4,Chr5' CENH3_in_bodies 4 1.00 0.2 '%1.0f' '%1.2f' '%1.0f' 0.65

for i in $(seq 1 5)
do
  /applications/R/R-3.5.0/bin/Rscript CEN180_quantile_stat_density_mean_95CI_plot.R wSNV SNVs Chr${i} CENH3_in_bodies 4 1.00 0.2 '%2.0f' '%3.1f' '%2.0f' 0.65
  /applications/R/R-3.5.0/bin/Rscript CEN180_quantile_stat_density_mean_95CI_plot.R map_K150_E4_in_bodies 'Mappability (k=150 e=4)' Chr${i} CENH3_in_bodies 4 1.00 0.2 '%1.2f' '%3.1f' '%1.2f' 0.65
  /applications/R/R-3.5.0/bin/Rscript CEN180_quantile_stat_density_mean_95CI_plot.R CENH3_in_bodies 'CENH3' Chr${i} CENH3_in_bodies 4 1.00 0.2 '%1.1f' '%3.1f' '%1.1f' 0.65
  /applications/R/R-3.5.0/bin/Rscript CEN180_quantile_stat_density_mean_95CI_plot.R H3K9me2_in_bodies 'H3K9me2' Chr${i} CENH3_in_bodies 4 1.00 0.2 '%1.1f' '%3.1f' '%1.1f' 0.65
  /applications/R/R-3.5.0/bin/Rscript CEN180_quantile_stat_density_mean_95CI_plot.R H3K27me1_in_bodies 'H3K27me1' Chr${i} CENH3_in_bodies 4 1.00 0.2 '%1.1f' '%3.1f' '%1.2f' 0.65
  /applications/R/R-3.5.0/bin/Rscript CEN180_quantile_stat_density_mean_95CI_plot.R HORlengthsSum 'Activity' Chr${i} CENH3_in_bodies 4 1.00 0.2 '%1.0f' '%1.2f' '%1.0f' 0.65
  /applications/R/R-3.5.0/bin/Rscript CEN180_quantile_stat_density_mean_95CI_plot.R array_size 'Array size' Chr${i} CENH3_in_bodies 4 1.00 0.2 '%1.0f' '%1.2f' '%1.0f' 0.65
done


# map_K150_E4_in_bodies 4 quantiles
/applications/R/R-3.5.0/bin/Rscript CEN180_quantile_stat_density_mean_95CI_plot.R wSNV SNVs 'Chr1,Chr2,Chr3,Chr4,Chr5' map_K150_E4_in_bodies 4 1.00 0.2 '%2.0f' '%3.1f' '%2.0f' 0.65
/applications/R/R-3.5.0/bin/Rscript CEN180_quantile_stat_density_mean_95CI_plot.R map_K150_E4_in_bodies 'Mappability (k=150 e=4)' 'Chr1,Chr2,Chr3,Chr4,Chr5' map_K150_E4_in_bodies 4 1.00 0.2 '%1.2f' '%3.1f' '%1.2f' 0.65
/applications/R/R-3.5.0/bin/Rscript CEN180_quantile_stat_density_mean_95CI_plot.R CENH3_in_bodies 'CENH3' 'Chr1,Chr2,Chr3,Chr4,Chr5' map_K150_E4_in_bodies 4 1.00 0.2 '%1.1f' '%3.1f' '%1.1f' 0.65
/applications/R/R-3.5.0/bin/Rscript CEN180_quantile_stat_density_mean_95CI_plot.R H3K9me2_in_bodies 'H3K9me2' 'Chr1,Chr2,Chr3,Chr4,Chr5' map_K150_E4_in_bodies 4 1.00 0.2 '%1.1f' '%3.1f' '%1.1f' 0.65
/applications/R/R-3.5.0/bin/Rscript CEN180_quantile_stat_density_mean_95CI_plot.R H3K27me1_in_bodies 'H3K27me1' 'Chr1,Chr2,Chr3,Chr4,Chr5' map_K150_E4_in_bodies 4 1.00 0.2 '%1.1f' '%3.1f' '%1.2f' 0.65
/applications/R/R-3.5.0/bin/Rscript CEN180_quantile_stat_density_mean_95CI_plot.R HORlengthsSum 'Activity' 'Chr1,Chr2,Chr3,Chr4,Chr5' map_K150_E4_in_bodies 4 1.00 0.2 '%1.0f' '%1.2f' '%1.0f' 0.65
/applications/R/R-3.5.0/bin/Rscript CEN180_quantile_stat_density_mean_95CI_plot.R array_size 'Array size' 'Chr1,Chr2,Chr3,Chr4,Chr5' map_K150_E4_in_bodies 4 1.00 0.2 '%1.0f' '%1.2f' '%1.0f' 0.65

for i in $(seq 1 5)
do
  /applications/R/R-3.5.0/bin/Rscript CEN180_quantile_stat_density_mean_95CI_plot.R wSNV SNVs Chr${i} map_K150_E4_in_bodies 4 1.00 0.2 '%2.0f' '%3.1f' '%2.0f' 0.65
  /applications/R/R-3.5.0/bin/Rscript CEN180_quantile_stat_density_mean_95CI_plot.R map_K150_E4_in_bodies 'Mappability (k=150 e=4)' Chr${i} map_K150_E4_in_bodies 4 1.00 0.2 '%1.2f' '%3.1f' '%1.2f' 0.65
  /applications/R/R-3.5.0/bin/Rscript CEN180_quantile_stat_density_mean_95CI_plot.R CENH3_in_bodies 'CENH3' Chr${i} map_K150_E4_in_bodies 4 1.00 0.2 '%1.1f' '%3.1f' '%1.1f' 0.65
  /applications/R/R-3.5.0/bin/Rscript CEN180_quantile_stat_density_mean_95CI_plot.R H3K9me2_in_bodies 'H3K9me2' Chr${i} map_K150_E4_in_bodies 4 1.00 0.2 '%1.1f' '%3.1f' '%1.1f' 0.65
  /applications/R/R-3.5.0/bin/Rscript CEN180_quantile_stat_density_mean_95CI_plot.R H3K27me1_in_bodies 'H3K27me1' Chr${i} map_K150_E4_in_bodies 4 1.00 0.2 '%1.1f' '%3.1f' '%1.2f' 0.65
  /applications/R/R-3.5.0/bin/Rscript CEN180_quantile_stat_density_mean_95CI_plot.R HORlengthsSum 'Activity' Chr${i} map_K150_E4_in_bodies 4 1.00 0.2 '%1.0f' '%1.2f' '%1.0f' 0.65
  /applications/R/R-3.5.0/bin/Rscript CEN180_quantile_stat_density_mean_95CI_plot.R array_size 'Array size' Chr${i} map_K150_E4_in_bodies 4 1.00 0.2 '%1.0f' '%1.2f' '%1.0f' 0.65
done
