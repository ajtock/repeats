#!/bin/bash

# wSNV 4 quantiles
/applications/R/R-4.0.0/bin/Rscript CEN180_quantile_frequency_chrProfiles.R 'Chr1,Chr2,Chr3,Chr4,Chr5' wSNV 4 10000 101
/applications/R/R-4.0.0/bin/Rscript CEN180_quantile_frequency_chrProfilesPlot.R 'Chr1,Chr2,Chr3,Chr4,Chr5' wSNV 4 10000

for i in $(seq 1 5)
do
  /applications/R/R-4.0.0/bin/Rscript CEN180_quantile_frequency_chrProfiles.R Chr${i} wSNV 4 10000 101
  /applications/R/R-4.0.0/bin/Rscript CEN180_quantile_frequency_chrProfilesPlot.R Chr${i} wSNV 4 10000
done

# HORlengthsSum 4 quantiles
/applications/R/R-4.0.0/bin/Rscript CEN180_quantile_frequency_chrProfiles.R 'Chr1,Chr2,Chr3,Chr4,Chr5' HORlengthsSum 4 10000 101
/applications/R/R-4.0.0/bin/Rscript CEN180_quantile_frequency_chrProfilesPlot.R 'Chr1,Chr2,Chr3,Chr4,Chr5' HORlengthsSum 4 10000

for i in $(seq 1 5)
do
  /applications/R/R-4.0.0/bin/Rscript CEN180_quantile_frequency_chrProfiles.R Chr${i} HORlengthsSum 4 10000 101
  /applications/R/R-4.0.0/bin/Rscript CEN180_quantile_frequency_chrProfilesPlot.R Chr${i} HORlengthsSum 4 10000
done

# CENH3_in_bodies 4 quantiles
/applications/R/R-4.0.0/bin/Rscript CEN180_quantile_frequency_chrProfiles.R 'Chr1,Chr2,Chr3,Chr4,Chr5' CENH3_in_bodies 4 10000 101
/applications/R/R-4.0.0/bin/Rscript CEN180_quantile_frequency_chrProfilesPlot.R 'Chr1,Chr2,Chr3,Chr4,Chr5' CENH3_in_bodies 4 10000

for i in $(seq 1 5)
do
  /applications/R/R-4.0.0/bin/Rscript CEN180_quantile_frequency_chrProfiles.R Chr${i} CENH3_in_bodies 4 10000 101
  /applications/R/R-4.0.0/bin/Rscript CEN180_quantile_frequency_chrProfilesPlot.R Chr${i} CENH3_in_bodies 4 10000
done

# map_K150_E4_in_bodies 4 quantiles
/applications/R/R-4.0.0/bin/Rscript CEN180_quantile_frequency_chrProfiles.R 'Chr1,Chr2,Chr3,Chr4,Chr5' map_K150_E4_in_bodies 4 10000 101
/applications/R/R-4.0.0/bin/Rscript CEN180_quantile_frequency_chrProfilesPlot.R 'Chr1,Chr2,Chr3,Chr4,Chr5' map_K150_E4_in_bodies 4 10000

for i in $(seq 1 5)
do
  /applications/R/R-4.0.0/bin/Rscript CEN180_quantile_frequency_chrProfiles.R Chr${i} map_K150_E4_in_bodies 4 10000 101
  /applications/R/R-4.0.0/bin/Rscript CEN180_quantile_frequency_chrProfilesPlot.R Chr${i} map_K150_E4_in_bodies 4 10000
done
