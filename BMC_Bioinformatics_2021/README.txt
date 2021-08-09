GSE14468.Rmd is the code that was used to process the raw CEL files and combine the CEL files with the phenotypic data that should be downloaded from Gene Expression Omnibus from GSE14468_series_matrix.txt.gz using the Bioconductor GEOqery package. 
To replicate this code, download the CEL files from Gene Expression Omnibus and stored in the baseDir in a folder named GSE14468_RAW. A "GSE14468.RData" object is created from running GSE14468.Rmd.
Required packages include digest, purrr, stringr, affyio, lubridate, GEOquery, affy.


Code for fitting models:
Note 1: change your working directory by specifying your directory that houses the GSE14468.RData object in the line: setwd("your-directory-that-stores-GSE14468.RData").
Note 2: The code assumes the user has access to a machine that can leverage three nodes for running the three chains using parallel processing.
Note 3: Requires packages include affy, HDInterval, MASS, Matrix, clone  

MODEL I
GSE14468_Bayesian_Model1.Rmd

MODEL II
Beta(1,19) prior      GSE14468_Bayesian_Model2_Beta_1_19.Rmd
Beta(0.01,0.19) prior GSE14468_Bayesian_Model2_Beta_1_19.Rmd
Fixed pi=0.05 prior   GSE14468_Bayesian_Model2_Fixed_05.Rmd
Fixed pi=0.50 prior   GSE14468_Bayesian_Model2_Fixed_50.Rmd

MODEL III
Beta(1,19) prior      GSE14468_Bayesian_Model3_Beta_1_19.Rmd
Beta(0.01,0.19) prior GSE14468_Bayesian_Model3_Beta_1_19.Rmd
Fixed pi=0.05 prior   GSE14468_Bayesian_Model3_Fixed_05.Rmd
Fixed pi=0.50 prior   GSE14468_Bayesian_Model3_Fixed_50.Rmd

MODEL IV
Beta(1,19) prior      GSE14468_Bayesian_Model4_Beta_1_19.Rmd
Beta(0.01,0.19) prior GSE14468_Bayesian_Model4_Beta_1_19.Rmd
Fixed pi=0.05 prior   GSE14468_Bayesian_Model4_Fixed_05.Rmd
Fixed pi=0.50 prior   GSE14468_Bayesian_Model4_Fixed_50.Rmd

Summarizing results from all 12 models appears in GSE14468_Summarize_Results.R

NOTE: The following version of software was used
R/RStudio 4.1.0;  dclone 2.3-0
version 4.1.0 (2021-05-18)
Platform: x86_64-apple-darwin17.0 (64-bit)
Running under: macOS Big Sur 11.2.3

Matrix products: default
LAPACK: /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRlapack.dylib

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] parallel  stats     graphics  grDevices utils     datasets
[7] methods   base

other attached packages:
[1] dclone_2.3-0        coda_0.19-4         Matrix_1.3-4
[4] MASS_7.3-54         HDInterval_0.2.2    affy_1.70.0
[7] Biobase_2.52.0      BiocGenerics_0.38.0

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.6            lattice_0.20-44       prettyunits_1.1.1
 [4] ps_1.6.0              digest_0.6.27         utf8_1.2.1
 [7] V8_3.4.2              R6_2.5.0              stats4_4.1.0
[10] evaluate_0.14         ggplot2_3.3.5         pillar_1.6.1
[13] zlibbioc_1.38.0       rlang_0.4.11          curl_4.3.2
[16] callr_3.7.0           rjags_4-10            preprocessCore_1.54.0
[19] rmarkdown_2.9         loo_2.4.1             munsell_0.5.0
[22] compiler_4.1.0        xfun_0.24             rstan_2.21.2
[25] pkgconfig_2.0.3       pkgbuild_1.2.0        htmltools_0.5.1.1
[28] tidyselect_1.1.1      tibble_3.1.2          gridExtra_2.3
[31] R2WinBUGS_2.1-21      codetools_0.2-18      matrixStats_0.59.0
[34] fansi_0.5.0           crayon_1.4.1          dplyr_1.0.7
[37] withr_2.4.2           grid_4.1.0            jsonlite_1.7.2
[40] gtable_0.3.0          lifecycle_1.0.0       DBI_1.1.1
[43] magrittr_2.0.1        StanHeaders_2.21.0-7  scales_1.1.1
[46] RcppParallel_5.1.4    cli_3.0.0             affyio_1.62.0
[49] ellipsis_0.3.2        vctrs_0.3.8           generics_0.1.0
[52] boot_1.3-28           tools_4.1.0           glue_1.4.2
[55] purrr_0.3.4           processx_3.5.2        yaml_2.2.1
[58] inline_0.3.19         colorspace_2.0-2      BiocManager_1.30.16
[61] knitr_1.33
