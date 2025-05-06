# Installing essential CRAN packages
install.packages(c(
  "stringr",    # Advanced string manipulation
  "readxl",     # Reading Excel files (.xlsx)
  "openxlsx",   # Creating and modifying Excel files
  "dplyr",      # Data manipulation
  "tidyr"       # Data cleaning and transformation
))

# Install ggplot2 at specific version
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_version("ggplot2", version = "3.3.4", repos = "https://cran.r-project.org")

# Installing Bioconductor packages required for SpliceLauncher
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

BiocManager::install(c(
  "WriteXLS",   # Exporting data in Excel format
  "Cairo"       # Advanced graphics rendering
  "Rsamtools",
  "GenomicRanges",
  "IRanges",
  "S4Vectors",
  "SummarizedExperiment",
  "BiocGenerics"
))