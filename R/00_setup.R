# ============================================================
# 00_setup.R
# Install and load all packages needed for the replication.
# ============================================================

required_packages <- c(
  "tidyverse",   # dplyr, ggplot2, tidyr, readr, purrr, etc.
  "plm",         # panel linear models (within, fd, pooled)
  "sandwich",    # cluster-robust / HC standard errors
  "lmtest",      # coeftest() with robust SEs
  "broom",       # tidy model summaries
  "knitr",       # table rendering
  "kableExtra",  # enhanced kable tables
  "haven",       # read Stata .dta files
  "here",        # portable file paths
  "ggplot2",     # graphics (already via tidyverse, listed for clarity)
  "scales",      # axis formatting in ggplot2
  "patchwork"    # combine ggplot2 panels
)

# Install any packages that are not yet available
new_pkgs <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]
if (length(new_pkgs) > 0) {
  message("Installing missing packages: ", paste(new_pkgs, collapse = ", "))
  install.packages(new_pkgs, repos = "https://cloud.r-project.org")
}

# Load all packages (suppress startup messages for cleanliness)
invisible(lapply(required_packages, function(pkg) {
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}))

message("All required packages loaded successfully.")
