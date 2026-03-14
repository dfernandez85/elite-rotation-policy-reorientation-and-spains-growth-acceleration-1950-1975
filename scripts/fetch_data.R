# Utility script to download the pinned PWT input and write its md5 checksum.
# Usage: Rscript scripts/fetch_data.R

source("R/config/constants.R")
source("R/utils/helpers.R")
source("R/data/downloaders.R")

fetch_pinned_data(force = TRUE)
