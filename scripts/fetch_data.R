# Utility script to download pinned data files and write their md5 checksums.
# Usage: Rscript scripts/fetch_data.R

source("R/config/constants.R")
source("R/utils/helpers.R")
source("R/data/downloaders.R")

fetch_pinned_data(force = TRUE)
