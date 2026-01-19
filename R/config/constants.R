MADDISON_URL <- "https://dataverse.nl/api/access/datafile/421303" # MPD 2023 Stata
# PWT 11.0 hosted in dataverse.nl (public file id 554105)
PWT11_URL <- "https://dataverse.nl/api/access/datafile/554105"
PWT11_LOCAL <- file.path("data", "raw", "pwt110.xlsx")
MADDISON_LOCAL <- file.path("data", "raw", "maddison_mpd2023.dta")

PWT11_CHECKSUM_KEY <- "pwt110"
MADDISON_CHECKSUM_KEY <- "maddison_mpd2023"
CHECKSUM_FILE <- file.path("data", "raw", "checksums.txt")

# Output folders
OUTPUT_BASE   <- "output"
OUTPUT_SESSIONS <- file.path(OUTPUT_BASE, "sessions")
OUTPUT_PLOTS  <- "plots"
OUTPUT_TABLES <- "tables"
