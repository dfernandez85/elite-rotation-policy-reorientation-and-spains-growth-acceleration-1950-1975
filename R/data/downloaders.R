load_maddison_data <- function(local_path = MADDISON_LOCAL,
                               checksum_key = MADDISON_CHECKSUM_KEY) {
  validate_file_hash(local_path, checksum_key)
  rio::import(local_path)
}

load_pwt_data <- function(local_path = PWT11_LOCAL,
                          checksum_key = PWT11_CHECKSUM_KEY) {
  validate_file_hash(local_path, checksum_key)

  read_pwt11 <- function(path) {
    readxl::read_excel(path, sheet = "Data") |>
      dplyr::rename(isocode = countrycode)
  }

  tibble::as_tibble(read_pwt11(local_path))
}
