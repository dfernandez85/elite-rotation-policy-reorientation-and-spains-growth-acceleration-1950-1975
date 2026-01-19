source("R/config/constants.R")
source("R/config/packages.R")
source("R/config/settings.R")
source("R/utils/helpers.R")
source("R/data/downloaders.R")
source("R/data/processors.R")
source("R/visualization/plots.R")
load_required_packages()

raw <- load_pwt_data()
panel <- prepare_pwt11_panel(raw)
panel$region <- as.numeric(panel$region)
panel_df <- as.data.frame(panel)
panel_df$country <- as.character(panel_df$country)

coverage_df <- panel_df |>
  dplyr::filter(country != "Average") |>
  dplyr::group_by(year) |>
  dplyr::summarise(
    n_countries = dplyr::n(),
    dplyr::across(where(is.numeric), \(x) sum(is.na(x)), .names = "na_{.col}"),
    dplyr::across(where(is.numeric), \(x) mean(is.na(x)), .names = "na_share_{.col}"),
    .groups = "drop"
  )

key_vars <- c("pop", "gdpcap", "hc", "csh_c", "csh_i", "csh_g", "csh_x", "csh_m", "rknacapita", "pl_gdpo")
pre_window <- range(pre_treatment_years)
pre_years <- pre_window[1]:pre_window[2]
pre_len <- length(pre_years)
coverage_pre <- panel_df |>
  dplyr::filter(
    country != "Average",
    year >= pre_window[1], year <= pre_window[2]
  ) |>
  dplyr::group_by(country) |>
  dplyr::summarise(
    dplyr::across(dplyr::all_of(key_vars), \(x) sum(!is.na(x)), .names = "n_{.col}"),
    .groups = "drop"
  )

excluded_pre <- coverage_pre |>
  tidyr::pivot_longer(cols = dplyr::starts_with("n_"), names_to = "var", values_to = "non_na") |>
  dplyr::mutate(var = gsub("^n_", "", var)) |>
  dplyr::filter(non_na < pre_len) |>
  dplyr::mutate(reason = paste0("Insufficient pre data for ", var, " (", non_na, "/", pre_len, ")")) |>
  dplyr::select(country, reason) |>
  dplyr::distinct()

post_window <- c(pre_window[2] + 1, post_cutoff_year - 1)
excluded_post <- dplyr::tibble()
coverage_post <- NULL
if (post_window[1] <= post_window[2]) {
  post_len <- post_window[2] - post_window[1] + 1
  coverage_post <- panel_df |>
    dplyr::filter(
      country != "Average",
      year >= post_window[1], year <= post_window[2]
    ) |>
    dplyr::group_by(country) |>
    dplyr::summarise(
      dplyr::across(dplyr::all_of(key_vars), \(x) sum(!is.na(x)), .names = "n_{.col}"),
      .groups = "drop"
    )

  excluded_post <- coverage_post |>
    tidyr::pivot_longer(cols = dplyr::starts_with("n_"), names_to = "var", values_to = "non_na") |>
    dplyr::mutate(var = gsub("^n_", "", var)) |>
    dplyr::filter(non_na < post_len) |>
    dplyr::mutate(reason = paste0("Insufficient post data for ", var, " (", non_na, "/", post_len, ")")) |>
    dplyr::select(country, reason) |>
    dplyr::distinct()
}

excluded <- dplyr::bind_rows(excluded_pre, excluded_post) |>
  dplyr::group_by(country) |>
  dplyr::summarise(reason = paste(unique(reason), collapse = "; "), .groups = "drop")

if (nrow(excluded) > 0) {
  panel_df <- dplyr::filter(panel_df, !country %in% excluded$country)
}

panel_no_avg <- panel_df |>
  dplyr::filter(country != "Average")
if (nrow(panel_no_avg) == 0) {
  stop("No donor countries remain after coverage filters (pre/post).")
}
region_avg <- unique(panel$region[panel$isocode == "AVG"])[1]
if (length(region_avg) == 0 || is.na(region_avg)) {
  region_avg <- max(panel_no_avg$region, na.rm = TRUE) + 1
}
if (!is.finite(region_avg)) {
  region_avg <- 1
}
numeric_cols <- setdiff(names(dplyr::select(panel_no_avg, where(is.numeric))), c("region", "year"))
avg_clean <- panel_no_avg |>
  dplyr::group_by(year) |>
  dplyr::summarise(
    dplyr::across(dplyr::all_of(numeric_cols), \(x) mean(x, na.rm = TRUE)),
    .groups = "drop"
  ) |>
  dplyr::mutate(country = "Average", isocode = "AVG", region = region_avg)
panel_df <- dplyr::bind_rows(avg_clean, panel_no_avg) |>
  dplyr::mutate(country = forcats::fct_relevel(country, "Average")) |>
  dplyr::arrange(country, year)
panel_df$region <- as.integer(panel_df$region)

print(sapply(panel_df, class))
print(head(panel_df$region))
print(table(panel_df$country)[1:5])

# Try listFromLong
foo <- MSCMT::listFromLong(
  panel_df,
  unit.variable = "region",
  time.variable = "year",
  unit.names.variable = "country"
)
print(str(foo))
