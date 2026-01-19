prepare_growth_panel <- function(raw_df,
                                 countries = growth_countries,
                                 start_year = start_year_growth) {
  df <- raw_df |>
    dplyr::mutate(
      year = year_to_date(year),
      country = as.character(country)
    ) |>
    dplyr::filter(
      country %in% countries,
      year >= as.Date(paste0(start_year, "-01-01"))
    ) |>
    dplyr::select(country, year, pop, gdppc)
  list(
    panel = df,
    spain = dplyr::filter(df, country == "Spain")
  )
}

augment_spain_growth <- function(spain_df,
                                 windows = regime_windows,
                                 levels = regime_levels,
                                 ma_window = moving_average_window) {
  spain_df |>
    dplyr::mutate(
      growth = as.numeric(safe_delt(gdppc, k = 1)),
      mm = as.numeric(moving_average(growth, ma_window)), # trailing MA (sin look-ahead)
      loess_trend = dplyr::coalesce(
        one_sided_loess(year, growth),
        mm,
        growth
      ),
      regimes = apply_regime_labels(year, windows, levels)
    )
}

build_regime_windows <- function(windows = regime_windows) {
  windows$regimes <- factor(windows$regimes, levels = regime_levels)
  windows
}

build_western_europe_share <- function(panel_df, spain_df) {
  income <- panel_df |>
    dplyr::group_by(year) |>
    dplyr::summarise(
      Western_Europe = stats::weighted.mean(
        gdppc[country != "Spain"],
        pop[country != "Spain"],
        na.rm = TRUE
      ),
      .groups = "drop"
    )

  res <- dplyr::tibble(
    year = spain_df$year,
    Spain = spain_df$gdppc
  ) |>
    dplyr::left_join(income, by = "year")

  missing_years <- res$year[is.na(res$Western_Europe)]
  if (length(missing_years) > 0) {
    years_str <- paste(sort(unique(format(missing_years, "%Y"))), collapse = ", ")
    stop(sprintf("Western Europe income missing for years: %s", years_str))
  }

  res |>
    dplyr::mutate(Share = Spain / Western_Europe)
}

prepare_pwt11_panel <- function(raw_df,
                                countries = synthetic_countries,
                                code_map = region_codes,
                                max_year = NULL) {
  df <- raw_df

  if (!is.null(countries)) {
    df <- dplyr::filter(df, country %in% countries)
  }

  df <- df |>
    dplyr::mutate(
      country = as.character(country),
      country = textclean::mgsub(
        country,
        c("Bolivia (Plurinational State of)", "Venezuela (Bolivarian Republic of)"),
        c("Bolivia", "Venezuela")
      ),
      gdpcap = rgdpna / pop,
      rknacapita = rnna / pop
    )

  if (!is.null(max_year)) {
    df <- dplyr::filter(df, year <= max_year)
  }

  vars_keep <- c(
    "country", "isocode", "year", "pop", "hc", "gdpcap",
    "csh_c", "csh_i", "csh_g", "csh_x", "csh_m",
    "rknacapita", "pl_gdpo"
  )

  df <- dplyr::select(df, dplyr::any_of(vars_keep))

  avg <- df |>
    dplyr::group_by(year) |>
    dplyr::summarise(
      dplyr::across(where(is.numeric), \(x) mean(x, na.rm = TRUE)),
      .groups = "drop"
    ) |>
    dplyr::mutate(country = "Average", isocode = "AVG")

  df <- dplyr::bind_rows(avg, df) |>
    dplyr::mutate(country = forcats::fct_relevel(country, "Average")) |>
    dplyr::arrange(country, year)

  if (is.null(code_map)) {
    codes <- sort(unique(df$isocode))
    code_map <- seq_along(codes)
    names(code_map) <- codes
  }

  df$region <- as.integer(unname(code_map[df$isocode]))
  df
}

