analyze_growth_regimes <- function(session_dir) {
  plots_dir <- file.path(session_dir, OUTPUT_PLOTS)
  tables_dir <- file.path(session_dir, OUTPUT_TABLES)
  ensure_dir(plots_dir)
  ensure_dir(tables_dir)

  raw <- load_maddison_data()
  prepared <- prepare_growth_panel(raw)
  windows_df <- build_regime_windows()

  spain <- augment_spain_growth(prepared$spain, windows = regime_windows, levels = regime_levels)
  share_df <- build_western_europe_share(prepared$panel, prepared$spain)

  regime_growth <- spain |>
    dplyr::filter(!is.na(regimes), regimes != "Unclassified", !is.na(gdppc)) |>
    dplyr::arrange(regimes, year) |>
    dplyr::group_by(regimes) |>
    dplyr::summarise(
      start_year_date = dplyr::first(year),
      end_year_date = dplyr::last(year),
      start_rgdpnapc = dplyr::first(gdppc),
      end_rgdpnapc = dplyr::last(gdppc),
      start_year = as.integer(format(start_year_date, "%Y")),
      end_year = as.integer(format(end_year_date, "%Y")),
      duration_years = end_year - start_year,
      CAGR = dplyr::if_else(
        duration_years > 0 & start_rgdpnapc > 0,
        (end_rgdpnapc / start_rgdpnapc)^(1 / duration_years) - 1,
        NA_real_
      ),
      .groups = "drop"
    ) |>
    dplyr::select(regimes, start_year, end_year, duration_years, CAGR)

  p1 <- plot_growth_trend(spain, windows_df)
  civil_war_highlight <- data.frame(
    start = as.Date("1936-01-01"),
    end = as.Date("1939-12-31")
  )

  p1_hp <- plot_growth_trend(
    spain,
    windows_df,
    value_col = "loess_trend",
    title = "Figure 1b Economic Growth (LOESS one-sided) in Spain during Different Regimes",
    y_label = "Growth in real GDP per capita (LOESS one-sided)",
    y_limits = c(-0.1, 0.1),
    y_breaks = seq(-0.1, 0.1, 0.02),
    highlight_ranges = civil_war_highlight
  )
  p2 <- plot_income_share(share_df, windows_df)

  ggplot2::ggsave(file.path(plots_dir, "Figure_1_growth.png"), p1, width = 20, height = 11, units = "cm", dpi = 600)
  ggplot2::ggsave(file.path(plots_dir, "Figure_1b_growth_loess.png"), p1_hp, width = 20, height = 11, units = "cm", dpi = 600)
  ggplot2::ggsave(file.path(plots_dir, "Figure_2_income_share.png"), p2, width = 20, height = 11, units = "cm", dpi = 600)

  readr::write_csv(spain, file.path(tables_dir, "spain_growth_series.csv"))
  readr::write_csv(share_df, file.path(tables_dir, "income_share_western_europe.csv"))
  readr::write_csv(regime_growth, file.path(tables_dir, "growth_by_regime.csv"))
  save_table_jpg(spain, file.path(tables_dir, "spain_growth_series.jpg"), title = "Spain growth series (truncated)")
  save_table_jpg(share_df, file.path(tables_dir, "income_share_western_europe.jpg"), title = "Income share vs Western Europe (truncated)")
  save_table_jpg(regime_growth, file.path(tables_dir, "growth_by_regime.jpg"), title = "Average GDP pc growth by regime")

  list(
    growth_plot = p1,
    growth_plot_hp = p1_hp,
    income_share_plot = p2,
    growth_data = spain,
    share_data = share_df,
    growth_by_regime = regime_growth
  )
}
