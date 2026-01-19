plot_growth_trend <- function(spain_df,
                              windows_df,
                              value_col = "mm",
                              title = "Figure 1 Economic Growth in Spain during Different Regimes",
                              y_label = "Growth in real GDP per capita (10-year moving average, trailing)",
                              y_limits = NULL,
                              y_breaks = NULL,
                              highlight_ranges = NULL) {
  line_df <- spain_df |>
    dplyr::mutate(value = .data[[value_col]])
  x_min <- as.Date(paste0(plot_start_year, "-01-01"))
  x_max <- max(c(line_df$year, windows_df$end), na.rm = TRUE)

  p <- ggplot2::ggplot() +
    ggplot2::geom_rect(
      data = windows_df, alpha = 0.9,
      ggplot2::aes(
        xmin = start, xmax = end,
        ymin = -Inf, ymax = Inf,
        fill = regimes
      )
    )

  if (!is.null(highlight_ranges) && nrow(highlight_ranges) > 0) {
    p <- p +
      ggplot2::geom_rect(
        data = highlight_ranges,
        ggplot2::aes(
          xmin = start, xmax = end,
          ymin = -Inf, ymax = Inf
        ),
        inherit.aes = FALSE,
        fill = "black",
        alpha = 0.18
      )
  }

  p <- p +
    ggplot2::geom_line(data = line_df, ggplot2::aes(year, value), linewidth = 1.2, na.rm = TRUE) +
    ggplot2::labs(
      title = title,
      y = y_label,
      x = "Year"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.ticks = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 0.5, vjust = 0.5),
      axis.text = ggplot2::element_text(face = "bold"),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      legend.title = ggplot2::element_blank()
    ) +
    ggplot2::scale_y_continuous(
      limits = y_limits,
      breaks = if (is.null(y_breaks)) scales::breaks_width(0.01) else y_breaks,
      minor_breaks = NULL,
      labels = scales::percent
    ) +
    ggplot2::scale_x_date(
      expand = ggplot2::expansion(mult = c(0, 0.01)),
      date_breaks = "5 years",
      minor_breaks = NULL,
      date_labels = "%Y"
    ) +
    ggplot2::coord_cartesian(xlim = c(x_min, x_max))

  p +
    ggplot2::geom_hline(
      yintercept = 0,
      linetype = "dashed",
      colour = "black"
    )
}

plot_income_share <- function(share_df,
                              windows_df,
                              value_col = "Share",
                              title = "Figure 2 GDP per Capita (Spain / Rest of Western Europe)",
                              y_label = "Spain income as % of the rest of Western Europe") {
  line_df <- share_df |>
    dplyr::mutate(value = .data[[value_col]])
  x_min <- as.Date(paste0(plot_start_year, "-01-01"))
  x_max <- max(c(line_df$year, windows_df$end), na.rm = TRUE)

  ggplot2::ggplot() +
    ggplot2::geom_rect(
      data = windows_df, alpha = 0.9,
      ggplot2::aes(
        xmin = start, xmax = end,
        ymin = -Inf, ymax = Inf,
        fill = regimes
      )
    ) +
    ggplot2::geom_line(data = line_df, ggplot2::aes(year, value), linewidth = 1.2, na.rm = TRUE) +
    ggplot2::labs(
      title = title,
      y = y_label,
      x = "Year"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.ticks = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 0.5, vjust = 0.5),
      axis.text = ggplot2::element_text(face = "bold"),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      legend.title = ggplot2::element_blank()
    ) +
    ggplot2::scale_y_continuous(
      breaks = seq(0.0, 1, 0.1),
      minor_breaks = NULL,
      labels = scales::percent
    ) +
    ggplot2::scale_x_date(
      expand = ggplot2::expansion(mult = c(0, 0.01)),
      date_breaks = "5 years",
      minor_breaks = NULL,
      date_labels = "%Y"
    ) +
    ggplot2::coord_cartesian(xlim = c(x_min, x_max))
}

plot_synth_comparison <- function(res,
                                  x_limits = c(1950, 1975),
                                  y_limits = NULL,
                                  y_breaks = waiver(),
                                  y_label = "Real GDP per capita",
                                  y_formatter = scales::dollar_format()) {
  xlim_dates <- NULL
  if (!is.null(x_limits)) {
    xlim_dates <- as.Date(paste0(x_limits, "-01-01"))
  }

  p <- suppressWarnings(
    ggplot2::ggplot(res,
      type = "comparison",
      ylab = y_label,
      xlab = "Year",
      main = "Figure 3 GDP per Capita Spain and Synthetic Spain",
      labels = c("Spain", "Synthetic Spain")
    )
  )

  p +
    ggplot2::theme(
      axis.ticks = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 0.5, vjust = 0.5),
      axis.text = ggplot2::element_text(face = "bold")
    ) +
    ggplot2::scale_x_date(
      expand = ggplot2::expansion(mult = c(0, 0.01)),
      date_breaks = "5 years",
      minor_breaks = NULL,
      date_labels = "%Y"
    ) +
    ggplot2::coord_cartesian(xlim = xlim_dates) +
    ggplot2::scale_y_continuous(
      limits = y_limits,
      breaks = y_breaks,
      minor_breaks = waiver(),
      labels = y_formatter,
      expand = expansion(mult = c(0, 0.02))
    )
}

plot_synth_gaps <- function(res,
                            band_df = NULL,
                            x_limits = c(1950, 1975),
                            y_limits = NULL,
                            y_breaks = waiver(),
                            y_label = "Differential GDP per capita",
                            y_formatter = scales::dollar_format(),
                            treatment_year = NULL) {
  xlim_dates <- NULL
  if (!is.null(x_limits)) {
    xlim_dates <- as.Date(paste0(x_limits, "-01-01"))
  }

  p <- suppressWarnings(
    ggplot2::ggplot(res,
      type = "gaps",
      ylab = y_label,
      xlab = "Year",
      main = "Figure 4 GDP per Capita Spain over Synthetic Spain",
      labels = c("Spain", "Synthetic Spain")
    )
  ) +
    ggplot2::theme(
      axis.ticks = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 0.5, vjust = 0.5),
      axis.text = ggplot2::element_text(face = "bold")
    ) +
    ggplot2::scale_x_date(
      expand = ggplot2::expansion(mult = c(0, 0.01)),
      date_breaks = "5 years",
      minor_breaks = NULL,
      date_labels = "%Y"
    ) +
    ggplot2::coord_cartesian(xlim = xlim_dates) +
    ggplot2::scale_y_continuous(
      limits = y_limits,
      breaks = y_breaks,
      minor_breaks = waiver(),
      labels = y_formatter,
      expand = expansion(mult = c(0, 0.02))
    )

  if (!is.null(treatment_year)) {
    treat_date <- as.Date(paste0(treatment_year, "-01-01"))
    p <- p + ggplot2::geom_vline(xintercept = treat_date, linetype = "dashed", colour = "darkred", linewidth = 0.6)
  }

  if (!is.null(band_df)) {
    band_df <- band_df |>
      dplyr::mutate(year = as.Date(paste0(year, "-01-01")))

    lower_col <- dplyr::case_when(
      "q05" %in% names(band_df) ~ "q05",
      "lower" %in% names(band_df) ~ "lower",
      TRUE ~ NA_character_
    )
    upper_col <- dplyr::case_when(
      "q95" %in% names(band_df) ~ "q95",
      "upper" %in% names(band_df) ~ "upper",
      TRUE ~ NA_character_
    )
    median_col <- dplyr::case_when(
      "q50" %in% names(band_df) ~ "q50",
      "median" %in% names(band_df) ~ "median",
      "treated_gap" %in% names(band_df) ~ "treated_gap",
      TRUE ~ NA_character_
    )

    if (!is.na(lower_col) && !is.na(upper_col)) {
      band_df$lower <- band_df[[lower_col]]
      band_df$upper <- band_df[[upper_col]]
    }
    if (!is.na(median_col)) {
      band_df$median <- band_df[[median_col]]
    }

    p <- p +
      ggplot2::geom_ribbon(
        data = band_df,
        ggplot2::aes(x = year, ymin = lower, ymax = upper),
        fill = "grey70",
        alpha = 0.3,
        inherit.aes = FALSE
      ) +
      ggplot2::geom_line(
        data = band_df,
        ggplot2::aes(x = year, y = median),
        linetype = "dashed",
        colour = "black",
        linewidth = 0.7,
        inherit.aes = FALSE
      )
  }

  p
}

plot_placebo <- function(resplacebo,
                         x_limits = c(1950, 1975),
                         y_limits = NULL,
                         y_breaks = waiver(),
                         y_label = "Differential GDP per capita",
                         y_formatter = scales::dollar_format()) {
  xlim_dates <- NULL
  if (!is.null(x_limits)) {
    xlim_dates <- as.Date(paste0(x_limits, "-01-01"))
  }

  p <- suppressWarnings(
    ggplot2::ggplot(resplacebo,
      exclude.ratio = 5,
      ratio.type = "mspe",
      ylab = y_label,
      xlab = "Year",
      main = "Figure 5 Placebo Test"
    )
  )

  p +
    ggplot2::theme(
      axis.ticks = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 0.5, vjust = 0.5),
      axis.text = ggplot2::element_text(face = "bold")
    ) +
    ggplot2::scale_x_date(
      expand = ggplot2::expansion(mult = c(0, 0.01)),
      date_breaks = "5 years",
      minor_breaks = NULL,
      date_labels = "%Y"
    ) +
    ggplot2::coord_cartesian(xlim = xlim_dates) +
    ggplot2::scale_y_continuous(
      limits = y_limits,
      breaks = y_breaks,
      minor_breaks = waiver(),
      labels = y_formatter,
      expand = expansion(mult = c(0, 0.02))
    )
}

plot_post_pre_ratio <- function(ppratio_df) {
  ggplot2::ggplot() +
    ggplot2::geom_col(
      data = ppratio_df,
      ggplot2::aes(ppratio, y = reorder(Country, ppratio)),
      fill = "black"
    ) +
    ggplot2::geom_col(
      data = ppratio_df[1, ],
      ggplot2::aes(x = ppratio, y = reorder(Country, ppratio)),
      fill = "red"
    ) +
    ggplot2::labs(
      x = NULL,
      y = "Country",
      title = "Figure 6 post-pre MSPE ratio"
    ) +
    ggplot2::coord_cartesian(xlim = c(0, 700), expand = FALSE) +
    ggplot2::scale_x_continuous(
      breaks = seq(0, 700, 100),
      minor_breaks = NULL,
      expand = c(0, 0)
    )
}

