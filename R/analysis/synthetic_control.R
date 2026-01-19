run_synthetic_control <- function(session_dir) {
  plots_dir <- file.path(session_dir, OUTPUT_PLOTS)
  tables_dir <- file.path(session_dir, OUTPUT_TABLES)
  ensure_dir(plots_dir)
  ensure_dir(tables_dir)

  safe_ggsave <- function(path, plot_obj, width = 20, height = 11, units = "cm", dpi = 600) {
    tryCatch(
      ggplot2::ggsave(path, plot_obj, width = width, height = height, units = units, dpi = dpi),
      error = function(e) message(sprintf("ggsave failed for %s: %s", path, e$message))
    )
  }

  write_run_metadata <- function(session_dir, metadata) {
    val <- function(x) {
      if (length(x) == 0 || is.na(x) || is.null(x)) return("")
      as.character(x)
    }
    to_block <- function(name, items) {
      if (length(items) == 0) return(character(0))
      c(sprintf("%s:", name), sprintf(" - %s", items))
    }
    lines <- c(
      sprintf("timestamp: \"%s\"", val(metadata$timestamp)),
      sprintf("session_dir: \"%s\"", val(metadata$session_dir)),
      sprintf("r_version: \"%s\"", val(metadata$r_version)),
      sprintf("platform: \"%s\"", val(metadata$platform)),
      sprintf("renv_version: \"%s\"", val(metadata$renv_version)),
      sprintf("placebo_mspe_ratio_cutoff: \"%s\"", val(metadata$placebo_mspe_ratio_cutoff)),
      "env_flags:",
      sprintf(" OUTCOME_ONLY: \"%s\"", val(metadata$env_flags$OUTCOME_ONLY)),
      sprintf(" SPEC_ONLY: \"%s\"", val(metadata$env_flags$SPEC_ONLY)),
      sprintf(" SKIP_DROP_ONE: \"%s\"", val(metadata$env_flags$SKIP_DROP_ONE)),
      sprintf(" placebo_mspe_ratio_cutoff_env: \"%s\"", val(metadata$env_flags$placebo_mspe_ratio_cutoff))
    )
    lines <- c(
      lines,
      to_block("outcomes", metadata$outcomes),
      to_block("specs", metadata$specs)
    )
    if (!is.null(metadata$data_files) && length(metadata$data_files) > 0) {
      lines <- c(lines, "data_files:")
      df_lines <- unlist(lapply(metadata$data_files, function(df) {
        c(
          sprintf(" - key: \"%s\"", val(df$key)),
          sprintf(" path: \"%s\"", val(df$path)),
          sprintf(" expected_md5: \"%s\"", val(df$expected_md5))
        )
      }))
      lines <- c(lines, df_lines)
    }
    ensure_dir(session_dir)
    writeLines(lines, file.path(session_dir, "run_config.yml"))
  }

  extract_weights <- function(res, controls_identifier) {
    safe_get <- function(expr) tryCatch(expr, error = function(e) NULL)
    candidates <- list(
      safe_get(res$solution.w),
      safe_get(res$solution$w),
      safe_get(res$w),
      safe_get(res$weights),
      safe_get(res$W),
      safe_get(res$Results$W)
    )
    for (w in candidates) {
      if (is.null(w)) next
      w_mat <- as.matrix(w)
      if (nrow(w_mat) == length(controls_identifier)) {
        w_vec <- w_mat[, 1]
        nms <- rownames(w_mat)
      } else if (ncol(w_mat) == length(controls_identifier)) {
        w_vec <- w_mat[1, ]
        nms <- colnames(w_mat)
      } else if (length(w_mat) == length(controls_identifier)) {
        w_vec <- as.numeric(w_mat)
        nms <- names(w_mat)
      } else {
        next
      }
      if (is.null(nms) || any(is.na(nms)) || any(nms == "")) {
        nms <- controls_identifier
      }
      names(w_vec) <- nms
      return(w_vec)
    }
    NULL
  }

  extract_gap_matrix <- function(obj, dep_var = NULL) {
    candidates <- list(
      obj$gaps,
      obj$gap,
      obj$res_gaps,
      obj$Results$gaps,
      obj$Results$gap
    )
    for (g in candidates) {
      if (is.null(g)) next
      if (is.vector(g) && !is.list(g)) {
        return(matrix(as.numeric(g), ncol = 1, dimnames = list(names(g), NULL)))
      }
      if (is.matrix(g)) return(g)
      if (is.data.frame(g)) return(as.matrix(g))
    }
    # mscmt objects often store gaps as a list per variable; pick the dependent var
    if (inherits(obj, "mscmt") && is.list(obj$gaps)) {
      dep <- dep_var
      if (is.null(dep) && !is.null(obj$dependent)) dep <- obj$dependent
      if (!is.null(dep) && !is.null(obj$gaps[[dep]])) {
        g <- as.numeric(obj$gaps[[dep]])
        return(matrix(g, ncol = 1, dimnames = list(names(obj$gaps[[dep]]), NULL)))
      }
    }
    # MSCMT placebo objects may come back as a list of mscmt objects (one per unit)
    if (is.list(obj) && any(vapply(obj, inherits, logical(1), "mscmt"))) {
      valid_idx <- which(vapply(obj, inherits, logical(1), "mscmt"))
      if (length(valid_idx) == 0) return(NULL)
      dep <- dep_var
      if (is.null(dep)) {
        dep <- obj[[valid_idx[1]]]$dependent
      }
      if (is.null(dep)) return(NULL)
      # Collect gaps for the dependent variable from each mscmt element
      gap_list <- lapply(obj[valid_idx], function(el) {
        if (!inherits(el, "mscmt")) return(NULL)
        if (is.null(el$gaps)) return(NULL)
        g <- el$gaps[[dep]]
        if (is.null(g)) return(NULL)
        as.numeric(g)
      })
      # Remove NULLs and align lengths
      valid_idx <- which(!vapply(gap_list, is.null, logical(1)))
      if (length(valid_idx) == 0) return(NULL)
      gap_list <- gap_list[valid_idx]
      lens <- vapply(gap_list, length, integer(1))
      common_len <- min(lens)
      gap_mat <- do.call(cbind, lapply(gap_list, function(g) g[seq_len(common_len)]))
      colnames(gap_mat) <- names(obj)[valid_idx]
      return(gap_mat)
    }
    NULL
  }

  resolve_years <- function(gap_mat, fallback_years) {
    yrs <- fallback_years
    if (!is.null(rownames(gap_mat))) {
      suppressWarnings({
        rn_years <- as.numeric(rownames(gap_mat))
      })
      if (!any(is.na(rn_years))) {
        yrs <- rn_years
      }
    }
    if (length(yrs) != nrow(gap_mat)) {
      yrs <- seq_len(nrow(gap_mat))
    }
    yrs
  }

  gap_matrix_to_long <- function(gap_mat, years, unit_names = NULL) {
    if (is.null(colnames(gap_mat)) || any(colnames(gap_mat) == "")) {
      if (!is.null(unit_names) && length(unit_names) == ncol(gap_mat)) {
        colnames(gap_mat) <- unit_names
      } else {
        colnames(gap_mat) <- paste0("unit_", seq_len(ncol(gap_mat)))
      }
    }
    dplyr::as_tibble(gap_mat) |>
      dplyr::mutate(year = years) |>
      tidyr::pivot_longer(cols = -year, names_to = "Country", values_to = "gap")
  }

  summarize_gap_metrics <- function(gaps_long, pre_window, post_window, treatment_identifier, spec_name, pool_label, outcome) {
    if (is.null(gaps_long)) return(NULL)
    treated <- dplyr::filter(gaps_long, Country == treatment_identifier)
    if (nrow(treated) == 0) return(NULL)
    pre_mask <- treated$year >= pre_window[1] & treated$year <= pre_window[2]
    post_mask <- treated$year >= post_window[1] & treated$year <= post_window[2]
    tibble::tibble(
      outcome = outcome,
      spec = spec_name,
      pool = pool_label,
      Country = treatment_identifier,
      pre_mspe = mean(treated$gap[pre_mask]^2, na.rm = TRUE),
      post_mspe = mean(treated$gap[post_mask]^2, na.rm = TRUE),
      avg_gap_post = mean(treated$gap[post_mask], na.rm = TRUE),
      n_pre = sum(pre_mask & !is.na(treated$gap)),
      n_post = sum(post_mask & !is.na(treated$gap))
    )
  }

  rank_pvalue <- function(treated_value, placebo_values, direction = "greater") {
    placebo_values <- placebo_values[is.finite(placebo_values)]
    if (is.null(treated_value) || length(treated_value) == 0) return(NA_real_)
    treated_value <- treated_value[1]
    if (!is.finite(treated_value)) return(NA_real_)
    if (length(placebo_values) == 0) {
      return(NA_real_)
    }
    if (direction == "greater") {
      (1 + sum(placebo_values >= treated_value)) / (length(placebo_values) + 1)
    } else {
      (1 + sum(placebo_values <= treated_value)) / (length(placebo_values) + 1)
    }
  }

  build_placebo_outputs <- function(resplacebo, years, unit_names, treatment_identifier, pre_period, post_period, tables_dir, dep_var = NULL) {
    gap_mat <- extract_gap_matrix(resplacebo, dep_var = dep_var)
    if (is.null(gap_mat)) {
      return(list(gaps_long = NULL, band_df = NULL, mspe_df = NULL, placebo_only = NULL))
    }
    years <- resolve_years(gap_mat, years)
    colnames(gap_mat) <- if (is.null(colnames(gap_mat))) unit_names else colnames(gap_mat)
    gaps_long <- gap_matrix_to_long(gap_mat, years, unit_names)
    readr::write_csv(gaps_long, file.path(tables_dir, "placebo_gaps_long.csv"))

    placebo_only <- dplyr::filter(gaps_long, Country != treatment_identifier, Country != "Average")

    band_df <- placebo_only |>
      dplyr::group_by(year) |>
      dplyr::summarise(
        q05 = stats::quantile(gap, 0.05, na.rm = TRUE),
        q50 = stats::quantile(gap, 0.50, na.rm = TRUE),
        q95 = stats::quantile(gap, 0.95, na.rm = TRUE),
        .groups = "drop"
      )
    readr::write_csv(band_df, file.path(tables_dir, "placebo_gap_quantiles.csv"))

    placebo_periods <- placebo_only |>
      dplyr::mutate(
        period = dplyr::case_when(
          year >= pre_period[1] & year <= pre_period[2] ~ "pre",
          year >= post_period[1] & year <= post_period[2] ~ "post",
          TRUE ~ NA_character_
        )
      ) |>
      dplyr::filter(!is.na(period))

    mspe_df <- NULL
    if (nrow(placebo_periods) > 0) {
      mspe_df <- placebo_periods |>
        dplyr::group_by(Country, period) |>
        dplyr::summarise(mspe = mean(gap^2, na.rm = TRUE), .groups = "drop") |>
        tidyr::pivot_wider(names_from = period, values_from = mspe)
      # Si falta pre o post, completa con NA para evitar errores downstream
      if (!"pre" %in% names(mspe_df)) mspe_df$pre <- NA_real_
      if (!"post" %in% names(mspe_df)) mspe_df$post <- NA_real_
      mspe_df <- dplyr::mutate(mspe_df, mspe_ratio = post / pre)
      readr::write_csv(mspe_df, file.path(tables_dir, "placebo_mspe.csv"))
    }

    list(gaps_long = gaps_long, band_df = band_df, mspe_df = mspe_df, placebo_only = placebo_only)
  }

  conformal_bands_from_placebos <- function(treated_gaps, placebo_long, alpha = 0.1, treatment_identifier = "Spain") {
    if (is.null(treated_gaps) || is.null(placebo_long)) {
      return(list(band_df = NULL, pvals_df = NULL))
    }
    treated <- dplyr::filter(treated_gaps, Country == treatment_identifier)
    placebo_only <- placebo_long |>
      dplyr::filter(Country != treatment_identifier, Country != "Average")
    if (nrow(treated) == 0 || nrow(placebo_only) == 0) {
      return(list(band_df = NULL, pvals_df = NULL))
    }
    calc <- lapply(seq_len(nrow(treated)), function(i) {
      g_t <- treated$gap[i]
      y_t <- treated$year[i]
      pool <- placebo_only$gap[placebo_only$year == y_t]
      pool <- pool[!is.na(pool)]
      if (length(pool) == 0 || is.na(g_t)) {
        return(NULL)
      }
      q <- stats::quantile(abs(pool), probs = 1 - alpha, na.rm = TRUE, names = FALSE)
      p_val <- (1 + sum(abs(pool) >= abs(g_t)))/(length(pool) + 1)
      dplyr::tibble(
        year = y_t,
        treated_gap = g_t,
        conformal_q = q,
        lower = g_t - q,
        upper = g_t + q,
        placebo_n = length(pool),
        p_value = p_val
      )
    })
    out <- dplyr::bind_rows(calc)
    if (nrow(out) == 0) {
      return(list(band_df = NULL, pvals_df = NULL))
    }
    band_df <- out |>
      dplyr::select(year, lower, upper, median = treated_gap, conformal_q, placebo_n, p_value)
    pvals_df <- out |>
      dplyr::select(year, p_value, placebo_n)
    list(band_df = band_df, pvals_df = pvals_df)
  }

  prepare_panel_for_spec <- function(panel_full, pre_window, post_window, key_vars, outcome_var, spec_name, pool_label, tables_dir, treatment_identifier, donor_filter = NULL, write_outputs = TRUE) {
    df <- panel_full
    if (!is.null(donor_filter)) {
      df <- dplyr::filter(
        df,
        country %in% donor_filter | country %in% treatment_identifier | country == "Average"
      )
    }
    base_pool <- unique(df$country[df$country != treatment_identifier & df$country != "Average"])
    df <- dplyr::filter(df, year >= pre_window[1], year <= post_window[2])
    pre_len <- pre_window[2] - pre_window[1] + 1
    post_len <- post_window[2] - post_window[1] + 1

    coverage_pre <- df |>
      dplyr::filter(country != "Average", year >= pre_window[1], year <= pre_window[2]) |>
      dplyr::group_by(country) |>
      dplyr::summarise(
        dplyr::across(dplyr::all_of(key_vars), \(x) sum(!is.na(x)), .names = "n_{.col}"),
        .groups = "drop"
      )

    coverage_post <- df |>
      dplyr::filter(country != "Average", year >= post_window[1], year <= post_window[2]) |>
      dplyr::group_by(country) |>
      dplyr::summarise(
        n_outcome = sum(!is.na(.data[[outcome_var]])),
        .groups = "drop"
      )

    if (write_outputs) {
      readr::write_csv(coverage_pre, file.path(tables_dir, "pwt_coverage_pre.csv"))
      save_table_jpg(coverage_pre, file.path(tables_dir, "pwt_coverage_pre.jpg"), title = sprintf("PWT coverage pre (%s-%s)", pre_window[1], pre_window[2]))
      readr::write_csv(coverage_post, file.path(tables_dir, "pwt_coverage_post_outcome.csv"))
      save_table_jpg(coverage_post, file.path(tables_dir, "pwt_coverage_post_outcome.jpg"), title = sprintf("PWT coverage post outcome (%s-%s)", post_window[1], post_window[2]))
    }

    excluded_pre <- coverage_pre |>
      tidyr::pivot_longer(cols = dplyr::starts_with("n_"), names_to = "var", values_to = "non_na") |>
      dplyr::mutate(var = gsub("^n_", "", var)) |>
      dplyr::filter(non_na < pre_len) |>
      dplyr::mutate(reason = paste0("Insufficient pre data for ", var, " (", non_na, "/", pre_len, ")")) |>
      dplyr::select(country, reason) |>
      dplyr::distinct()

    excluded_post <- coverage_post |>
      dplyr::filter(n_outcome < post_len) |>
      dplyr::mutate(
        var = outcome_var,
        reason = paste0("Insufficient post data for ", outcome_var, " (", n_outcome, "/", post_len, ")")
      ) |>
      dplyr::select(country, reason) |>
      dplyr::distinct()

    excluded <- dplyr::bind_rows(excluded_pre, excluded_post) |>
      dplyr::group_by(country) |>
      dplyr::summarise(reason = paste(unique(reason), collapse = "; "), .groups = "drop")

    if (nrow(excluded) > 0) {
      if (write_outputs) {
        readr::write_csv(excluded, file.path(tables_dir, "excluded_countries.csv"))
        save_table_jpg(excluded, file.path(tables_dir, "excluded_countries.jpg"), title = "Excluded countries (coverage)")
        if (nrow(excluded_pre) > 0) {
          readr::write_csv(excluded_pre, file.path(tables_dir, "excluded_countries_pre.csv"))
          save_table_jpg(excluded_pre, file.path(tables_dir, "excluded_countries_pre.jpg"), title = "Excluded (pre, predictors)")
        }
        if (nrow(excluded_post) > 0) {
          readr::write_csv(excluded_post, file.path(tables_dir, "excluded_countries_post_outcome.csv"))
          save_table_jpg(excluded_post, file.path(tables_dir, "excluded_countries_post_outcome.jpg"), title = "Excluded (post, outcome)")
        }
      }
      df <- dplyr::filter(df, !country %in% excluded$country)
    }

    panel_no_avg <- df |>
      dplyr::filter(country != "Average")
    panel_donors_only <- panel_no_avg |>
      dplyr::filter(country != treatment_identifier)
    if (nrow(panel_donors_only) == 0) stop(sprintf("No donor countries remain after coverage filters for spec %s / pool %s.", spec_name, pool_label))

    # Pool estable de donantes: se fija una sola vez tras los filtros de cobertura
    remaining_pool <- sort(unique(panel_donors_only$country))

    region_avg <- max(panel_donors_only$region, na.rm = TRUE) + 1
    if (!is.finite(region_avg)) region_avg <- 1

    panel_included <- panel_no_avg |>
      dplyr::filter(country == treatment_identifier | country %in% remaining_pool)

    numeric_cols <- setdiff(names(dplyr::select(panel_donors_only, where(is.numeric))), c("region", "year", "pop"))
    avg_clean <- panel_donors_only |>
      dplyr::filter(country %in% remaining_pool) |>
      dplyr::group_by(year) |>
      dplyr::summarise(
        pop = sum(pop, na.rm = TRUE),
        dplyr::across(
          dplyr::all_of(numeric_cols),
          ~ {
            wts <- dplyr::cur_data_all()$pop
            valid <- !is.na(.x) & !is.na(wts)
            if (!any(valid)) return(NA_real_)
            stats::weighted.mean(.x[valid], w = wts[valid])
          }
        ),
        .groups = "drop"
      ) |>
      dplyr::mutate(country = "Average", isocode = "AVG", region = region_avg)

    panel_df <- dplyr::bind_rows(avg_clean, panel_included) |>
      dplyr::mutate(country = forcats::fct_relevel(country, "Average")) |>
      dplyr::arrange(country, year)

    panel_df$region <- as.integer(panel_df$region)
    panel_df <- as.data.frame(panel_df)
    panel_df$country <- as.character(panel_df$country)
    panel_df$region <- as.numeric(panel_df$region)

    pool_status <- tibble::tibble(
      country = sort(unique(c(base_pool, remaining_pool, excluded$country))),
      status = dplyr::if_else(country %in% remaining_pool, "included", "excluded"),
      phase = dplyr::case_when(
        country %in% excluded_pre$country ~ "pre_predictors",
        country %in% excluded_post$country ~ "post_outcome",
        TRUE ~ NA_character_
      )
    ) |>
      dplyr::left_join(
        dplyr::bind_rows(
          dplyr::mutate(excluded_pre, phase = "pre_predictors"),
          dplyr::mutate(excluded_post, phase = "post_outcome")
        ),
        by = c("country", "phase")
      )

    if (write_outputs) {
      final_pool <- tibble::tibble(country = remaining_pool)
      readr::write_csv(pool_status, file.path(tables_dir, "pool_status.csv"))
      save_table_jpg(pool_status, file.path(tables_dir, "pool_status.jpg"), title = "Pool status (included/excluded)")
      readr::write_csv(final_pool, file.path(tables_dir, "final_pool.csv"))
      save_table_jpg(final_pool, file.path(tables_dir, "final_pool.jpg"), title = "Final donor pool")
    }

    pool_status_rows[[length(pool_status_rows) + 1]] <<- dplyr::mutate(
      pool_status,
      outcome = outcome$id,
      spec = spec$name,
      pool = pool_label
    ) |>
      dplyr::relocate(outcome, spec, pool)

    pool_summary_rows[[length(pool_summary_rows) + 1]] <<- tibble::tibble(
      outcome = outcome$id,
      spec = spec$name,
      pool = pool_label,
      donors_included = sum(pool_status$status == "included", na.rm = TRUE),
      donors_excluded = sum(pool_status$status != "included", na.rm = TRUE),
      exclusion_reasons = paste(sort(unique(stats::na.omit(pool_status$reason))), collapse = "; ")
    )

    list(
      panel = panel_df,
      coverage_pre = coverage_pre,
      coverage_post = coverage_post,
      excluded = excluded,
      pool_status = pool_status,
      key_vars = key_vars
    )
  }

  specs <- list(
    list(name = "baseline", label = "Baseline (1950-1959 pre, 1960-1975 post)", pre_window = c(1950, 1959), post_window = c(1960, 1975), save_legacy = TRUE),
    list(name = "treat_1970", label = "Treatment 1970 (1960-1969 pre, 1970-1985 post)", pre_window = c(1960, 1969), post_window = c(1970, 1985), save_legacy = FALSE),
    list(name = "treat_1980", label = "Treatment 1980 (1970-1979 pre, 1980-1995 post)", pre_window = c(1970, 1979), post_window = c(1980, 1995), save_legacy = FALSE)
  )

  spec_env <- Sys.getenv("SPEC_ONLY", "")
  if (nzchar(spec_env)) {
    wanted_specs <- trimws(strsplit(spec_env, ",")[[1]])
    specs <- Filter(function(s) s$name %in% wanted_specs, specs)
  }

  outcomes <- list(
    list(
      id = "gdpcap",
      dep_var = "gdpcap",
      y_label = "Real GDP per capita",
      y_formatter = scales::dollar_format(),
      legacy = TRUE
    ),
    list(
      id = "rknacapita",
      dep_var = "rknacapita",
      y_label = "Real capital stock per capita (rknacapita)",
      y_formatter = scales::dollar_format(),
      legacy = FALSE
    ),
    list(
      id = "hc",
      dep_var = "hc",
      y_label = "Human capital index (PWT hc)",
      y_formatter = scales::number_format(accuracy = 0.01),
      legacy = FALSE
    )
  )

  outcome_env <- Sys.getenv("OUTCOME_ONLY", "")
  if (nzchar(outcome_env)) {
    wanted <- strsplit(outcome_env, ",")[[1]]
    wanted <- trimws(wanted)
    outcomes <- Filter(function(o) o$id %in% wanted, outcomes)
  }

  env_flags <- Sys.getenv(c("OUTCOME_ONLY", "SPEC_ONLY", "SKIP_DROP_ONE", "placebo_mspe_ratio_cutoff"))
  checksum_info <- tryCatch(read_checksum_file(), error = function(e) e)
  checksum_lookup <- function(key) {
    if (inherits(checksum_info, "error") || is.null(checksum_info[[key]])) return(NA_character_)
    checksum_info[[key]]
  }

  metadata <- list(
    timestamp = basename(session_dir),
    session_dir = normalizePath(session_dir, winslash = "/", mustWork = FALSE),
    r_version = paste(R.version$major, R.version$minor, sep = "."),
    platform = R.version$platform,
    renv_version = tryCatch(as.character(utils::packageVersion("renv")), error = function(e) NA_character_),
    placebo_mspe_ratio_cutoff = placebo_mspe_ratio_cutoff,
    env_flags = as.list(env_flags),
    outcomes = vapply(outcomes, `[[`, character(1), "id"),
    specs = vapply(specs, `[[`, character(1), "name"),
    data_files = list(
      list(key = "pwt11", path = PWT11_LOCAL, expected_md5 = checksum_lookup(PWT11_CHECKSUM_KEY)),
      list(key = "maddison_mpd2023", path = MADDISON_LOCAL, expected_md5 = checksum_lookup(MADDISON_CHECKSUM_KEY))
    )
  )
  write_run_metadata(session_dir, metadata)

  max_post_year <- max(vapply(specs, function(s) s$post_window[2], numeric(1)))
  raw <- load_pwt_data()
  panel_full <- prepare_pwt11_panel(raw, max_year = max_post_year)
  panel_full$region <- as.numeric(panel_full$region)
  panel_full <- as.data.frame(panel_full)
  panel_full$country <- as.character(panel_full$country)
  panel_full$region <- as.numeric(panel_full$region)

  treatment_identifier <- "Spain"
  base_key_vars <- c("pop", "gdpcap", "hc", "csh_c", "csh_i", "csh_g", "csh_x", "csh_m", "rknacapita", "pl_gdpo")

  metrics_accum <- list()
  predictor_tables <- list()
  donor_weights_summary <- list()
  donor_weights_rows <- list()
  pool_status_rows <- list()
  pool_summary_rows <- list()
  placebo_filter_rows <- list()
  spec_outputs <- list()
  pval_accum <- list()

  for (outcome in outcomes) {
    key_vars <- base_key_vars
    outcome_root_tables <- file.path(tables_dir, paste0("outcome_", outcome$id))
    outcome_root_plots <- file.path(plots_dir, paste0("outcome_", outcome$id))
    ensure_dir(outcome_root_tables)
    ensure_dir(outcome_root_plots)

    spec_outputs[[outcome$id]] <- list()

    for (spec in specs) {
      spec_tables_root <- file.path(outcome_root_tables, "specs", spec$name)
      spec_plots_root <- file.path(outcome_root_plots, "specs", spec$name)
      ensure_dir(spec_tables_root)
      ensure_dir(spec_plots_root)

      pool_results <- list()
      top_donor <- NULL
      dep_var <- outcome$dep_var

      run_spec_pool <- function(pool_label, donor_filter, save_legacy_figs = FALSE) {
        pool_tables <- file.path(spec_tables_root, pool_label)
        pool_plots <- file.path(spec_plots_root, pool_label)
        ensure_dir(pool_tables)
        ensure_dir(pool_plots)

        panel_res <- prepare_panel_for_spec(
          panel_full = panel_full,
          pre_window = spec$pre_window,
          post_window = spec$post_window,
          key_vars = key_vars,
          outcome_var = dep_var,
          spec_name = spec$name,
          pool_label = pool_label,
          tables_dir = pool_tables,
          treatment_identifier = treatment_identifier,
          donor_filter = donor_filter
        )

        key_vars <- panel_res$key_vars
        panel_df <- panel_res$panel
        keep_vars <- unique(c("country", "isocode", "region", "year", key_vars))
        panel_df <- panel_df[, intersect(names(panel_df), keep_vars), drop = FALSE]

        datprep <- MSCMT::listFromLong(
          panel_df,
          unit.variable = which(names(panel_df) == "region"),
          time.variable = which(names(panel_df) == "year"),
          unit.names.variable = which(names(panel_df) == "country")
        )

        controls_identifier <- setdiff(colnames(datprep[[1]]), c(treatment_identifier, "Average"))
        pre_period <- range(spec$pre_window)
        post_period <- range(spec$post_window)
        pred_vars <- key_vars

        times_dep <- matrix(pre_period, ncol = 1)
        colnames(times_dep) <- dep_var
        times_pred <- do.call(cbind, lapply(pred_vars, function(x) pre_period))
        colnames(times_pred) <- pred_vars
        agg_fns <- rep("mean", ncol(times_pred))

        res <- suppressWarnings(
          MSCMT::mscmt(
            datprep,
            treatment.identifier = treatment_identifier,
            controls.identifier = controls_identifier,
            times.dep = times_dep,
            times.pred = times_pred,
            agg.fns = agg_fns,
            seed = 1,
            outer.optim = "DEoptim",
            verbose = FALSE
          )
        )

        predictor_table <- as.data.frame(res$predictor.table)
        predictor_table <- tibble::rownames_to_column(predictor_table, var = "Predictor")
        colnames(predictor_table) <- c("Predictor", "Spain", "Synthetic Spain", "Sample")
        predictor_table <- predictor_table |>
          dplyr::mutate(
            diff_synth = Spain - `Synthetic Spain`,
            diff_sample = Spain - Sample,
            improvement_pct = dplyr::if_else(
              abs(diff_sample) < .Machine$double.eps,
              NA_real_,
              (abs(diff_sample) - abs(diff_synth)) / abs(diff_sample) * 100
            )
          ) |>
          dplyr::mutate(dplyr::across(where(is.numeric), ~round(.x, 4)))
        predictor_tables[[paste(outcome$id, spec$name, pool_label, sep = "_")]] <<- predictor_table
        readr::write_csv(predictor_table, file.path(pool_tables, "predictor_table.csv"))
        save_table_jpg(predictor_table, file.path(pool_tables, "predictor_table.jpg"), title = sprintf("Predictor balance (%s / %s / %s)", outcome$id, spec$name, pool_label))

        w_vec <- extract_weights(res, controls_identifier)
        donor_weights <- NULL
        if (!is.null(w_vec)) {
          donor_weights <- tibble::tibble(
            Country = names(w_vec),
            Weight = as.numeric(w_vec)
          ) |>
            dplyr::filter(Weight > 0) |>
            dplyr::arrange(dplyr::desc(Weight))
          readr::write_csv(donor_weights, file.path(pool_tables, "donor_weights.csv"))
          save_table_jpg(donor_weights, file.path(pool_tables, "donor_weights.jpg"), title = sprintf("Donor weights (%s / %s / %s)", outcome$id, spec$name, pool_label))
          donor_weights_summary[[paste(outcome$id, spec$name, pool_label, sep = "_")]] <<- donor_weights
          donor_weights_rows[[length(donor_weights_rows) + 1]] <<- dplyr::mutate(
            donor_weights,
            outcome = outcome$id,
            spec = spec$name,
            pool = pool_label
          ) |>
            dplyr::relocate(outcome, spec, pool)
        }

        if (pool_label == "all" && !is.null(donor_weights) && nrow(donor_weights) > 0) {
          top_donor <<- donor_weights$Country[1]
        }

        gap_mat_main <- extract_gap_matrix(res, dep_var = dep_var)
        if (!is.null(gap_mat_main) && ncol(gap_mat_main) == 1) {
          colnames(gap_mat_main) <- treatment_identifier
        }
        gaps_long_main <- NULL
        if (!is.null(gap_mat_main)) {
          years_guess <- suppressWarnings(as.numeric(rownames(datprep[[1]])))
          gaps_long_main <- gap_matrix_to_long(
            gap_mat_main,
            years = resolve_years(gap_mat_main, years_guess),
            unit_names = colnames(datprep[[1]])
          )
        }

        metrics_main <- summarize_gap_metrics(
          gaps_long = gaps_long_main,
          pre_window = pre_period,
          post_window = post_period,
          treatment_identifier = treatment_identifier,
          spec_name = spec$name,
          pool_label = pool_label,
          outcome = outcome$id
        )
        if (!is.null(metrics_main)) metrics_accum[[length(metrics_accum) + 1]] <<- metrics_main

        x_limits <- c(spec$pre_window[1], spec$post_window[2])
        gap_ylim <- NULL
        if (!is.null(gaps_long_main)) {
          g_range <- range(gaps_long_main$gap, na.rm = TRUE)
          g_pad <- diff(g_range) * 0.08
          if (!is.finite(g_pad)) g_pad <- 100
          gap_ylim <- c(g_range[1] - g_pad, g_range[2] + g_pad)
        }
        comp_ylim <- NULL
        if (!is.null(gap_mat_main) && !is.null(res$Y1plot)) {
          comp_range <- range(res$Y1plot, res$Y1plot + gap_mat_main, na.rm = TRUE)
          comp_pad <- diff(comp_range) * 0.05
          if (!is.finite(comp_pad)) comp_pad <- 500
          comp_ylim <- c(comp_range[1] - comp_pad, comp_range[2] + comp_pad)
        }

        p3 <- plot_synth_comparison(res, x_limits = x_limits, y_limits = comp_ylim, y_label = outcome$y_label, y_formatter = outcome$y_formatter)
        p4 <- plot_synth_gaps(res, x_limits = x_limits, y_limits = gap_ylim, y_label = paste0(outcome$y_label, " gap"), y_formatter = outcome$y_formatter, treatment_year = post_period[1])

        if (save_legacy_figs && outcome$legacy) {
          safe_ggsave(file.path(plots_dir, "Figure_3_synth_comparison.png"), p3)
          safe_ggsave(file.path(plots_dir, "Figure_4_synth_gaps.png"), p4)
        }
        safe_ggsave(file.path(pool_plots, "Figure_3_synth_comparison.png"), p3)
        safe_ggsave(file.path(pool_plots, "Figure_4_synth_gaps.png"), p4)

        cl <- tryCatch(parallel::makeCluster(max(1, min(4, parallel::detectCores()))), error = function(e) NULL)
        resplacebo <- suppressWarnings(
          MSCMT::mscmt(
            datprep,
            treatment.identifier = treatment_identifier,
            controls.identifier = controls_identifier,
            times.dep = times_dep,
            times.pred = times_pred,
            agg.fns = agg_fns,
            cl = cl,
            placebo = TRUE,
            seed = 1,
            outer.optim = "DEoptim"
          )
        )
        if (!is.null(cl)) parallel::stopCluster(cl)

        time_index <- suppressWarnings(as.numeric(rownames(datprep[[1]])))
        if (any(is.na(time_index))) time_index <- sort(unique(panel_df$year))

        placebo_outputs <- build_placebo_outputs(
          resplacebo,
          years = time_index,
          unit_names = colnames(datprep[[1]]),
          treatment_identifier = treatment_identifier,
          pre_period = pre_period,
          post_period = post_period,
          tables_dir = pool_tables,
          dep_var = dep_var
        )

        # Filtro ex ante: descarta placebos con RMSPE pre muy superior al tratado (Abadie et al. 2010)
        if (!is.null(placebo_outputs$gaps_long)) {
          pre_mspe_df <- placebo_outputs$gaps_long |>
            dplyr::filter(year >= pre_period[1], year <= pre_period[2]) |>
            dplyr::group_by(Country) |>
            dplyr::summarise(pre_mspe = mean(gap^2, na.rm = TRUE), .groups = "drop")
          treated_mspe <- pre_mspe_df$pre_mspe[pre_mspe_df$Country == treatment_identifier]
          cutoff <- placebo_mspe_ratio_cutoff
          use_filter <- length(treated_mspe) > 0 && is.finite(treated_mspe) && treated_mspe > 0 &&
            is.finite(cutoff) && cutoff > 0
          pre_mspe_df <- pre_mspe_df |>
            dplyr::mutate(
              ratio_vs_treated = if (use_filter) pre_mspe / treated_mspe else NA_real_,
              keep = if (use_filter) (is.na(ratio_vs_treated) | ratio_vs_treated <= cutoff) else TRUE
            )
          good_units <- pre_mspe_df$Country[pre_mspe_df$keep | pre_mspe_df$Country == treatment_identifier | pre_mspe_df$Country == "Average"]
          good_units <- unique(c(good_units, treatment_identifier))
          placebos_total <- sum(pre_mspe_df$Country != treatment_identifier & pre_mspe_df$Country != "Average", na.rm = TRUE)
          placebos_kept <- sum(pre_mspe_df$Country %in% good_units & pre_mspe_df$Country != treatment_identifier & pre_mspe_df$Country != "Average", na.rm = TRUE)
          placebo_excluded <- pre_mspe_df |>
            dplyr::filter(!keep, Country != treatment_identifier, Country != "Average")
          if (nrow(placebo_excluded) > 0) {
            readr::write_csv(placebo_excluded, file.path(pool_tables, "placebo_excluded_badfit.csv"))
            save_table_jpg(placebo_excluded, file.path(pool_tables, "placebo_excluded_badfit.jpg"), title = sprintf("Placebos excluded (pre RMSPE > %sx treated)", cutoff))
          }
          message(sprintf(
            "Placebos retenidos: %d de %d (cutoff pre-RMSPE x tratado = %s, tratado=%.3f)",
            placebos_kept,
            placebos_total,
            if (!use_filter) "disabled" else format(cutoff, digits = 4),
            ifelse(length(treated_mspe) > 0, treated_mspe[1], NA_real_)
          ))
          placebo_filter_rows[[length(placebo_filter_rows) + 1]] <<- tibble::tibble(
            outcome = outcome$id,
            spec = spec$name,
            pool = pool_label,
            filter_applied = isTRUE(use_filter),
            cutoff = if (!use_filter) NA_real_ else cutoff,
            treated_pre_mspe = ifelse(length(treated_mspe) > 0, treated_mspe[1], NA_real_),
            placebos_total = placebos_total,
            placebos_kept = placebos_kept
          )
          gaps_long_filt <- placebo_outputs$gaps_long |>
            dplyr::filter(Country %in% good_units)
          placebo_only_filt <- gaps_long_filt |>
            dplyr::filter(Country != treatment_identifier, Country != "Average")
          band_df_filt <- NULL
          if (nrow(placebo_only_filt) > 0) {
            band_df_filt <- placebo_only_filt |>
              dplyr::group_by(year) |>
              dplyr::summarise(
                q05 = stats::quantile(gap, 0.05, na.rm = TRUE),
                q50 = stats::quantile(gap, 0.50, na.rm = TRUE),
                q95 = stats::quantile(gap, 0.95, na.rm = TRUE),
                .groups = "drop"
              )
          }
          mspe_df_filt <- NULL
          placebo_periods <- placebo_only_filt |>
            dplyr::mutate(
              period = dplyr::case_when(
                year >= pre_period[1] & year <= pre_period[2] ~ "pre",
                year >= post_period[1] & year <= post_period[2] ~ "post",
                TRUE ~ NA_character_
              )
            ) |>
            dplyr::filter(!is.na(period))
          if (nrow(placebo_periods) > 0) {
            mspe_df_filt <- placebo_periods |>
              dplyr::group_by(Country, period) |>
              dplyr::summarise(mspe = mean(gap^2, na.rm = TRUE), .groups = "drop") |>
              tidyr::pivot_wider(names_from = period, values_from = mspe)
            if (!"pre" %in% names(mspe_df_filt)) mspe_df_filt$pre <- NA_real_
            if (!"post" %in% names(mspe_df_filt)) mspe_df_filt$post <- NA_real_
            mspe_df_filt <- dplyr::mutate(mspe_df_filt, mspe_ratio = post / pre)
          }
          placebo_outputs$gaps_long <- gaps_long_filt
          placebo_outputs$band_df <- band_df_filt
          placebo_outputs$mspe_df <- mspe_df_filt
          placebo_outputs$placebo_only <- placebo_only_filt
          placebo_outputs$good_units <- good_units
          readr::write_csv(gaps_long_filt, file.path(pool_tables, "placebo_gaps_long.csv"))
          if (!is.null(band_df_filt)) readr::write_csv(band_df_filt, file.path(pool_tables, "placebo_gap_quantiles.csv"))
          if (!is.null(mspe_df_filt)) readr::write_csv(mspe_df_filt, file.path(pool_tables, "placebo_mspe.csv"))
        }

        if (is.null(placebo_outputs$gaps_long)) {
          placebo_filter_rows[[length(placebo_filter_rows) + 1]] <<- tibble::tibble(
            outcome = outcome$id,
            spec = spec$name,
            pool = pool_label,
            filter_applied = FALSE,
            cutoff = NA_real_,
            treated_pre_mspe = NA_real_,
            placebos_total = NA_real_,
            placebos_kept = NA_real_
          )
        }

        conformal <- conformal_bands_from_placebos(
          treated_gaps = gaps_long_main,
          placebo_long = placebo_outputs$gaps_long,
          alpha = 0.1,
          treatment_identifier = treatment_identifier
        )

        # Escala del eje Y en placebo plot: mantiene la escala fija para gdpcap,
        # pero usa escala dinamica para otros outcomes para evitar que se aplasten.
        p5_ylim <- NULL
        p5_breaks <- waiver()
        if (outcome$id == "gdpcap") {
          p5_range <- range(placebo_outputs$gaps_long$gap, na.rm = TRUE)
          max_abs <- suppressWarnings(max(abs(p5_range), na.rm = TRUE))
          if (!is.finite(max_abs)) max_abs <- 0
          target <- max(7500, ceiling(max_abs / 2500) * 2500)
          if (!is.finite(target) || target <= 0) target <- 7500
          p5_ylim <- c(-target, target)
          p5_breaks <- seq(-target, target, by = 2500)
        } else {
          combined_gaps <- c(
            if (!is.null(placebo_outputs$gaps_long)) placebo_outputs$gaps_long$gap else NA_real_,
            if (!is.null(gaps_long_main)) gaps_long_main$gap else NA_real_
          )
          combined_gaps <- combined_gaps[is.finite(combined_gaps)]
          if (length(combined_gaps) > 0) {
            g_range <- range(combined_gaps, na.rm = TRUE)
            pad <- diff(g_range) * 0.12
            if (!is.finite(pad) || pad == 0) pad <- max(abs(g_range), na.rm = TRUE) * 0.1
            if (!is.finite(pad) || pad == 0) pad <- 1
            p5_ylim <- c(g_range[1] - pad, g_range[2] + pad)
            p5_breaks <- scales::pretty_breaks(n = 8)(p5_ylim)
          }
        }

        p5 <- plot_placebo(
          resplacebo,
          x_limits = x_limits,
          y_limits = p5_ylim,
          y_breaks = p5_breaks,
          y_label = paste0(outcome$y_label, " gap"),
          y_formatter = outcome$y_formatter
        )
        safe_ggsave(file.path(pool_plots, "Figure_5_placebo.png"), p5)

        ppratio_df <- as.data.frame(MSCMT::ppratio(resplacebo, type = "mspe"))
        colnames(ppratio_df) <- c("ppratio")
        ppratio_df <- tibble::rownames_to_column(ppratio_df, var = "Country")
        readr::write_csv(ppratio_df, file.path(pool_tables, "post_pre_mspe_ratio.csv"))
        save_table_jpg(ppratio_df, file.path(pool_tables, "post_pre_mspe_ratio.jpg"), title = "Post/Pre MSPE ratio")

        pval_df <- tryCatch({
          vals <- MSCMT::pvalue(resplacebo, exclude.ratio = 5, ratio.type = "mspe", alternative = "less")
          tibble::enframe(vals, name = "Country", value = "p_value_less")
        }, error = function(e) NULL)
        if (!is.null(pval_df)) {
          readr::write_csv(pval_df, file.path(pool_tables, "placebo_pvalues.csv"))
          save_table_jpg(pval_df, file.path(pool_tables, "placebo_pvalues.jpg"), title = "Placebo p-values")
        }

        if (!is.null(conformal$band_df)) {
          readr::write_csv(conformal$band_df, file.path(pool_tables, "conformal_bands.csv"))
          save_table_jpg(conformal$band_df, file.path(pool_tables, "conformal_bands.jpg"), title = "Conformal bands (placebos espaciales)")
        }

        # Para graficar, priorizamos la banda de placebos (q05/q50/q95);
        # las bandas conformal se guardan en CSV/JPG pero no reemplazan la mediana placebo en la figura.
        band_for_plot <- placebo_outputs$band_df
        if (is.null(band_for_plot)) band_for_plot <- conformal$band_df
        if (!is.null(band_for_plot)) {
          band_for_plot <- dplyr::filter(band_for_plot, year >= post_period[1])
        }
        if (!is.null(band_for_plot)) {
          p4_band <- plot_synth_gaps(
            res,
            band_df = band_for_plot,
            x_limits = x_limits,
            y_limits = gap_ylim,
            y_label = paste0(outcome$y_label, " gap"),
            y_formatter = outcome$y_formatter,
            treatment_year = post_period[1]
          )
          safe_ggsave(file.path(pool_plots, "Figure_4_synth_gaps.png"), p4_band)
        }

        if (!is.null(placebo_outputs$band_df)) {
          save_table_jpg(placebo_outputs$band_df, file.path(pool_tables, "placebo_gap_quantiles.jpg"), title = "Placebo gap quantiles")
        }
        if (!is.null(placebo_outputs$mspe_df)) {
          save_table_jpg(placebo_outputs$mspe_df, file.path(pool_tables, "placebo_mspe.jpg"), title = "Placebo MSPE pre/post")
        }
        if (!is.null(placebo_outputs$gaps_long)) {
          save_table_jpg(placebo_outputs$gaps_long, file.path(pool_tables, "placebo_gaps_long.jpg"), title = "Placebo gaps (truncated)")
        }

        if (!is.null(placebo_outputs$gaps_long)) {
          post_avg_gaps <- placebo_outputs$gaps_long |>
            dplyr::filter(Country != "Average") |>
            dplyr::filter(year >= post_period[1], year <= post_period[2]) |>
            dplyr::group_by(Country) |>
            dplyr::summarise(avg_gap_post = mean(gap, na.rm = TRUE), .groups = "drop") |>
            dplyr::mutate(
              is_treated = Country == treatment_identifier,
              is_spain = Country == "Spain"
            )
          readr::write_csv(post_avg_gaps, file.path(pool_tables, "placebo_post_avg_gaps.csv"))
          save_table_jpg(post_avg_gaps, file.path(pool_tables, "placebo_post_avg_gaps.jpg"), title = "Average post-treatment gaps (all units)")

          p_gap <- NA_real_
          if (any(post_avg_gaps$is_treated)) {
            treated_gap <- post_avg_gaps$avg_gap_post[post_avg_gaps$is_treated]
            placebo_gaps <- post_avg_gaps$avg_gap_post[!post_avg_gaps$is_treated]
            p_gap <- rank_pvalue(treated_gap, placebo_gaps, direction = "greater")
            p_gap_df <- tibble::tibble(metric = "avg_post_gap", p_value = p_gap, treated = treated_gap, n_placebos = length(placebo_gaps))
            readr::write_csv(p_gap_df, file.path(pool_tables, "pvalue_post_gap.csv"))
            save_table_jpg(p_gap_df, file.path(pool_tables, "pvalue_post_gap.jpg"), title = "p-value (avg post gap, ranking placebo)")
            pval_accum[[length(pval_accum) + 1]] <<- tibble::tibble(
              outcome = outcome$id,
              spec = spec$name,
              pool = pool_label,
              metric = "avg_post_gap",
              p_value_raw = p_gap
            )
          }

          if (nrow(post_avg_gaps) > 0) {
            p_gap_hist <- ggplot2::ggplot(post_avg_gaps, ggplot2::aes(x = avg_gap_post, fill = is_spain)) +
              ggplot2::geom_histogram(alpha = 0.7, bins = 30, colour = "white") +
              ggplot2::geom_vline(
                data = dplyr::filter(post_avg_gaps, is_spain),
                ggplot2::aes(xintercept = avg_gap_post),
                colour = "red",
                linewidth = 1
              ) +
              ggplot2::scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "grey60")) +
              ggplot2::labs(x = "Average gap (post period)", y = "Count", title = "Distribution of post-treatment gaps (placebos vs treated)") +
              ggplot2::theme_minimal()
            safe_ggsave(file.path(pool_plots, "placebo_post_gap_hist.png"), p_gap_hist)

            p_gap_ecdf <- ggplot2::ggplot(post_avg_gaps, ggplot2::aes(x = avg_gap_post, colour = is_spain)) +
              ggplot2::stat_ecdf(linewidth = 0.9) +
              ggplot2::scale_colour_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
              ggplot2::labs(x = "Average gap (post period)", y = "ECDF", title = "ECDF of post-treatment gaps (placebos vs treated)") +
              ggplot2::theme_minimal()
            safe_ggsave(file.path(pool_plots, "placebo_post_gap_ecdf.png"), p_gap_ecdf)

            rank_df <- post_avg_gaps |>
              dplyr::arrange(avg_gap_post) |>
              dplyr::mutate(rank = dplyr::row_number())
            p_rank <- ggplot2::ggplot(rank_df, ggplot2::aes(x = rank, y = avg_gap_post, colour = is_spain)) +
              ggplot2::geom_point(size = 2) +
              ggplot2::geom_hline(yintercept = 0, linetype = "dotted") +
              ggplot2::scale_colour_manual(values = c("TRUE" = "red", "FALSE" = "grey30")) +
              ggplot2::labs(x = "Rank (ascending)", y = "Average gap (post)", title = "Ranked post-treatment gaps", subtitle = if (!is.na(p_gap)) sprintf("Placebo p-value (greater): %.3f", p_gap) else NULL) +
              ggplot2::theme_minimal()
            safe_ggsave(file.path(pool_plots, "placebo_post_gap_rank.png"), p_rank)
          }
        }

        if (!is.null(placebo_outputs$mspe_df)) {
          mspe_df_plot <- placebo_outputs$mspe_df |>
            dplyr::mutate(is_treated = Country == treatment_identifier) |>
            dplyr::filter(is.finite(mspe_ratio))
          readr::write_csv(mspe_df_plot, file.path(pool_tables, "placebo_mspe_ratios.csv"))
          if (nrow(mspe_df_plot) > 0) {
            treated_mspe_ratio <- mspe_df_plot$mspe_ratio[mspe_df_plot$is_treated]
            placebo_mspe_ratio <- mspe_df_plot$mspe_ratio[!mspe_df_plot$is_treated]
            p_mspe <- rank_pvalue(treated_mspe_ratio, placebo_mspe_ratio, direction = "greater")
            p_mspe_df <- tibble::tibble(metric = "mspe_ratio", p_value = p_mspe, treated = treated_mspe_ratio, n_placebos = length(placebo_mspe_ratio))
            readr::write_csv(p_mspe_df, file.path(pool_tables, "pvalue_mspe_ratio.csv"))
            save_table_jpg(p_mspe_df, file.path(pool_tables, "pvalue_mspe_ratio.jpg"), title = "p-value (post/pre MSPE ratio)")
            pval_accum[[length(pval_accum) + 1]] <<- tibble::tibble(
              outcome = outcome$id,
              spec = spec$name,
              pool = pool_label,
              metric = "mspe_ratio",
              p_value_raw = p_mspe
            )
            p_mspe_hist <- ggplot2::ggplot(mspe_df_plot, ggplot2::aes(x = mspe_ratio, fill = is_treated)) +
              ggplot2::geom_histogram(alpha = 0.7, bins = 30, colour = "white") +
              ggplot2::geom_vline(
                data = dplyr::filter(mspe_df_plot, is_treated),
                ggplot2::aes(xintercept = mspe_ratio),
                colour = "red",
                linewidth = 1
              ) +
              ggplot2::scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "grey60")) +
              ggplot2::labs(x = "Post/Pre MSPE ratio", y = "Count", title = "Distribution of post/pre MSPE ratios") +
              ggplot2::theme_minimal()
            safe_ggsave(file.path(pool_plots, "placebo_mspe_ratio_hist.png"), p_mspe_hist)
            p_mspe_ecdf <- ggplot2::ggplot(mspe_df_plot, ggplot2::aes(x = mspe_ratio, colour = is_treated)) +
              ggplot2::stat_ecdf(linewidth = 0.9) +
              ggplot2::scale_colour_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
              ggplot2::labs(x = "Post/Pre MSPE ratio", y = "ECDF", title = "ECDF of post/pre MSPE ratios") +
              ggplot2::theme_minimal()
            safe_ggsave(file.path(pool_plots, "placebo_mspe_ratio_ecdf.png"), p_mspe_ecdf)
            rank_mspe <- mspe_df_plot |>
              dplyr::arrange(mspe_ratio) |>
              dplyr::mutate(rank = dplyr::row_number())
            p_mspe_rank <- ggplot2::ggplot(rank_mspe, ggplot2::aes(x = rank, y = mspe_ratio, colour = is_treated)) +
              ggplot2::geom_point(size = 2) +
              ggplot2::scale_colour_manual(values = c("TRUE" = "red", "FALSE" = "grey30")) +
              ggplot2::labs(x = "Rank (ascending)", y = "Post/Pre MSPE ratio", title = "Ranked post/pre MSPE ratios", subtitle = if (!is.na(p_mspe)) sprintf("Placebo p-value (greater): %.3f", p_mspe) else NULL) +
              ggplot2::theme_minimal()
            safe_ggsave(file.path(pool_plots, "placebo_mspe_ratio_rank.png"), p_mspe_rank)
          }
        }

        # Resumen de RMSPE y tendencias pre
        if (!is.null(placebo_outputs$gaps_long)) {
          # Bandas placebo en el pre
          pre_band <- placebo_outputs$placebo_only |>
            dplyr::filter(year >= pre_period[1], year <= pre_period[2]) |>
            dplyr::group_by(year) |>
            dplyr::summarise(
              q05 = stats::quantile(gap, 0.05, na.rm = TRUE),
              q50 = stats::quantile(gap, 0.50, na.rm = TRUE),
              q95 = stats::quantile(gap, 0.95, na.rm = TRUE),
              .groups = "drop"
            )
          if (nrow(pre_band) > 0 && !is.null(gaps_long_main)) {
            readr::write_csv(pre_band, file.path(pool_tables, "pre_placebo_gap_quantiles.csv"))
            save_table_jpg(pre_band, file.path(pool_tables, "pre_placebo_gap_quantiles.jpg"), title = "Pre-treatment placebo bands (gaps)")
            treated_pre <- gaps_long_main |>
              dplyr::filter(year >= pre_period[1], year <= pre_period[2], Country == treatment_identifier)
            if (nrow(treated_pre) > 0) {
              pre_band$year_num <- as.integer(round(pre_band$year))
              treated_pre$year_num <- as.integer(round(treated_pre$year))
              p_pre <- ggplot2::ggplot() +
                ggplot2::geom_ribbon(
                  data = pre_band,
                  ggplot2::aes(x = year_num, ymin = q05, ymax = q95),
                  fill = "grey70",
                  alpha = 0.3
                ) +
                ggplot2::geom_line(
                  data = pre_band,
                  ggplot2::aes(x = year_num, y = q50),
                  linetype = "dashed",
                  colour = "black"
                ) +
                ggplot2::geom_line(
                  data = treated_pre,
                  ggplot2::aes(x = year_num, y = gap),
                  colour = "red",
                  linewidth = 1
                ) +
                ggplot2::geom_hline(yintercept = 0, linetype = "dotted") +
                ggplot2::labs(title = "Pre-treatment gaps with placebo bands", x = "Year", y = "Gap") +
                ggplot2::theme_minimal()
              safe_ggsave(file.path(pool_plots, "pre_gap_bands.png"), p_pre)
            }
          }

          # Bandas placebo en el post
          post_band <- placebo_outputs$placebo_only |>
            dplyr::filter(year >= post_period[1], year <= post_period[2]) |>
            dplyr::group_by(year) |>
            dplyr::summarise(
              q05 = stats::quantile(gap, 0.05, na.rm = TRUE),
              q50 = stats::quantile(gap, 0.50, na.rm = TRUE),
              q95 = stats::quantile(gap, 0.95, na.rm = TRUE),
              .groups = "drop"
            )
          if (nrow(post_band) > 0 && !is.null(gaps_long_main)) {
            readr::write_csv(post_band, file.path(pool_tables, "post_placebo_gap_quantiles.csv"))
            save_table_jpg(post_band, file.path(pool_tables, "post_placebo_gap_quantiles.jpg"), title = "Post-treatment placebo bands (gaps)")
            treated_post <- gaps_long_main |>
              dplyr::filter(year >= post_period[1], year <= post_period[2], Country == treatment_identifier)
            if (nrow(treated_post) > 0) {
              post_band$year_num <- as.integer(round(post_band$year))
              treated_post$year_num <- as.integer(round(treated_post$year))
              p_post <- ggplot2::ggplot() +
                ggplot2::geom_ribbon(
                  data = post_band,
                  ggplot2::aes(x = year_num, ymin = q05, ymax = q95),
                  fill = "grey70",
                  alpha = 0.3
                ) +
                ggplot2::geom_line(
                  data = post_band,
                  ggplot2::aes(x = year_num, y = q50),
                  linetype = "dashed",
                  colour = "black"
                ) +
                ggplot2::geom_line(
                  data = treated_post,
                  ggplot2::aes(x = year_num, y = gap),
                  colour = "red",
                  linewidth = 1
                ) +
                ggplot2::geom_hline(yintercept = 0, linetype = "dotted") +
                ggplot2::labs(title = "Post-treatment gaps with placebo bands", x = "Year", y = "Gap") +
                ggplot2::theme_minimal()
              safe_ggsave(file.path(pool_plots, "post_gap_bands.png"), p_post)
            }
          }

          # Test de tendencia pre (pendiente)
          slope_fun <- function(df) {
            if (nrow(df) < 2) return(NA_real_)
            fit <- try(stats::lm(gap ~ year, data = df), silent = TRUE)
            if (inherits(fit, "try-error")) return(NA_real_)
            as.numeric(coef(fit)[["year"]])
          }
          slopes <- placebo_outputs$gaps_long |>
            dplyr::filter(year >= pre_period[1], year <= pre_period[2]) |>
            dplyr::group_by(Country) |>
            dplyr::summarise(
              n_pre = dplyr::n(),
              slope = slope_fun(dplyr::cur_data_all()),
              .groups = "drop"
            )
          if (nrow(slopes) > 0) {
            readr::write_csv(slopes, file.path(pool_tables, "pre_gap_slopes.csv"))
            treated_slope <- slopes$slope[slopes$Country == treatment_identifier]
            placebo_slopes <- slopes$slope[slopes$Country != treatment_identifier & slopes$Country != "Average"]
            p_slope <- rank_pvalue(abs(treated_slope), abs(placebo_slopes), direction = "greater")
            slope_summary <- tibble::tibble(
              metric = "pre_gap_slope",
              treated_slope = treated_slope,
              p_value = p_slope,
              n_placebos = length(placebo_slopes)
            )
            readr::write_csv(slope_summary, file.path(pool_tables, "pre_gap_slope_summary.csv"))
            save_table_jpg(slope_summary, file.path(pool_tables, "pre_gap_slope_summary.jpg"), title = "Pre-treatment gap slope (rank placebo p-value)")
          }
        }

        if (!is.null(placebo_outputs$mspe_df)) {
          # Resumen RMSPE pre/post con percentiles placebo
          mspe_df_plot <- placebo_outputs$mspe_df |>
            dplyr::mutate(is_treated = Country == treatment_identifier) |>
            dplyr::filter(is.finite(pre) | is.finite(post))
          treated_pre_mspe <- mspe_df_plot$pre[mspe_df_plot$is_treated]
          treated_post_mspe <- mspe_df_plot$post[mspe_df_plot$is_treated]
          treated_ratio <- mspe_df_plot$mspe_ratio[mspe_df_plot$is_treated]
          treated_pre_mspe_val <- if (length(treated_pre_mspe) > 0) treated_pre_mspe[1] else NA_real_
          treated_post_mspe_val <- if (length(treated_post_mspe) > 0) treated_post_mspe[1] else NA_real_
          treated_ratio_val <- if (length(treated_ratio) > 0) treated_ratio[1] else NA_real_
          placebo_ratios <- mspe_df_plot$mspe_ratio[!mspe_df_plot$is_treated]
          placebo_pre_vals <- mspe_df_plot$pre[!mspe_df_plot$is_treated]
          safe_quant <- function(x) {
            if (length(x) == 0 || all(is.na(x))) return(rep(NA_real_, 3))
            stats::quantile(x, probs = c(0.05, 0.5, 0.95), na.rm = TRUE, names = FALSE)
          }
          ratio_quant <- safe_quant(placebo_ratios)
          pre_quant <- safe_quant(placebo_pre_vals)
          rmspe_summary <- tibble::tibble(
            row = c("treated", "placebo_q05", "placebo_q50", "placebo_q95"),
            pre_mspe = c(treated_pre_mspe_val, pre_quant),
            post_mspe = c(treated_post_mspe_val, NA, NA, NA),
            mspe_ratio = c(treated_ratio_val, ratio_quant)
          )
          p_mspe_val <- rank_pvalue(treated_ratio_val, placebo_ratios, direction = "greater")
          p_mspe_summary <- tibble::tibble(
            metric = "mspe_ratio_rank",
            treated = treated_ratio_val,
            p_value = p_mspe_val,
            n_placebos = length(placebo_ratios)
          )
          readr::write_csv(rmspe_summary, file.path(pool_tables, "rmspe_summary.csv"))
          save_table_jpg(rmspe_summary, file.path(pool_tables, "rmspe_summary.jpg"), title = "RMSPE pre/post summary")
          readr::write_csv(p_mspe_summary, file.path(pool_tables, "rmspe_pvalue_summary.csv"))
          save_table_jpg(p_mspe_summary, file.path(pool_tables, "rmspe_pvalue_summary.jpg"), title = "RMSPE ratio p-value (placebo rank)")
        }

        p6 <- plot_post_pre_ratio(ppratio_df)
        safe_ggsave(file.path(pool_plots, "Figure_6_mspe_ratio.png"), p6)
        saveRDS(resplacebo, file.path(pool_tables, "resplacebo_object.rds"))

        metrics_placebo <- summarize_gap_metrics(
          gaps_long = placebo_outputs$gaps_long,
          pre_window = pre_period,
          post_window = post_period,
          treatment_identifier = treatment_identifier,
          spec_name = spec$name,
          pool_label = pool_label,
          outcome = outcome$id
        )
        if (!is.null(metrics_placebo)) {
          metrics_placebo$metric_source <- "placebo_gaps"
          metrics_accum[[length(metrics_accum) + 1]] <<- metrics_placebo
        }

        list(
          mscmt = res,
          placebo = resplacebo,
          predictor_table = predictor_table,
          donor_weights = donor_weights,
          gaps_main = gaps_long_main,
          placebo_outputs = placebo_outputs
        )
      }

      pool_results[["all"]] <- run_spec_pool(pool_label = "all", donor_filter = NULL, save_legacy_figs = spec$save_legacy && outcome$legacy)
      donor_weights_all <- pool_results[["all"]]$donor_weights
      drop_env <- tolower(Sys.getenv("SKIP_DROP_ONE", "0"))
      skip_drop_one <- drop_env %in% c("1", "true", "yes")
      donor_pos <- NULL
      if (!is.null(donor_weights_all) && nrow(donor_weights_all) > 0) {
        donor_pos <- unique(donor_weights_all$Country[donor_weights_all$Weight > 0])
      }
      if (!skip_drop_one && !is.null(donor_pos) && length(donor_pos) > 0) {
        for (donor in donor_pos) {
          drop_label <- sprintf("drop_%s", donor)
          if (!is.null(pool_results[[drop_label]])) next
          drop_countries <- setdiff(unique(panel_full$country), donor)
          pool_results[[drop_label]] <- run_spec_pool(pool_label = drop_label, donor_filter = drop_countries, save_legacy_figs = FALSE)
        }
      }

      spec_outputs[[outcome$id]][[spec$name]] <- pool_results
    }
  }

  if (length(metrics_accum) > 0) {
    metrics_df <- dplyr::bind_rows(metrics_accum)
    metrics_path <- file.path(tables_dir, "specs", "gap_metrics_summary.csv")
    ensure_dir(dirname(metrics_path))
    readr::write_csv(metrics_df, metrics_path)
    save_table_jpg(metrics_df, file.path(tables_dir, "specs", "gap_metrics_summary.jpg"), title = "Gap/MSPE summary by outcome/spec/pool")
  }

  # Resumen de p-values crudos y ajustados (familia principal: baseline/pool all)
  if (length(pval_accum) > 0) {
    pvals_df <- dplyr::bind_rows(pval_accum)
    pvals_df$p_value_adj <- NA_real_
    family_mask <- pvals_df$spec == "baseline" & pvals_df$pool == "all"
    if (any(family_mask, na.rm = TRUE)) {
      pvals_df$p_value_adj[family_mask] <- stats::p.adjust(pvals_df$p_value_raw[family_mask], method = "holm")
    }
    pvals_path <- file.path(tables_dir, "specs", "pvalues_summary.csv")
    ensure_dir(dirname(pvals_path))
    readr::write_csv(pvals_df, pvals_path)
    save_table_jpg(pvals_df, file.path(tables_dir, "specs", "pvalues_summary.jpg"), title = "P-values (raw and Holm-adjusted for baseline/all)")
  }

  specs_dir <- file.path(tables_dir, "specs")
  ensure_dir(specs_dir)

  if (length(pool_status_rows) > 0) {
    pool_status_df <- dplyr::bind_rows(pool_status_rows)
    pool_status_path <- file.path(specs_dir, "pool_status_all.csv")
    readr::write_csv(pool_status_df, pool_status_path)
    save_table_jpg(pool_status_df, file.path(specs_dir, "pool_status_all.jpg"), title = "Pool status (all outcomes/specs)")
  }
  if (length(pool_summary_rows) > 0) {
    pool_summary_df <- dplyr::bind_rows(pool_summary_rows)
    pool_summary_path <- file.path(specs_dir, "pool_summary.csv")
    readr::write_csv(pool_summary_df, pool_summary_path)
    save_table_jpg(pool_summary_df, file.path(specs_dir, "pool_summary.jpg"), title = "Pool summary by outcome/spec/pool")
  }
  if (length(donor_weights_rows) > 0) {
    donor_weights_df <- dplyr::bind_rows(donor_weights_rows)
    donor_weights_path <- file.path(specs_dir, "donor_weights_summary.csv")
    readr::write_csv(donor_weights_df, donor_weights_path)
    save_table_jpg(donor_weights_df, file.path(specs_dir, "donor_weights_summary.jpg"), title = "Donor weights (all specs/pools)")
  }
  if (length(placebo_filter_rows) > 0) {
    placebo_filter_df <- dplyr::bind_rows(placebo_filter_rows)
    placebo_filter_path <- file.path(specs_dir, "placebo_filter_summary.csv")
    readr::write_csv(placebo_filter_df, placebo_filter_path)
    save_table_jpg(placebo_filter_df, file.path(specs_dir, "placebo_filter_summary.jpg"), title = "Placebo filter summary (RMSPE cutoffs)")
  }

  list(
    specs = spec_outputs,
    predictor_tables = predictor_tables,
    donor_weights = donor_weights_summary
  )
}

