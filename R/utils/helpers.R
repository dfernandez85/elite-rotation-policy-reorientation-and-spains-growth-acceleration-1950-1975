ensure_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  invisible(path)
}

read_checksum_file <- function(path = CHECKSUM_FILE) {
  if (!file.exists(path)) {
    stop(sprintf("Checksum file not found at %s. Create it or run scripts/fetch_data.R to populate it.", path))
  }

  lines <- readLines(path, warn = FALSE)
  lines <- trimws(lines)
  lines <- lines[lines != "" & !grepl("^#", lines)]

  if (length(lines) == 0) {
    stop(sprintf("Checksum file %s is empty. Add lines like '<key> <md5sum>'.", path))
  }

  parts <- strsplit(lines, "\\s+")
  keys <- vapply(parts, function(x) x[1], character(1))
  hashes <- vapply(parts, function(x) if (length(x) >= 2) x[2] else NA_character_, character(1))
  names(hashes) <- keys
  hashes
}

is_placeholder_hash <- function(hash_value) {
  if (is.null(hash_value) || length(hash_value) == 0) return(TRUE)
  hv <- trimws(hash_value[1])
  hv == "" || is.na(hv) || grepl("FILL", hv, ignore.case = TRUE) || hv == "ADD_MD5"
}

validate_file_hash <- function(path, checksum_key, expected_hash = NULL) {
  if (!file.exists(path)) {
    stop(sprintf(
      "Missing data file: %s. Place the pinned file and its md5 in %s (see scripts/fetch_data.R).",
      path, CHECKSUM_FILE
    ))
  }

  if (is.null(expected_hash)) {
    checks <- read_checksum_file()
    expected_hash <- checks[[checksum_key]]
  }

  if (is_placeholder_hash(expected_hash)) {
    stop(sprintf(
      "No valid checksum found for key '%s' in %s. Add the md5 for %s.",
      checksum_key, CHECKSUM_FILE, basename(path)
    ))
  }

  file_hash <- unname(tools::md5sum(path))

  if (!identical(tolower(file_hash), tolower(expected_hash))) {
    stop(sprintf(
      "Checksum mismatch for %s. Found %s but expected %s. Replace the file or update %s if you intentionally changed the data version.",
      basename(path), file_hash, expected_hash, CHECKSUM_FILE
    ))
  }

  invisible(path)
}

year_to_date <- function(year_numeric) {
  as.Date(zoo::as.yearmon(year_numeric))
}

apply_regime_labels <- function(dates, windows, levels) {
  labels <- rep("Unclassified", length(dates))
  for (i in seq_len(nrow(windows))) {
    in_range <- dates >= windows$start[i] & dates <= windows$end[i]
    labels[in_range] <- as.character(windows$regimes[i])
  }
  factor(labels, levels = c(levels, "Unclassified"))
}

moving_average <- function(series, window) {
  # Requiere ventanas completas: si faltan observaciones en la ventana, devuelve NA
  zoo::rollapply(
    data = series,
    width = window,
    align = "right",
    fill = NA,
    FUN = function(x) {
      if (any(is.na(x))) return(NA_real_)
      mean(x)
    }
  )
}

safe_delt <- function(x, k = 1) {
  quantmod::Delt(x, k = k)
}

one_sided_loess <- function(x, y,
                            window_years = loess_window_years,
                            span = loess_span,
                            min_points = loess_min_points) {
  x_num <- as.numeric(x)
  n <- length(y)
  smooth <- rep(NA_real_, n)
  for (i in seq_len(n)) {
    cutoff <- x_num[i] - window_years * 365.25
    idx <- seq_len(i)
    keep <- idx[!is.na(y[idx]) & !is.na(x_num[idx]) & x_num[idx] >= cutoff]
    if (length(keep) < min_points) next
    fit <- try(stats::loess(
      y[keep] ~ x_num[keep],
      span = span,
      control = stats::loess.control(surface = "direct", statistics = "exact")
    ), silent = TRUE)
    if (inherits(fit, "try-error")) next
    smooth[i] <- stats::predict(fit, newdata = x_num[i])
  }
  smooth
}

save_table_jpg <- function(df, path, title = NULL, max_rows = 40) {
  if (!is.data.frame(df)) {
    df <- as.data.frame(df)
  }
  if (nrow(df) == 0) {
    warning(sprintf("save_table_jpg skipped for %s: zero rows", path))
    return(invisible(NULL))
  }
  if (nrow(df) > max_rows) {
    df <- df[seq_len(max_rows), , drop = FALSE]
  }

  font_size <- 0.7
  tbl <- gridExtra::tableGrob(
    df,
    rows = NULL,
    theme = gridExtra::ttheme_minimal(
      core = list(fg_params = list(cex = font_size)),
      colhead = list(fg_params = list(fontface = "bold", cex = font_size))
    )
  )

  if (!is.null(title)) {
    tbl <- gridExtra::arrangeGrob(
      tbl,
      top = grid::textGrob(title, gp = grid::gpar(fontsize = 12, fontface = "bold"))
    )
  }

  # Permite mayor ancho para tablas anchas y reduce altura por fila para menos corte
  width_cm <- min(60, 3 + 2.8 * ncol(df))
  height_cm <- min(60, 2 + 0.6 * nrow(df))

  ggplot2::ggsave(
    filename = path,
    plot = tbl,
    width = width_cm,
    height = height_cm,
    units = "cm",
    dpi = 300
  )
}

