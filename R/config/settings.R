growth_countries <- c(
  "Austria", "Belgium", "Denmark", "Finland", "France", "Germany",
  "Greece", "Italy", "Netherlands", "Norway", "Portugal", "Spain",
  "Sweden", "Switzerland", "United Kingdom"
)

start_year_growth <- 1868 # anos disponibles para calculo (incluye sexenio)
plot_start_year <- 1875 # ano minimo mostrado en las figuras
moving_average_window <- 10
loess_window_years <- 25 # ventana hacia atras en anos para LOESS unilateral
loess_span <- 0.9 # fraccion de puntos dentro de la ventana usada por loess
loess_min_points <- 5 # minimo de puntos para arrancar el suavizado

regime_windows <- data.frame(
  start = as.Date(c(
    "1875-01-01", "1923-01-01",
    "1930-01-01", "1935-01-01", "1939-01-01", "1975-01-01"
  )),
  end = as.Date(c(
    "1923-12-31", "1930-12-31",
    "1935-12-31", "1939-12-31", "1975-12-31", "2022-12-31"
  )),
  regimes = factor(c(
    "Bourbon Restoration",
    "Primo Rivera Dictatorship",
    "Second Republic",
    "Civil War",
    "Francoism",
    "Democracy"
  ))
)
regime_levels <- as.character(regime_windows$regimes)

synthetic_countries <- c(
  # Europa
  "Albania", "Austria", "Belgium", "Cyprus", "Denmark", "Finland",
  "France", "Germany", "Iceland", "Ireland", "Italy", "Luxembourg",
  "Netherlands", "Norway", "Portugal", "Spain", "Sweden",
  "Switzerland", "Turkey", "United Kingdom",
  # America
  "Argentina", "Bolivia", "Brazil", "Canada", "Chile", "Colombia",
  "Costa Rica", "Ecuador", "El Salvador", "Guatemala", "Honduras",
  "Mexico", "Nicaragua", "Panama", "Peru", "Trinidad and Tobago",
  "United States of America", "Uruguay", "Venezuela"
)

region_codes <- NULL # se generan codigos numericos automaticamente

pre_treatment_years <- c(1950, 1959) # Incluye 1959 en el periodo pre
post_cutoff_year <- 1976 # Regla ex ante: descarta placebos cuyo RMSPE pre supere este multiple del RMSPE pre de Espana
# (criterio de Abadie et al. 2010). Si es NA o <= 0, no se filtra.
placebo_mspe_ratio_cutoff <- 10
