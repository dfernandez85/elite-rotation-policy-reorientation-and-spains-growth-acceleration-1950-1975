# Estructura del proyecto

```
Franquismo/
  main.R                         # Punto de entrada del replication package
  README.md                      # Instrucciones y vision general
  PROJECT_STRUCTURE.md           # Este documento
  scripts/
    fetch_data.R                 # Fuerza una descarga limpia de PWT/Maddison y reescribe checksums
  R/
    config/
      constants.R                # Rutas base de output y claves de checksum
      packages.R                 # Paquetes requeridos y loader
      settings.R                 # Paises, ventanas y parametros del pipeline
    data/
      downloaders.R              # Bootstrap de datos fijados, validacion de hash y carga local
      processors.R               # Transformaciones y paneles listos para estimar
    analysis/
      synthetic_control.R        # Pipeline de estimacion
      journal_exports.R          # Export de tablas/figuras finales del paper
    utils/
      helpers.R                  # Utilidades generales
    visualization/
      plots.R                    # Funciones de graficos del control sintetico
  data/
    raw/
      checksums.txt              # md5 exigido para cada insumo crudo
      maddison_mpd2023.dta       # Maddison Project Database 2023 (no se versiona en git)
      pwt110.xlsx                # Penn World Table 11.0 (no se versiona en git)
    processed/                   # Derivados intermedios (reservado)
  output/
    sessions/
      <timestamp>/
        plots/
          main/
          appendix/
        tables/
          main/
          appendix/
        run_config.yml
        session_info.txt
  renv/
  renv.lock
```

## Flujo resumido
1. `main.R` carga config y paquetes; si faltan los insumos crudos o sus checksums, los descarga automaticamente.
2. `R/analysis/synthetic_control.R` estima las especificaciones necesarias del manuscrito y del apendice.
3. `R/analysis/journal_exports.R` exporta unicamente las figuras y tablas finales y limpia los arboles intermedios de `plots/` y `tables/`.
4. La sesion deja `session_info.txt` y `run_config.yml` como metadatos de reproducibilidad.
