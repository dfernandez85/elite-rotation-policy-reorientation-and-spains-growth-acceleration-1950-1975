# Estructura del proyecto

```
Franquismo/
  main.R                         # Punto de entrada: ejecuta ambos analisis, guarda salidas y session_info
  README.md                      # Instrucciones y vision general
  PROJECT_STRUCTURE.md           # Este documento
  scripts/
    fetch_data.R                 # Descarga Maddison y PWT fijados y escribe checksums
  R/
    config/
      constants.R                # Rutas/base de output y claves de checksum
      packages.R                 # Paquetes requeridos y loader
      settings.R                 # Listas de paises, ventanas, parametros
    data/
      downloaders.R              # Carga de fuentes locales (Maddison, PWT) con validacion de hash
      processors.R               # Transformaciones y paneles listos para analisis
    analysis/
      growth_regimes.R           # Orquestador de crecimiento y regimenes
      synthetic_control.R        # Orquestador de control sintetico
    utils/
      helpers.R                  # Utilidades (paths, fechas, promedios, etiquetas, checksums)
    visualization/
      plots.R                    # Funciones de graficos (Figuras 1-6)
  data/
    raw/
      checksums.txt              # md5 exigido para cada insumo crudo
      maddison_mpd2023.dta       # Maddison Project Database 2023 (no se versiona en git)
      pwt110.csv                 # Penn World Table 11.0 (no se versiona en git)
    processed/                   # Derivados intermedios (reservado)
  output/
    sessions/                    # Cada corrida crea una carpeta con plots/tables y session_info
      <timestamp>/
        plots/
        tables/
        session_info.txt
  renv/                          # Entorno reproducible (no editar a mano)
  renv.lock                      # Snapshot de dependencias
  Program Statistics Spanish Economic Miracle.R                     # Script original (referencia)
  Program Statistics Spanish Economi Miracle2 (Synth Control).R     # Script original (referencia)
```

## Flujo resumido
1) `main.R` carga config, paquetes y verifica hashes.  
2) `analysis/growth_regimes.R` produce Figuras 1-2 y CSV de series.  
3) `analysis/synthetic_control.R` ejecuta MSCMT, placebos y produce Figuras 3-6 + tablas.  
4) Se registra `sessionInfo()` en `output/sessions/<timestamp>/session_info.txt`.
