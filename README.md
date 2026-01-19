# Proyecto Franquismo - Analisis Modular

Este proyecto replica el flujo de los scripts originales en una estructura modular similar a `ropke_index`. Incluye dos piezas principales:

1) Crecimiento y regimenes historicos (Maddison Project).  
2) Control sintetico para el "milagro economico" espanol (PWT 11.0, donantes Europa+Americas).

## Estructura
```
R/
  config/         # Constantes, paquetes y parametros
  data/           # Descarga y preparacion de datos
  analysis/       # Orquestadores de cada estudio
  utils/          # Utilidades generales
  visualization/  # Funciones de graficos
data/
  raw/ processed/ # Insumos/derivados locales
output/
  sessions/<timestamp>/
    plots/        # PNG generados
    tables/       # CSV (series y tablas del control sintetico)
main.R            # Punto de entrada
```

## Dependencias
- R 4.x
- `renv` inicializado: ejecuta `renv::restore()` antes de correr `main.R`.
- Paquetes listados en `R/config/packages.R` (se cargan via `load_required_packages()`).

## Datos y reproducibilidad
- No hay descargas automaticas en el pipeline. Los insumos deben existir en `data/raw/` y con hash validado.
- Archivos requeridos:
  - `data/raw/pwt110.xlsx` (PWT 11.0, GGDC, hoja "Data"), clave de checksum `pwt110`.
  - `data/raw/maddison_mpd2023.dta` (Maddison Project Database 2023), clave de checksum `maddison_mpd2023`.
- Cada archivo debe tener su md5 en `data/raw/checksums.txt`, con formato `<key> <md5sum>`. El pipeline falla si falta el archivo, si no hay hash o si el hash no coincide.
- Para obtener los datos una vez y generar hashes, usa `Rscript scripts/fetch_data.R` (usa las URLs fijas declaradas en `R/config/constants.R`). No se ejecuta dentro de `main.R`.
- Si actualizas la version de un insumo, reemplaza el archivo y su md5 en `data/raw/checksums.txt` y documenta la fuente/fecha.

## Como ejecutar
```r
renv::restore()        # una sola vez por entorno
source("main.R")       # desde R/RStudio
# o
Rscript main.R         # desde terminal
```
Los resultados se guardan en `output/sessions/<timestamp>/plots` y `.../tables`, y la informacion de la sesion (paquetes/versiones) en `output/sessions/<timestamp>/session_info.txt`.

## Que hace cada modulo
- `R/analysis/growth_regimes.R`: lee Maddison local, prepara panel de Europa Occidental, calcula crecimiento (media movil de 10 anos sin look-ahead) y un suavizado LOESS unilateral para Espana, anade regimenes historicos y genera las figuras 1-2 + una Figura 1b (LOESS).
- `R/analysis/synthetic_control.R`: usa PWT 11.0 desde XLSX local (hoja "Data"), arma panel restringido a Europa+Americas, ejecuta MSCMT (control sintetico y placebo), exporta tabla de predictores y genera figuras 3-6; usa ventana pretratamiento 1950-1959, calcula el promedio "Average" con `na.rm=TRUE`, filtra donantes por cobertura pre (1950-1959) y post (1960-1975) antes de estimar, recalcula el promedio sobre el pool limpio y guarda trayectorias placebo, cuantiles (bandas), ratios MSPE, p-values, cobertura anual/pre/post (NAs por variable), pesos de donantes, donantes incluidos y el objeto RDS de los placebos.

## Scripts originales
Se conservan `Program Statistics Spanish Economic Miracle.R` y `Program Statistics Spanish Economi Miracle2 (Synth Control).R` como referencia; su logica ahora esta en los modulos descritos arriba.
