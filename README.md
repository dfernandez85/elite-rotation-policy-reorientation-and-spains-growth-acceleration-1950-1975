# Proyecto Franquismo - Replication Package

Este repositorio reproduce los resultados econometricos del manuscrito y su apendice mediante un unico pipeline de control sintetico. La version actual del proyecto esta pensada como paquete de replicacion para journal: genera solo las tablas y figuras finales citadas en el paper y el appendix.

## Inicio rapido
```r
install.packages("renv")  # solo si renv no esta instalado
renv::restore()
source("main.R")
```

O desde terminal:
```powershell
Rscript main.R
```

En el primer arranque hace falta conexion a internet: si `data/raw/` no contiene los insumos crudos fijados, `main.R` los descarga automaticamente y actualiza `data/raw/checksums.txt`.

## Estructura
```
R/
  config/         # Constantes, paquetes y parametros
  data/           # Carga y preparacion de datos
  analysis/       # Pipeline principal y exports finales del journal
  utils/          # Utilidades generales
  visualization/  # Funciones de graficos
scripts/
  fetch_data.R    # Descarga insumos y actualiza checksums
output/
  sessions/<timestamp>/
    plots/
      main/       # Figuras del manuscrito principal
      appendix/   # Figuras del apendice
    tables/
      main/       # Tablas del manuscrito principal
      appendix/   # Tablas del apendice
    run_config.yml
    session_info.txt
main.R            # Punto de entrada
```

## Dependencias
- R 4.x
- `renv` inicializado: ejecuta `renv::restore()` antes de correr `main.R`.
- Paquetes listados en `R/config/packages.R` (se cargan via `load_required_packages()`).

## Datos y reproducibilidad
- `main.R` arranca de forma autosuficiente: si faltan los insumos crudos o el checksum no es valido, descarga automaticamente las versiones fijadas y reescribe `data/raw/checksums.txt`.
- El repositorio no versiona binarios crudos ni outputs generados: los datos se recuperan automaticamente desde sus fuentes fijadas y los resultados se regeneran ejecutando `main.R`.
- Archivos requeridos:
  - `data/raw/pwt110.xlsx` (PWT 11.0, hoja `Data`), clave de checksum `pwt110`.
  - `data/raw/maddison_mpd2023.dta` (Maddison Project Database 2023), clave de checksum `maddison_mpd2023`.
- Cada archivo debe tener su md5 en `data/raw/checksums.txt`, con formato `<key> <md5sum>`.
- Si quieres forzar una descarga limpia de los insumos, usa `Rscript scripts/fetch_data.R`.
- El repositorio incluye un workflow de GitHub Actions en `/.github/workflows/reproducibility-check.yml` que verifica en Windows un `renv::restore()` limpio y la ejecucion completa de `main.R`, y sube la sesion final como artefacto.

## Como ejecutar
```r
renv::restore()        # una sola vez por entorno
source("main.R")       # desde R/RStudio
# o
Rscript main.R         # desde terminal
```

Tras `renv::restore()`, no hace falta preparar manualmente `data/raw/`: el propio pipeline lo bootstrappea si es necesario.

Cada corrida crea una sesion nueva en `output/sessions/<timestamp>/`. La salida final queda reducida a cuatro carpetas: `plots/main`, `plots/appendix`, `tables/main` y `tables/appendix`, mas `run_config.yml` y `session_info.txt`.

## Que genera
- `plots/main`: figuras del manuscrito principal.
- `plots/appendix`: figuras del apendice.
- `tables/main`: tablas del manuscrito principal.
- `tables/appendix`: tablas del apendice.
- `run_config.yml` y `session_info.txt`: metadatos de reproducibilidad de cada sesion.

## Que hace cada modulo
- `R/analysis/synthetic_control.R`: prepara el panel, ejecuta MSCMT, corre placebos, filtros de ajuste, compuerta de estabilidad y especificaciones de robustez necesarias para el paper.
- `R/analysis/journal_exports.R`: toma los objetos estimados y exporta solo las figuras y tablas canonicamente usadas en el manuscrito y el apendice; al final elimina los arboles de trabajo de `plots/` y `tables/` para dejar la sesion limpia.

La version reproducible y entregable del proyecto es `main.R`.
