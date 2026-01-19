source("R/config/constants.R")
source("R/config/packages.R")
source("R/config/settings.R")
source("R/utils/helpers.R")
source("R/data/downloaders.R")
source("R/data/processors.R")
source("R/visualization/plots.R")
load_required_packages()
raw <- load_pwt_data()
panel <- prepare_pwt11_panel(raw)
panel_df <- as.data.frame(panel)
panel_df$country <- as.character(panel_df$country)

# replicate filters briefly to keep focus on listFromLong without exclusions
panel_df$region <- as.numeric(panel_df$region)
res <- MSCMT::listFromLong(panel_df,
  unit.variable = which(names(panel_df)=="region"),
  time.variable = which(names(panel_df)=="year"),
  unit.names.variable = which(names(panel_df)=="country")
)
str(res)
