setwd(here::here())

# Load functions
source(here::here('R/01_functions.R'))

# Reconstruct time series and prepare data for visualisation
source(here::here('R/02_analysis-data-prep.R'))
source(here::here('R/02_gam-interpolation.R'))
source(here::here('R/03_decadal_data_prep.R'))

# Produce figures and tables
source(here::here('R/04_bioregional_GAM_trends.R'))
source(here::here('R/04_descriptive_stats.R'))
source(here::here('R/04_driver-attributions.R'))
source(here::here('R/04_net-change.R'))
source(here::here('R/04_site-map.R'))
source(here::here('R/04_status_plots.R'))


# Produce results files
rmarkdown::render(here::here('R/05_Figures.Rmd'), output_dir = here::here("results_summaries"))
rmarkdown::render(here::here('R/05_Tables.Rmd'), output_dir = here::here("results_summaries"))

# Supplementary results
rmarkdown::render(here::here('R/06_Supplementary_figures.Rmd'), output_dir = here::here("results_summaries"))
rmarkdown::render(here::here('R/06_Supplementary_tables.Rmd'), output_dir = here::here("results_summaries"))
