library(tidyverse)
library(mgcv)
library(visreg)
library(here)
library(ggforce)
library(imputeTS)
library(DataGLMRepeat)
library(beepr)
# ----------------- #
# GAM factor spline
# ----------------- #
message("\n     L12: Preparing data for GAM reconstructions.")
#source(here("R/00_functions.R"))

dat <- 
  read_csv(here("data_outputs/clean_area_ts_for_analysis.csv"), guess_max = 3000) %>% 
  select(study_site, site_groupings, studyid, site, year, area, log_area, 
         max_area, max_plus_0.1, depth_cat, first_year, last_year, duration, decade_time_points,
         num_unique_years, lat, lon, dom_species, all_species, lh, att_method, lat_zone2, 
         apriori_cat, tidied_drivers, measurement_cat, waycott_study, Lat_Zone, 
         ECOREGION, PROVINCE, bioregion, decade, site_duplication, 
         tidied_drivers) %>% 
    mutate_if(is.character, as.factor) %>%  # need factor for gam
    mutate(study_site = make.names(study_site)) %>% 
    mutate(study_site = as.factor(study_site))

# Set colour scheme for bioregional trends
bioregion_colours <- RColorBrewer::brewer.pal(7, "Dark2")
bioregion_levels  <- c("Temperate North Pacific", 
                       "Temperate North Atlantic West", 
                       "Temperate North Atlantic East", 
                       "Tropical Atlantic", "Mediterranean", 
                       "Tropical Indo-Pacific", 
                       "Temperate Southern Oceans")
bioregion_key <- tibble(bioregion_colours = bioregion_colours, bioregion = bioregion_levels)

# Set parameters for GAM
k_lower <- 2 #lower value for k (too small and it breaks gam)
k_upper <- 8 #upper value for k (too big and computations slow down)
k_diff <- 2 #how much lower than n do you want k to be? 

# Create dataframe of site by site sample sizes and k values
dreg <- 
  dat %>% 
    group_by(bioregion, study_site) %>%
    summarize(n = n()) %>%
    group_by(bioregion, study_site) %>%
    mutate(k = min(max(k_lower, n-k_diff), k_upper)) %>% 
    ungroup()

# Make site by site smooths for the GAM, using k from dataframe above
funform <- function(i, dregions){
  thisreg <- dregions$study_site[i]
  paste0("s(year, by = as.numeric(study_site == ",shQuote(thisreg),"), m = 1, k = ", dregions$k[i], ",bs = 'tp')")
}

# List of formulae
get_bioregion_form <- function(bioregion_df) {
  #bioregion_df <- filter(full_df, bioregion == region)
  x <- lapply(1:nrow(bioregion_df), funform, dregions = bioregion_df)
  form <- 
    as.formula(paste("log_area ~ s(year, k = 5, m = 2, bs = 'tp') +",
                      do.call("paste", c(x, sep = " + ")),
                      " + s(site_groupings, bs = 're')"))
  return(form)
}

# ------------------------------
# Commented out to save time
# ------------------------------
# # Fit the models
# ------------------------------
# message("\n     L116: Fitting and saving GAM objects")
# #takes about 1 minute
# north_pacific_mod <- 
#   bam(get_bioregion_form(filter(dreg, bioregion == "Temperate North Pacific")), 
#       data = filter(dat, bioregion == "Temperate North Pacific") %>% droplevels())

# atlantic_west_mod <- 
#   bam(get_bioregion_form(filter(dreg, bioregion == "Temperate North Atlantic West")), 
#       data = filter(dat, bioregion == "Temperate North Atlantic West"))

# atlantic_east_mod <- 
#   bam(get_bioregion_form(filter(dreg, bioregion == "Temperate North Atlantic East")), 
#       data = filter(dat, bioregion == "Temperate North Atlantic East"))

# atlantic_tropical_mod <- 
#   bam(get_bioregion_form(filter(dreg, bioregion == "Tropical Atlantic")), 
#       data = filter(dat, bioregion == "Tropical Atlantic"))

# indo_pacific_mod <- 
#   bam(get_bioregion_form(filter(dreg, bioregion == "Tropical Indo-Pacific")), 
#       data = filter(dat, bioregion == "Tropical Indo-Pacific"))

# southern_ocean_mod <- 
#   bam(get_bioregion_form(filter(dreg, bioregion == "Temperate Southern Oceans")), 
#       data = filter(dat, bioregion == "Temperate Southern Oceans"))

# # Troublesome bioregion
# mediterranean_mod <-
#   bam(formula = log_area ~ s(year, k = 3, bs = "tp", m = 2) +
#                            s(year, by = study_site, k = 3, bs = "tp", m = 1) +
#                            s(study_site, bs = "re"), 
#       data = filter(dat, num_unique_years >=3, bioregion == "Mediterranean"), 
#       method = 'REML')
# beep()

# saveRDS(north_pacific_mod, file = here::here("data_outputs/GAM_objects/Temperate_North_Pacific_varying-K_v2.RDS"))
# saveRDS(atlantic_west_mod, file = here::here("data_outputs/GAM_objects/Temperate_North_Atlantic_West_varying-K_v2.RDS"))
# saveRDS(atlantic_east_mod, file = here::here("data_outputs/GAM_objects/Temperate_North_Atlantic_East_varying-K_v2.RDS"))
# saveRDS(mediterranean_mod, file = here::here("data_outputs/GAM_objects/Mediterranean_varying-K_v2.RDS"))
# saveRDS(atlantic_tropical_mod, file = here::here("data_outputs/GAM_objects/Tropical_Atlantic_varying-K_v2.RDS"))
# saveRDS(indo_pacific_mod, file = here::here("data_outputs/GAM_objects/Tropical_Indo_Pacific_varying-K_v2.RDS"))
# saveRDS(southern_ocean_mod, file = here::here("data_outputs/GAM_objects/Southern_Ocean_varying-K_v2.RDS"))

# Load GAM objects
message("\n     L117: Loading saved GAM objects")
north_pacific_mod <- readRDS(file = here("data_outputs/GAM_objects/Temperate_North_Pacific_varying-K_v2.RDS"))
atlantic_west_mod <- readRDS(file = here("data_outputs/GAM_objects/Temperate_North_Atlantic_West_varying-K_v2.RDS"))
atlantic_east_mod <- readRDS(file = here("data_outputs/GAM_objects/Temperate_North_Atlantic_East_varying-K_v2.RDS"))
mediterranean_mod <- readRDS(file = here("data_outputs/GAM_objects/Mediterranean_varying-K_v2.RDS"))
atlantic_tropical_mod <- readRDS(file = here("data_outputs/GAM_objects/Tropical_Atlantic_varying-K_v2.RDS"))
indo_pacific_mod <- readRDS(file = here("data_outputs/GAM_objects/Tropical_Indo_Pacific_varying-K_v2.RDS"))
southern_ocean_mod <- readRDS(file = here("data_outputs/GAM_objects/Southern_Ocean_varying-K_v2.RDS"))

message("\n     L126: Getting Mediterranean 2 time-point data")
# Mediterranean gam couldn't handle all the sites with only 2 sampling points
med_2pt_list <-
  filter(dat, num_unique_years <3, bioregion == "Mediterranean") %>% 
    split(., .$study_site) %>% 
    rlist::list.clean(., fun = function(x) nrow(x) == 0)

mediterranean_2pts <- 
  lapply(med_2pt_list, FUN = function(study_site_df) {
    # Linear interpolation
    imputation_df <- 
      complete(year = full_seq(year, period = 1), study_site_df) %>% 
        select(year, study_site, decade, log_area) %>% 
        fill(study_site)
    imputations <- 
      do_na_interpolation(imputation_df, col = "log_area", new_col = "log_area_impute") %>% 
      mutate(area_impute = exp(log_area_impute),
             decade = floor(year / 10) * 10)
    return(imputations)
  }) %>% 
  rlist::list.clean() %>% 
  bind_rows()

med_2pt_rates <- 
  mediterranean_2pts %>% 
  group_by(study_site, decade) %>% 
  slice(1, n()) %>% 
  mutate(timepoint = c('first', 'last')) %>% 
  ungroup() %>% 
  mutate(year = ifelse(timepoint == "first", year, year + 0.99)) %>% 
  pivot_wider(., id_cols = c(study_site, decade), names_from = timepoint, values_from = c(year,  area_impute, log_area_impute))  %>% 
  group_by(study_site, decade) %>% 
  summarise(year_first = unique(year_first), 
    year_last = unique(year_last), 
    time_diff = year_last - year_first,
    log_area1 = log_area_impute_first,
    log_area2 = log_area_impute_last,
    area_1 = exp(log_area1), 
    area_2 = exp(log_area2),
    lr_trend = (log_area_impute_last - log_area_impute_first) / (year_last - year_first),
    trend = 100 * (exp((log_area_impute_last - log_area_impute_first) / (year_last - year_first)) - 1), 
    area_change = exp(log_area_impute_last) - exp(log_area_impute_first), 
    area_change2 = area_impute_last - area_impute_first) %>% 
  ungroup() %>% 
  select(-c(time_diff, lr_trend, area_change2))


# Plot all site fits to check
# ------------------------------------------------------------------------------
plot_all_sites <- function(data, model_obj, bioregion) {
  newdata <-   
    filter(data, bioregion == !!bioregion) %>% 
    select(year, study_site, site_groupings, decade) %>% 
    complete(year = full_seq(year, period = 1)) %>% 
    fill(study_site, site_groupings, decade)
  preds <- get_yearly_preds(model_obj = model_obj, newdata = newdata)
  raw_data_subset <- 
  filter(data, bioregion == !!bioregion) %>% 
  select(study_site, site_groupings, year, area, log_area, max_area, max_plus_0.1)
  pred_out <- left_join(preds, raw_data_subset)
  return(pred_out)
}

get_site_plots <- function(pred_out) {
  gout <- pred_out %>% 
    group_by(study_site) %>% 
    with_groups(., {
      ggplot() + 
        aes(x = year, y = `log_area_50%`) + 
        geom_line(size = 1) + 
        geom_line(aes(x = year, y = log_area), colour = 'red', size = 1.3)+
        geom_point(aes(x = year, y = log_area), colour = 'red', size = 3)+
        geom_ribbon(aes(ymin = `log_area_2.5%`, ymax = `log_area_97.5%`), alpha = 0.5, color = NA) +
        scale_color_identity() +
        scale_fill_identity() +
        mytheme(base_size = 16) +
        theme(axis.title.y = element_text(size = 12))+
        ylab(paste0("Log Area (ha)\n", study_site)) +
        guides(fill = FALSE, colour = FALSE)
    })
}

# # Temperate North Pacific
# pred_out <- 
#   plot_all_sites(data = dat, model_obj = north_pacific_mod, bioregion = "Temperate North Pacific")
# gout <- get_site_plots(pred_out)
# gout[[1]]
# ml <- marrangeGrob(gout, nrow=4, ncol=2)
# ggsave(path = "../data_outputs/site_level_test", filename = "Temperate_North_Pacific_varying-K_v2.pdf", ml, width = 8, height = 11)

# # Temperate North Atlantic West
# pred_out <- 
#   plot_all_sites(data = dat, model_obj = atlantic_west_mod, bioregion = "Temperate North Atlantic West")
# gout <- get_site_plots(pred_out)
# gout[[1]]
# ml <- marrangeGrob(gout, nrow=4, ncol=2)
# ggsave(path = "../data_outputs/site_level_test", filename = "Temperate_Northern_Atlantic_West_varying-K_v2.pdf", ml, width = 8, height = 11)

# # Temperate North Atlantic East
# pred_out <- 
#   plot_all_sites(data = dat, model_obj = atlantic_east_mod, bioregion = "Temperate North Atlantic East")
# gout <- get_site_plots(pred_out)
# gout[[1]]
# ml <- marrangeGrob(gout, nrow=4, ncol=2)
# ggsave(path = "../data_outputs/site_level_test", filename = "Temperate_Northern_Atlantic_East_varying-K_v2.pdf", ml, width = 8, height = 11)

# # Tropical Atlantic
# pred_out <- 
#   plot_all_sites(data = dat, model_obj = atlantic_tropical_mod, bioregion = "Tropical Atlantic")
# gout <- get_site_plots(pred_out)
# gout[[1]]
# ml <- marrangeGrob(gout, nrow=4, ncol=2)
# ggsave(path = "../data_outputs/site_level_test", filename = "Tropical_Atlantic_varying-K_v2.pdf", ml, width = 8, height = 11)

# # Mediterranean
# pred_out <- 
#   plot_all_sites(data = filter(dat, num_unique_years > 2), model_obj = mediterranean_mod, bioregion = "Mediterranean")
# gout <- get_site_plots(pred_out)
# gout[[1]]
# ml <- marrangeGrob(gout, nrow=4, ncol=2)
# ggsave(path = "../data_outputs/site_level_test", filename = "Mediterranean_varying-K_v2.pdf", ml, width = 8, height = 11)

# # Tropical Indo-Pacific
# pred_out <- 
#   plot_all_sites(data = dat, model_obj = indo_pacific_mod, bioregion = "Tropical Indo-Pacific")
# gout <- get_site_plots(pred_out)
# gout[[1]]
# ml <- marrangeGrob(gout, nrow=4, ncol=2)
# ggsave(path = "../data_outputs/site_level_test", filename = "Tropical_Indo_Pacific_varying-K_v2.pdf", ml, width = 8, height = 11)

# # Temperate Southern Oceans
# pred_out <- 
#   plot_all_sites(data = dat, model_obj = southern_ocean_mod, bioregion = "Temperate Southern Oceans")
# gout <- get_site_plots(pred_out)
# gout[[1]]
# ml <- marrangeGrob(gout, nrow=4, ncol=2)
# ggsave(path = "../data_outputs/site_level_test", filename = "Southern_Ocean_varying-K_v2.pdf", ml, width = 8, height = 11)

# ------------------------------------------------------------------------------
# Get bioregional trends
# ------------------------------------------------------------------------------
get_bioregion_fits <- function(gam_obj, bioregion) {
  bioregion_plot <- plot.gam(gam_obj, select = 1, scale = 0)
  bioregion_fits <- 
    tibble(bioregion = bioregion,
           x = bioregion_plot[[1]]$x, 
           fit = bioregion_plot[[1]]$fit[, 1], 
           se = bioregion_plot[[1]]$se) %>% 
    mutate(lwr = fit - se, upr = fit + se)
  return(bioregion_fits)
}

message("\n     L278: Getting and saving bioregional fits.")
north_pacific  <- get_bioregion_fits(north_pacific_mod, bioregion = "Temperate North Pacific")
atlantic_west  <- get_bioregion_fits(atlantic_west_mod, bioregion = "Temperate North Atlantic West")
atlantic_east  <- get_bioregion_fits(atlantic_east_mod, bioregion = "Temperate North Atlantic East")
atlantic_trop  <- get_bioregion_fits(atlantic_tropical_mod, bioregion = "Tropical Atlantic")
indo_pacific   <- get_bioregion_fits(indo_pacific_mod, bioregion = "Tropical Indo-Pacific")
southern_ocean <- get_bioregion_fits(southern_ocean_mod, bioregion = "Temperate Southern Oceans")
mediterranean  <- get_bioregion_fits(mediterranean_mod, bioregion = "Mediterranean")


all_bioregions <- 
  left_join(
    bind_rows(north_pacific, atlantic_west, atlantic_east, atlantic_trop, indo_pacific, southern_ocean, mediterranean), 
    bioregion_key
  )
saveRDS(all_bioregions, file = here::here("data_outputs/GAM_bioregion_fits.RDS"))
all_bioregions <-readRDS(file = here("data_outputs/GAM_bioregion_fits.RDS"))


# -------------------------------
get_decadal_rates <- function(data = dat, gam_obj, bioregion) {
  n_draws = 1000
  ci_probs = c(0.025, 0.5, 0.975)
  raw_data <- 
    filter(data, bioregion == !!bioregion) %>% 
    select(year, study_site, site_groupings)
  newdata <- 
    raw_data %>% 
      group_by(study_site) %>%
      complete(year = full_seq(year, period = 1)) %>% 
      fill(-year) %>% 
      ungroup() %>%
      mutate(decade = floor(year / 10) * 10) %>% 
      group_by(study_site, decade) %>% 
      slice(1, n()) %>% 
      mutate(timepoint = c('first', 'last')) %>%
      ungroup() %>% 
      mutate(year = ifelse(timepoint == "first", year, year +  0.99))
  # linear predictors
  lp1 <- predict(gam_obj, newdata = newdata, type = "lpmatrix")
  # n_draws of the parameters from covariance matrix
  br <- rmvn(n_draws, coef(gam_obj), gam_obj$Vp)
#
  pred_area  <- matrix(NA, nrow = nrow(newdata), ncol = n_draws)
  for (i in 1:n_draws) {
    pred_area[, i] <- lp1 %*% br[i, ]
  }
  pred_area_df <- bind_cols(newdata, as_tibble(pred_area))
#
  pred_out <- 
    pred_area_df %>% 
      pivot_longer(., V1:V1000, names_to = "sample", values_to = "log_area") %>%
      pivot_wider(., names_from = timepoint, values_from = c(year, log_area)) %>% 
      mutate(
        area1 = exp(log_area_first),
        area2 = exp(log_area_last), 
        trend = 100 * (exp((log_area_last - log_area_first) / (year_last - year_first)) - 1), 
        area_change = exp(log_area_last) - exp(log_area_first)) %>% 
      group_by(study_site, decade) %>%
      summarise(year_first = unique(year_first), 
                year_last = unique(year_last),
                log_area1_ci_2.5  = quantile(log_area_first, probs = 0.025), 
                log_area1_ci_50   = quantile(log_area_first, probs = 0.50), 
                log_area1_ci_97.5 = quantile(log_area_first, probs = 0.975), 
                log_area2_ci_2.5  = quantile(log_area_last, probs = 0.025), 
                log_area2_ci_50   = quantile(log_area_last, probs = 0.50), 
                log_area2_ci_97.5 = quantile(log_area_last, probs = 0.975), 
                trend_ci_2.5  = quantile(trend, probs = 0.025), 
                trend_ci_50   = quantile(trend, probs = 0.50), 
                trend_ci_97.5 = quantile(trend, probs = 0.975), 
                area_change_ci_2.5  = quantile(area_change, probs = 0.025), 
                area_change_ci_50   = quantile(area_change, probs = 0.50), 
                area_change_ci_97.5 = quantile(area_change, probs = 0.975)
                )  %>% 
      ungroup() %>% 
      mutate(est_area1 = exp(log_area1_ci_50), 
             est_area2 = exp(log_area2_ci_50))
}

message("\n     L357: Getting decadal rates for each bioregion.")
north_pacific_rates  <- 
  get_decadal_rates(data = dat, gam_obj = north_pacific_mod, 
                    bioregion = "Temperate North Pacific")
atlantic_west_rates  <- 
  get_decadal_rates(data = dat, gam_obj = atlantic_west_mod, 
                    bioregion = "Temperate North Atlantic West")
atlantic_east_rates  <- 
  get_decadal_rates(data = dat, gam_obj = atlantic_east_mod, 
                    bioregion = "Temperate North Atlantic East")
atlantic_trop_rates  <- 
  get_decadal_rates(data = dat, gam_obj = atlantic_tropical_mod, 
                    bioregion = "Tropical Atlantic")
indo_pacific_rates   <- 
  get_decadal_rates(data = dat, gam_obj = indo_pacific_mod, 
                    bioregion = "Tropical Indo-Pacific")
southern_ocean_rates <- 
  get_decadal_rates(data = dat, gam_obj = southern_ocean_mod, 
                    bioregion = "Temperate Southern Oceans")
mediterranean_rates  <- 
  get_decadal_rates(data = filter(dat, num_unique_years >=3), 
                    gam_obj = mediterranean_mod, 
                    bioregion = "Mediterranean")

all_rates <- 
  bind_rows(north_pacific_rates, atlantic_west_rates, atlantic_east_rates, atlantic_trop_rates, indo_pacific_rates, southern_ocean_rates, mediterranean_rates) %>% 
  # Add 2 point mediterranean rates. Missing CIs right now
  bind_rows(., 
      med_2pt_rates %>% 
      rename(log_area1_ci_50 = log_area1, 
             log_area2_ci_50 = log_area2, 
             trend_ci_50 = trend, 
             area_change_ci_50 = area_change, 
             est_area1 = area_1, 
             est_area2 = area_2)
  )
all_rates <- 
  left_join(all_rates, select(dat, study_site, site_groupings, studyid, site, year, area, log_area, max_area, max_plus_0.1, apriori_cat, depth_cat, lh, lat_zone2, duration, lat, lon, dom_species, bioregion, site_duplication, decade_time_points, measurement_cat, tidied_drivers, att_method)) %>% 
  left_join(., bioregion_key) %>% 
  mutate(bioregion = factor(bioregion, levels = bioregion_levels)) %>% 
  distinct(study_site, decade, .keep_all = TRUE)

# Save estimated rates of change for further analyses
write_csv(all_rates, path = here::here("data_outputs/decadal_rates.csv"))
