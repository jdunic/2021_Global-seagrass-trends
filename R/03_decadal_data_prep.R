# Status plots 
# ----------------- #
library(tidyverse)
library(scales)
library(patchwork)
library(here)

source(here::here('R/01_functions.R'))

# Set colour scheme for bioregional trends
bioregion_colours <- RColorBrewer::brewer.pal(7, "Dark2")
bioregion_levels  <- c("Temperate North Pacific", 
                       "Temperate North Atlantic West", 
                       "Temperate North Atlantic East",
                       "Tropical Atlantic", "Mediterranean", 
                       "Tropical Indo-Pacific", 
                       "Temperate Southern Oceans")
bioregion_key <- tibble(bioregion_colours = bioregion_colours, bioregion = bioregion_levels)

area_ts <- 
  read_csv(here::here("data_outputs/clean_area_ts_for_analysis.csv"), guess_max = 3000) %>% 
  left_join(., bioregion_key) %>% 
  mutate(bioregion = factor(bioregion, bioregion_levels))

provinces <- 
  read_csv(here::here("data_outputs/clean_area_ts_for_analysis.csv"), guess_max = 3000) %>% 
  select(study_site, PROVINCE, ECOREGION) %>% 
  rename(province = PROVINCE, ecoregion = ECOREGION)

decadal_rates <- 
  read_csv(here::here("data_outputs/decadal_rates.csv"), guess_max = 2000)

decade_levels <- c('pre-1940s', as.character(seq(1940, 2010, by = 10)))


decadal_df <- 
  decadal_rates %>% 
  mutate(max_area_km2 = convertr::convert(max_area, "ha", "km2")) %>% 
  mutate(study_site = paste(studyid, site, sep = "_")) %>% 
  mutate(fractional_change = exp(log_area2_ci_50) / (max_area)) %>% 
  mutate(decade_cat = ifelse(as.numeric(decade) < 1940, "pre-1940s", decade)) %>%
  mutate(decade_cat = factor(decade_cat, levels = decade_levels)) %>% 
  mutate(bioregion = factor(bioregion, bioregion_levels)) %>% 
  mutate(size_bins_ha = case_when(max_area <= 10 ~ "< 10 ha", 
                              (max_area > 10) & (max_area <= 100) ~ "10 - 100 ha", 
                              (max_area > 100) & (max_area <= 1000) ~ "100 - 1000 ha", 
                              (max_area > 1000) ~ "> 1000 ha")) %>% 
  mutate(size_bins_ha = 
           factor(size_bins_ha, levels = c("< 10 ha", "10 - 100 ha", 
                                        "100 - 1000 ha", "> 1000 ha"))) %>% 
  mutate(size_bins_km2 = case_when(max_area_km2 <= 1 ~ "< 1 km2", 
                              (max_area_km2 > 1) & (max_area_km2 <= 10) ~ "1 - 10 km2", 
                              (max_area_km2 > 10) & (max_area_km2 <= 100) ~ "10 - 100 km2", 
                              (max_area_km2 > 100) & (max_area_km2 <= 1000) ~ "100 - 1000 km2", 
                              (max_area_km2 > 1000) ~ "> 1000 km2")) %>% 
  mutate(size_bins_km2 = 
           factor(size_bins_km2, levels = c("< 1 km2", "1 - 10 km2", "10 - 100 km2", 
                                        "100 - 1000 km2", "> 1000 km2"))) %>% 
  mutate(lh = factor(lh, levels = c('C', 'O', 'P', "other"))) %>% 
  left_join(., provinces) %>% 
  distinct(study_site, decade, .keep_all = TRUE)

decadal_n_df <- 
  decadal_df %>% 
  group_by(decade) %>% 
  summarise(study_n = n_distinct(studyid), site_n = n_distinct(study_site))

bioregion_decadal_n_df <- 
  decadal_df %>% 
  group_by(bioregion, decade) %>% 
  summarise(study_n = n_distinct(studyid), site_n = n_distinct(study_site))

# Prepare data for status plots and net change
decadal_df <- mutate(decadal_df, frac_of_max = exp(log_area1_ci_50) / max_plus_0.1)
decadal_df_1930 <- decadal_df %>% filter(decade > 1920)

# Remove extreme values for status plot visualisation
decadal_df_ne <- decadal_df %>% filter(frac_of_max < 1.2)
decadal_df_ne_1930 <- 
  decadal_df_ne %>% 
  filter(decade > 1920) %>% 
  mutate(depth_cat = factor(depth_cat, levels = c('intertidal', 'subtidal', 'mixed', 'unspecified'))) %>% 
  mutate(waycott_study = ifelse(str_detect(studyid, "W\\d"), TRUE, FALSE))