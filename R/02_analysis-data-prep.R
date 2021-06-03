# Analysis prep

library(tidyverse)
library(sf)
library(lme4)
library(metafor)
library(broom)
library(gridExtra)
library(kableExtra)
library(beepr)

setwd(here::here())
source(here::here('R/01_functions.R'))

area_ts <- 
  read_csv(here::here('data/clean_areal_extent_ts.csv'), guess_max = 3000) %>% 
  mutate(absolute_change = A2 - A1) %>% 
  mutate(apriori_cat = replace_na(apriori_cat, "none")) %>% 
  mutate(lat_zone2 = ifelse(Lat_Zone == "Polar", "Temperate", Lat_Zone)) %>% 
  mutate(lat_zone2 = factor(lat_zone2, levels = c("Temperate", "Tropical"))) %>% 
  mutate(lat_zone2 = factor(lat_zone2, levels = c("Temperate", "Tropical"))) %>% 
  mutate(depth_cat = replace_na(depth_cat, "unspecified")) %>% 
  mutate(depth_cat = factor(depth_cat, levels = c('intertidal', 'subtidal', 'mixed', 'unspecified'))) %>% 
  mutate(lh_group2 = replace_na(lh_group2, "unspecified")) %>%
  mutate(lh_group2 = factor(lh_group2, levels = c("C", "CO", "O", "P", "CP", "COP", 'unspecified'))) %>% 
  mutate(lh = ifelse((is.na(lh_group3) | lh_group3 == "mixed"), "other", lh_group3)) %>% 
  mutate(lh = factor(lh, levels = c("other", "C", "O", "P"))) %>% 
  mutate(study_site_grouping = factor(paste0(site_groupings, studyid)))

# Preparations for rates of change analysis
# Calculate 10% of maximum size ever seen to add to 0 values to allow for log
# transformations
ten_percents <- 
  area_ts %>% 
    group_by(studyid, site) %>% 
    filter(area != 0) %>% 
    summarise(min_area = min(area), 
              max_area = max(area)) %>% 
    mutate(ten_percent = 0.1 * min_area)

area_ts <- 
  left_join(area_ts, ten_percents, by = c('studyid', 'site')) %>% 
  mutate(plus_10 = area + ten_percent) %>% 
  mutate(log_area = log(plus_10)) %>% 
  rename(max_area = max_area.y) %>% 
  select(-max_area.x) %>% 
  mutate(max_plus_0.1 = max_area + ten_percent)

# Decadal analysis prep --------------------------------------------------------
unique_year_count <- 
  area_ts %>% 
    group_by(study_site) %>% 
    summarise(num_unique_years = n_distinct(year))

area_ts <- area_ts %>% mutate(decade = floor(year / 10) * 10)

# Get number of time points for each site within each decade
area_ts <- 
  area_ts %>% 
  group_by(study_site, decade) %>% 
  summarise(decade_time_points = n()) %>% 
  ungroup() %>% 
  left_join(area_ts, .) %>%
  left_join(., unique_year_count)

# Best attempt at combining global seagrass bioregions with MEOWs
area_ts <- 
  area_ts %>% 
  mutate(bioregion = 
    case_when((ECO_CODE_X %in% c(1:3, 18:27)) & (lon > -30) ~
              "Temperate North Atlantic East", 
              (ECO_CODE_X %in% c(3:11, 37:41)) & (lon <= -30) ~
              "Temperate North Atlantic West", 
              ECO_CODE_X %in% c(42, 43, 62:77, 180, 79:86, 190) |
              ECO_CODE_X %in% c(166:170) ~ # Short also includes Pacific side
              "Tropical Atlantic", 
              ECO_CODE_X %in% c(28:36, 44) ~ 
              "Mediterranean", 
              ECO_CODE_X %in% c(45:52) & (lat > 32) | # East side with lat estimate from Fig 1 Short et al.
              ECO_CODE_X %in% c(12:14, 53:61) ~  # West side
              "Temperate North Pacific", 
              ECO_CODE_X %in% c(87:150, 152:163) ~
              "Tropical Indo-Pacific", 
              ECO_CODE_X %in% c(151, 177:179, 181:188, 191:214) ~  # a bit rough around australia
              "Temperate Southern Oceans"
              )
    )

# Simplify measurement categories for analysis
area_ts <- 
  area_ts %>% 
    mutate(measurement_cat = gsub(", ground mapping", "", measurement_cat)) %>% 
    mutate(measurement_cat = gsub(", ground truthing", "", measurement_cat)) %>% 
    mutate(measurement_cat = gsub(", other", "", measurement_cat))

# Keeping track of the excluded studies ----------------------------------------
area_ts <- 
  area_ts %>% 
    # Exclude the studies that don't have a lat/lon (Only two studies I think)
    filter(!is.na(lat))

write_csv(area_ts, here::here('data_outputs/clean_area_ts_for_analysis.csv'))