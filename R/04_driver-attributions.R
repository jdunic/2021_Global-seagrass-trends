library(fuzzyjoin)

att_df <- 
  area_ts %>% 
    gather(key = attribution, value = yes_no, none, inferential, visual, descriptive) %>% 
    mutate(attribution = factor(attribution, levels = c('none', 'descriptive', 'visual', 'inferential'))) %>% 
    filter(yes_no == TRUE) %>% 
    distinct(studyid, attribution, yes_no, .keep_all = TRUE) %>% 
    select(studyid, site, author_ts, attribution, tidied_drivers, `driver effect`:apriori_cat) %>% 
    # Deal with one missing attribution category
    mutate(attribution = replace(attribution, (studyid == "114" & site == "Obidos Lagoon"), "descriptive")) 

# ------------------------------------------------------------------------------
# Main text Figure 4. 
# ------------------------------------------------------------------------------
dev.off()
dev.new(width = 8 / cm(1), height = 7 / cm(1))

att_method_plot <- 
area_ts %>% 
  select(studyid, site, att_method, none:descriptive) %>% 
  gather(key = attribution, value = yes_no, none, inferential, visual, descriptive) %>% 
  mutate(attribution = factor(str_to_sentence(attribution), levels = str_to_sentence(c('none', 'descriptive', 'visual', 'inferential')))) %>% 
  filter(yes_no == TRUE) %>% 
  distinct(studyid, attribution, yes_no) %>% 
ggplot(., aes(x = attribution)) +
  geom_histogram(stat = "count", aes(y = ..prop.., group = 1)) + 
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) + 
  mytheme(base_size = 12) +
  ylab('Proportion of studies') + 
  xlab('Attribution method') +
  theme(axis.title.x = element_text(vjust = -1), 
        plot.margin = margin(t = 3, r = 1, b = 7, l = 6))
att_method_plot
ggsave(here::here('figures/attribution_methods.png'), width = 8 / cm(1), height = 7 / cm(1))
ggsave(plot = att_method_plot, filename = here::here('figures/main_text/Figure-4_attribution_methods.eps'), 
       device = cairo_ps,
       width = 8 / cm(1), height = 7 / cm(1))



### Drivers

driver_cats <- read_csv(here::here("data/driver_cat_desc_lookup.csv"))

driver_cat_lookup <- 
  driver_cats %>% 
  mutate(detection_strs = gsub("; ", "|", detection_strs))

driver_info <- 
    att_df %>% 
    separate_rows(tidied_drivers, sep = ";") %>% 
    mutate(tidied_drivers = trimws(tidied_drivers)) %>% 
    regex_left_join(., driver_cat_lookup, by = c("tidied_drivers" = "detection_strs")) %>% 
    distinct(studyid, site, driver_cat, .keep_all = TRUE) %>%
    mutate(driver_cat = replace_na(driver_cat, "none"))

driver_att_counts <-   
  driver_info %>% 
    select(studyid, site, attribution, tidied_drivers, driver_cat, `primary driver`) %>% 
    filter(driver_cat != "none", attribution != "none") %>% 
    distinct(studyid, driver_cat, .keep_all = TRUE) %>%  
    count(driver_cat, attribution, name = "driver_att_count")

driver_study_counts <- 
  driver_att_counts %>% 
    group_by(driver_cat) %>% 
    summarise(driver_study_count = sum(driver_att_count))

driver_cat_counts <- 
  left_join(driver_att_counts, driver_study_counts) %>% 
  mutate(driver_cat = str_to_sentence(driver_cat), 
         attribution = str_to_sentence(attribution)) %>%
  mutate(attribution = factor(attribution, levels = c("Descriptive", "Visual", "Inferential"))) %>%
  mutate(driver_cat = str_replace(driver_cat, "Management/restoration", "Management/Restoration")) %>% 
  mutate(driver_att_label = paste0(driver_cat, " (", driver_study_count, ")")) %>%
  mutate(driver_cat = fct_reorder(driver_cat, driver_study_count)) 

# ------------------------------------------------------------------------------
# Main text Figure 5. 
# ------------------------------------------------------------------------------
driver_freq_plot <-
  driver_cat_counts %>% 
  distinct(driver_cat, driver_study_count) %>% 
  ggplot(data = ., aes(y = driver_cat, x = driver_study_count)) + 
    geom_histogram(stat = "identity") + 
    mytheme(14) + 
    theme(plot.margin = margin(t = 2, r = 0, b = 10, l = 5), 
          axis.title.x = element_text(vjust = 0)) + 
    labs(x = "Number of studies", y = "")
dev.off()
dev.new(width = 16.9 / cm(1), height = 10 / cm(1))
driver_freq_plot
ggsave(here::here('figures/primary-driver-frequency-plot.png'), width = 16.9 / cm(1), height = 10 / cm(1))
ggsave(plot = driver_freq_plot, filename = here::here('figures/main_text/Figure-5_primary-driver-frequency-plot.eps'), 
       device = cairo_ps,
       width = 16.9 / cm(1), height = 10 / cm(1))

# Attribution method by driver
# ------------------------------------------------------------------------------
dev.off()
dev.new(width = 16.9 / cm(1), height = 14 / cm(1))
driver_att_freq_plot <-
  driver_cat_counts %>% 
  mutate(driver_att_label = fct_reorder(driver_att_label, -driver_study_count)) %>%
  ggplot(data = ., aes(x = attribution, y = driver_att_count)) + 
    geom_histogram(stat = "identity") + 
    mytheme(10) + 
    theme(plot.margin = margin(t = 2, r = 0, b = 10, l = 5), 
          axis.title.x = element_text(vjust = -2)) +
    labs(x = "Attribution method", y = "Number of studies") + 
    facet_wrap(. ~ driver_att_label, ncol = 4)
driver_att_freq_plot
ggsave(here::here('figures/SOM/attribution-method_primary-driver-plot.png'), plot = driver_att_freq_plot, width = 16.9 / cm(1), height = 14 / cm(1))


# ------------------------------------------------------------------------------
# Look at management sites separately
# ------------------------------------------------------------------------------
management_sites <- 
  decadal_df %>% 
  filter(str_detect(tidied_drivers, "manage|restor")) %>% 
  distinct(study_site, site_groupings) %>% 
  filter(!(study_site %in% c("41_North control", "41_South control")))

events <- 
  read_csv(here::here("data/management_events.csv")) %>% 
  unite(study_site, studyid, site, sep = "_", remove = FALSE)
write_csv(events, here::here('tables/events_and_dates.csv'))

first_event <- 
  events %>% 
    group_by(study_site) %>% 
    slice(which.min(event_year)) %>% 
    ungroup()

# Look at raw event data
event_ts <- 
  area_ts %>% 
    filter(study_site %in% management_sites$study_site) %>%
    # Remove redundant Tampa Bay data: 
    filter(study_site != "87_Tampa Bay") %>% 
    left_join(., select(events, studyid, site, name)) %>% 
    tidyr::separate(col = study_site, into = c('studyid', 'site'),  sep = "_", remove =FALSE) %>% 
    mutate(site_groupings = ifelse(study_site == "W8_Hillsborough Bay - test plantings", "Hillsborough - test plantings", site_groupings)) %>%
    mutate(name = fct_reorder(name, -duration))

# Temporary mediocre solution
int_breaks <- function(x, n = 4) {
  l <- pretty(x, n)
  l[abs(l %% 1) < .Machine$double.eps ^ 0.5] 
}

dev.off()
dev.new(width = 18 / cm(1), height = 16.9 / cm(1))

event_ts_plot <- 
  ggplot(data = event_ts) + 
    aes(x = year, y  = area, group = study_site) + 
    geom_smooth(se = FALSE, colour = "grey80", size = 0.5) + 
    geom_point(size = 2) + 
    geom_vline(data = events %>% mutate(name = factor(name, levels = levels(event_ts$name))), 
      aes(xintercept = event_year), colour = "grey50") + 
    mytheme(base_size = 12) + 
    theme(axis.ticks = element_line(), 
          axis.title.y = element_text(vjust = -2), 
          plot.margin = unit(c(10, 5, 10, 2), "mm")) + 
    scale_y_log10() +
    scale_x_continuous(breaks = int_breaks) + 
    guides(colour = "none") + 
    labs(x = "Year", y = "Area (ha)") + 
    facet_wrap(~ name, ncol = 3, scale = "free")
event_ts_plot

ggsave(event_ts_plot, file = here("figures/SOM/event_ts_plot.png"), width = 15.5, height = 10)

# ------------------------------------------------------------------------------

