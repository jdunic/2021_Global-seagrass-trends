#source(here::here('R/03_decadal_data_prep.R'))
# Study count
study_count <- area_ts %>% distinct(studyid) %>% nrow

# Site count
site_count <- area_ts %>% distinct(site) %>% nrow

# Study counts make the most sense for discussing apriori categories because
# apriori categories were most often made at the study level
apriori_study_count <- 
  area_ts %>% 
    distinct(studyid, apriori_cat) %>% 
    .['apriori_cat'] %>% 
    table()
apriori_study_count

# Apriori cats across bioregions
apriori_bioregion_count <- 
  area_ts %>% 
    distinct(bioregion, studyid, apriori_cat) %>% 
    .[c('bioregion', 'apriori_cat')] %>% 
    table()
apriori_bioregion_count
saveRDS(apriori_bioregion_count, here::here('data_outputs/apriori_bioregion_count.RDS'))


# Histogram of site/study durations
# dev.off()
# dev.new(width = 11, height = 6)
# area_ts %>% 
#   distinct(study_site, duration) %>% 
# ggplot(data = ., aes(x = duration)) + 
#   geom_histogram(fill = "grey40") + 
#   mytheme() + 
#   scale_x_continuous(breaks = scales::pretty_breaks(n = 8)) +
#   labs(x = "Duration (years)", y = "Number of sites")
# ggsave(here::here('figures/SOM/duration_site_count_histogram.png'))

dev.off()
dev.new(width = 11, height = 6)
area_ts %>% 
  distinct(studyid, duration) %>% 
ggplot(data = ., aes(x = duration)) + 
  geom_histogram(fill = "grey40") + 
  mytheme() + 
  scale_x_continuous(breaks = scales::pretty_breaks(n = 8)) +
  labs(x = "Duration (years)", y = "Number of studies")
ggsave(here::here('figures/SOM/duration_study_count_histogram.png'))


duration_df <- 
  area_ts %>% 
    mutate(duration_bin = case_when(duration > 40 ~ '> 40', (duration > 10 & duration <= 40) ~ '10 - 40', 
                   duration <= 10 ~ '<= 10')) %>% 
    mutate(duration_bin = factor(duration_bin, levels = c('<= 10', '10 - 40', '> 40')))

duration_df %>% 
  distinct(studyid, study_site, duration_bin) %>%
  janitor::tabyl(duration_bin) %>% 
  mutate(percent = round(percent*100)) %>% 
  rename('sites' = n, 'site_prop' = percent) %>% 
  left_join(., duration_df %>% distinct(studyid, duration_bin) %>% janitor::tabyl(duration_bin)) %>% 
  mutate(percent = round(percent*100)) %>% 
  rename('studies' = n, 'study_prop' = percent) %>% 
saveRDS(., 'data_outputs/duration_table.RDS')


# Relationship between rates of change and fractional size
decadal_df %>% 
  filter(fractional_change < 1) %>%
  distinct(study_site, decade, .keep_all = TRUE) %>% 
ggplot(data = .) + 
  aes(x = fractional_change, y = trend_ci_50) + 
  geom_point() + 
  mytheme() + 
  labs(y = "Annual rate of change (%)", 
       x = "Fraction of maximum observed area")
ggsave(here::here('figures/SOM/rates_of_change-fraction_of_max_scatterplot.png'))



dev.off()
dev.new(width = 16.9 / cm(1), height = 10 / cm(1))
x_labels <- c('1*m^2', 'phantom(0)', '100*m^2', 'phantom(0)', '1*ha^phantom(0)', 'phantom(0)', '100*ha^phantom(0)', 'phantom(0)', '100*km^2', 'phantom(0)', 'phantom(0)~~10000*km^2')
area_ts %>% 
  group_by(study_site) %>% 
  mutate(cv = sd(plus_10) / mean(plus_10)) %>% 
  distinct(study_site, .keep_all = TRUE) %>% 
ggplot(data = ., aes(x = max_area, y = cv)) + 
  geom_point(alpha = 0.5) + 
  scale_x_log10() + 
  mytheme(base_size = 14) + 
  labs(y = "Coefficient of variation", x = "Maximum observed area") + 
  scale_x_log10(breaks = lseq(from = 1e-4, to = 1e6, length.out = 11),
                labels = parse(text = x_labels), 
                limits = c(1e-4, 1.5e6))
ggsave('figures/SOM/cv_max_area.png', width = 16.9 / cm(1), height = 10 / cm(1))


# 
dev.off()
dev.new(width = 16.9 / cm(1), height = 10 / cm(1))

p1 <-
decadal_df_ne_1930 %>% 
filter(trend_ci_50 < 0) %>%
ggplot(data = .) +
  mytheme(base_size = 12) +
  aes(x = frac_of_max, y = abs(trend_ci_50), size = NULL) +
  geom_vline(xintercept = 1, colour = "grey70") + 
  geom_point(size = 1, alpha = 0.5) + 
  labs(y = "Absolute annual rate of change (%)", 
       x = "Fraction of maximum observed area") + 
  ylim(c(0, 85)) +
  xlim(c(0, 1.15)) +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_text(vjust = 4),
        plot.title = element_text(size = 12)) + 
  ggtitle("Rates of Loss")
p2 <-
decadal_df_ne_1930 %>% 
filter(trend_ci_50 > 0) %>%
filter(trend_ci_50 < 100) %>% 
ggplot(data = .) +
  mytheme(base_size = 12) +
  aes(x = frac_of_max, y = abs(trend_ci_50), size = NULL) +
    geom_vline(xintercept = 1, colour = "grey70") + 
  geom_point(size = 1, alpha = 0.5) + 
  theme(axis.title.y = element_blank(), 
        axis.title.x = element_blank(), 
        plot.title = element_text(size = 12)) + 
  labs(y = "Absolute annual rate of change (%)", 
       x = "Fraction of maximum observed area") + 
  ylim(c(0, 85)) +
  xlim(c(0, 1.15)) +
  ggtitle("Rates of Gain")
p1 + p2 +
plot_annotation(caption = "Fraction of maximum observed area", theme = theme(plot.title = element_blank(), plot.caption = element_text(size = 12, hjust = 0.5, vjust = 3)))
ggsave(here::here('figures/SOM/abs_loss_gain_funnel.png'))