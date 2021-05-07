# Cumulative change in area overall and per decade

# they calculated that. 
# Would it be fair to: 
# • Calculate all the losses and all the gains between time points. So if at a single site biomass went up from 1991 to 1992, we would count that as a gain, if it then went down from 1992 to 1993, we would count that as a loss
# • Plot cumulative losses and gains per decade (so just sum up all of the losses as one variable for the graph and sum up all of the gains as another variable). 
# • Plot cumulative losses/gains as a % of area surveyed in that year. 

# Is it ok to have things relative to maximum observed area? Or does this bias 
# things solely towards loss?

#
# Observed number of sites with loss and gain :over entire time series
# ------------------------------------------------------------------------------
#source(here::here('R/03_decadal_data_prep.R'))
library(flextable)

duplicated_sites <- 
  filter(area_ts, site_duplication == "yes")

bioregion_newlines <- str_replace(bioregion_levels, "North Atlantic", "North\nAtlantic")
bioregion_newlines <- str_replace(bioregion_newlines, "Temperate Southern", "Temperate\nSouthern")
bioregion_newlines <- str_replace(bioregion_newlines, "Temperate North Pacific", "Temperate\nNorth Pacific")

# To deal with some sites being replicated across studies, use only the studies
# with the most recent sampling date.
get_area_change_df <- function(df) {
  df %>% 
  distinct(study_site, .keep_all = TRUE) %>% 
  mutate(net_change = A2 - A1,
         net_change_sign = case_when(net_change > 0 ~ "gain", 
                                     net_change == 0 ~ "no change", 
                                     net_change < 0 ~ "loss")) %>% 
  mutate(decade = ifelse(as.numeric(decade) < 1940, "pre-1940s", decade)) %>% 
  mutate(decade = factor(decade, levels = decade_levels)) %>% 
  left_join(., bioregion_key) %>% 
  mutate(bioregion_levels = factor(bioregion, levels = bioregion_levels)) %>%
  mutate(bioregion_striptext = str_replace(bioregion_levels, "North Atlantic", "North\nAtlantic")) %>% 
  mutate(bioregion_striptext = factor(bioregion_striptext, levels = bioregion_newlines))
}

get_global_change <- function(area_change_df) {
    area_change_df %>% 
    group_by(net_change_sign) %>% 
    summarise(net_change = sum(net_change)) %>% 
    ungroup() %>% 
    spread(key = net_change_sign, value = net_change) %>% 
    select(-`no change`) %>% 
    mutate(bioregion = 'Global', 
           net_change = gain + loss, 
           total_area = sum(area_change_df$max_area)) %>% 
    mutate(percent_change = 100 * (net_change / total_area)) %>% 
    mutate(time_span = paste0(min(area_change_df$first_year), " - ", max(area_change_df$last_year)))
}

get_bioregion_change <- function(area_change_df) {
  area_change_df %>% 
      group_by(bioregion, net_change_sign) %>% 
      summarise(net_change = sum(net_change)) %>% 
      ungroup() %>% 
      spread(key = net_change_sign, value = net_change)
}

get_bioregion_total_area <- function(area_change_df)
  area_change_df %>% 
    distinct(study_site, .keep_all = TRUE) %>% 
    group_by(bioregion) %>% 
    summarise(total_area = sum(max_area), 
              min_year = min(first_year), 
              max_year = max(last_year)) %>% 
    mutate(time_span = paste(min_year, max_year, sep = " - ")) %>%
    ungroup()

get_area_table <- function(area_change_df) {
  bioregion_net_change <- get_bioregion_change(area_change_df)
  bioregion_total_area <- get_bioregion_total_area(area_change_df)
  total_change <- get_global_change(area_change_df)
  output <- 
  left_join(bioregion_net_change, bioregion_total_area) %>% 
    mutate(net_change = gain + loss) %>% 
    mutate(percent_change = 100 * (net_change / total_area)) %>% 
    arrange(percent_change) %>% 
    select(bioregion, percent_change, net_change, gain, loss, total_area, time_span) %>% 
    bind_rows(., total_change) %>%
    mutate(percent_change = round(percent_change, 1)) %>%
    mutate_at(c("net_change", "gain", "loss", "total_area"), round, 0) %>% 
  return(output)
}


# Most recent
# -----------------
most_recent_site_groupings <- 
  duplicated_sites %>%
    group_by(site_groupings, studyid) %>% 
    summarise(max_max_area = max(max_plus_0.1), most_recent = max(last_year)) %>%
    # Get the most recent last sampled date, or the largest max area if the most
    # recent dates are the same time
    ungroup() %>%
    group_by(site_groupings) %>%
    arrange(site_groupings, -most_recent, -max_max_area) %>% 
    slice(1) %>%
    ungroup() %>% 
    mutate(recent_max_area = TRUE)

recent_area_change <- 
left_join(area_ts, most_recent_site_groupings, by = c('site_groupings', 'studyid')) %>% 
  #left_join(., traj_lms %>% select(study_site, trajectory)) %>% 
  filter(site_duplication == "no" | recent_max_area == TRUE) %>% 
  get_area_change_df() %>% 
  get_area_table()

view(recent_area_change)


# Strategy 2 - look at only the largest sites
max_area_site_groupings <- 
  duplicated_sites %>%
    group_by(site_groupings, studyid) %>% 
    summarise(max_max_area = max(max_plus_0.1), most_recent = max(last_year)) %>%
    # Get the most recent last sampled date, or the largest max area if the most
    # recent dates are the same time
    ungroup() %>%
    group_by(site_groupings) %>%
    arrange(site_groupings, -max_max_area) %>% 
    slice(1) %>%
    ungroup() %>% 
    mutate(recent_max_area = TRUE)

max_area_area_change <- 
left_join(area_ts, max_area_site_groupings, by = c('site_groupings', 'studyid')) %>% 
  #left_join(., traj_lms %>% select(study_site, trajectory)) %>% 
  filter(site_duplication == "no" | recent_max_area == TRUE) %>% 
  get_area_change_df() %>% 
  get_area_table()

view(max_area_area_change)

saveRDS(max_area_area_change, "tables/net_change_max_area.RDS")

# Strategy 3 - look at longest duration sites
duration_site_groupings <- 
  duplicated_sites %>%
    group_by(site_groupings, studyid) %>% 
    summarise(max_max_area = max(max_plus_0.1), most_recent = max(last_year), max_duration = max(duration)) %>%
    # Get the most recent last sampled date, or the largest max area if the most
    # recent dates are the same time
    ungroup() %>%
    group_by(site_groupings) %>%
    arrange(site_groupings, -max_duration) %>% 
    slice(1) %>%
    ungroup() %>% 
    mutate(recent_max_area = TRUE)

max_duration_area_change <- 
left_join(area_ts, duration_site_groupings, by = c('site_groupings', 'studyid')) %>% 
  #left_join(., traj_lms %>% select(study_site, trajectory)) %>% 
  filter(site_duplication == "no" | recent_max_area == TRUE) %>% 
  get_area_change_df() %>% 
  get_area_table()


# -------------------------------------------------------------------------------
# Plot net change distributions
# -------------------------------------------------------------------------------
decade_levels <- c('pre-1940s', as.character(seq(1940, 2010, by = 10)))

area_change <- 
  left_join(area_ts, max_area_site_groupings, by = c('site_groupings', 'studyid')) %>% 
  #left_join(., traj_lms %>% select(study_site, trajectory)) %>% 
  filter(site_duplication == "no" | recent_max_area == TRUE) %>% 
  distinct(study_site, .keep_all = TRUE) %>% 
  mutate(net_change = A2 - A1,
         net_change_sign = case_when(net_change > 0 ~ "gain", 
                                     net_change == 0 ~ "no change", 
                                     net_change < 0 ~ "loss")) %>% 
  mutate(decade = ifelse(as.numeric(decade) < 1940, "pre-1940s", decade)) %>% 
  mutate(decade = factor(decade, levels = decade_levels)) %>% 
  left_join(., bioregion_key) %>% 
  mutate(bioregion_levels = factor(bioregion, levels = bioregion_levels)) %>%
  mutate(bioregion_striptext = str_replace(bioregion_levels, "North Atlantic", "North\nAtlantic")) %>% 
  mutate(bioregion_striptext = str_replace(bioregion_striptext, "Temperate Southern", "Temperate\nSouthern")) %>% 
  mutate(bioregion_striptext = str_replace(bioregion_striptext, "Temperate North Pacific", "Temperate\nNorth Pacific")) %>% 
  mutate(bioregion_striptext = factor(bioregion_striptext, levels = bioregion_newlines))

write_csv(area_change, here::here('data_outputs/area_change.csv'))


# https://stackoverflow.com/questions/23901907/create-a-log-sequence-across-multiple-orders-of-magnitude
lseq <- function(from=1, to=100000, length.out=10) {
  # logarithmic spaced sequence
  # blatantly stolen from library("emdbook"), because need only this
  exp(seq(log(from), log(to), length.out = length.out))
}

# Prepare base plot to get size bins for dot plot
# ------------------------------------------------------------------------------
net_change_plot <- 
  ggplot() + 
    geom_histogram(data = area_change %>% filter(net_change_sign == "gain"), aes(y = net_change, x = ..count.., fill = net_change_sign), alpha = 0.85) + 
    geom_histogram(data = area_change %>% filter(net_change_sign == "loss"), aes(y = abs(net_change), x = -..count.., fill = net_change_sign), alpha = 0.85) +
    #geom_text(data = tibble(x = 1e-4, y = c(30, -30), label = c("Gain", "Loss")), 
             #geom = "text", aes(x = x, y = y, label = label), size = 10) +
    guides(fill = "none") +
    scale_fill_manual(values = c("#2c7bb6", "#d7191c")) + 
    geom_hline(yintercept = 0, colour = "grey30") + 
    mytheme(base_size = 20)
dev.off()
dev.new(width = 10, height = 8)
y_labels <- c('1~m^2', '10~m^2', '100~m^2', '1000~m^2', '1~ha^phantom(0)', '10~ha^phantom(0)', '100~ha^phantom(0)', '10~km^2', '100~km^2', '1000~km^2', 'phantom(0)~~10000~km^2')
net_change_plot2 <- 
net_change_plot + 
  scale_y_log10("Net change in meadow extent", 
                breaks = lseq(from = 1e-4, to = 1e6, length.out = 11),
                labels = parse(text = y_labels)) +
  scale_x_continuous("Number of meadows", 
        breaks = seq(from = -30, to = 30, by = 10), 
        labels = abs(seq(from = -30, to = 30, by = 10)), 
        limits = c(-32, 32)) + 
  annotation_logticks(side = "l", colour = "grey70") + 
  theme(plot.margin = unit(c(10, 5, 10, 5), "mm"))

# Global net change plot -- Use dot plot with size = net change area
# ------------------------------------------------------------------------------
width_val <- ggplot_build(net_change_plot2)$data[[1]]$ymax[1] - ggplot_build(net_change_plot2)$data[[1]]$ymin[1]

lower_str1 <- paste0("\\(", signif(ggplot_build(net_change_plot2)$data[[1]]$ymin[1], 3))
lower_str2 <- paste0("\\[", signif(ggplot_build(net_change_plot2)$data[[1]]$ymin[1], 3))

dot_binwidths <- 
  ggplot_build(net_change_plot2)$data[[1]] %>% 
    select(y, ymin, ymax) %>% 
    mutate(dot_bins = paste0("(", signif(ymin, 3), ",", signif(ymax, 3), "]")) %>% 
    mutate(dot_bins = str_replace(dot_bins, lower_str1, lower_str2)) %>% 
    #mutate(bw = ymax - ymin) %>% 
    select(y, dot_bins)

dot_df <- 
area_change %>% 
  filter(net_change_sign != "no change") %>%
  mutate(dot_bins = cut_width(log10(abs(net_change)), width = width_val, dig.lab = 3)) %>%
  select(bioregion_levels, bioregion_striptext, net_change, net_change_sign, study_site, dot_bins) %>% 
  group_by(dot_bins, net_change_sign) %>% 
  arrange(net_change) %>% 
  mutate(row_number = row_number()) %>% 
  ungroup() %>% 
  left_join(., dot_binwidths) %>%
  mutate(net_change_sign = factor(net_change_sign, levels = c("loss", "no change", 'gain')))

dot_bioregion_df <- 
area_change %>% 
  filter(net_change_sign != "no change") %>%
  mutate(dot_bins = cut_width(log10(abs(net_change)), width = width_val, dig.lab = 3)) %>%
  select(bioregion_levels, bioregion_striptext, net_change, net_change_sign, study_site, dot_bins) %>%
  group_by(dot_bins, net_change_sign, bioregion_levels) %>% 
  arrange(net_change) %>% 
  mutate(row_number = row_number()) %>% 
  ungroup() %>%
  left_join(., dot_binwidths) %>% 
  mutate(net_change_sign = factor(net_change_sign, levels = c("loss", "no change", 'gain')))

dev.off()
dev.new(width = 11, height = 7)
ggplot() + 
  geom_point(
      data = filter(dot_df, net_change_sign == "gain"), 
      aes(y = 10^(y), x = row_number, size = abs(net_change), colour = net_change_sign), 
      alpha = 0.5) + 
  geom_point(
      data = filter(dot_df, net_change_sign == "loss"), 
      aes(y = 10^(y), x = -1 * row_number, size = (abs(net_change)), colour = net_change_sign), 
      alpha = 0.5) + 
  scale_size(range = c(1, 20)) +
  scale_size(trans = "sqrt", range = c(1, 10)) + 
  scale_y_log10() + 
  guides(size = "none", colour = "none") + 
  scale_colour_manual(values = c("#2c7bb6", "#d7191c")) + 
  geom_vline(xintercept = 0, colour = "grey30") + 
  mytheme(base_size = 20) + 
  scale_y_log10("Net change in meadow extent", 
                breaks = lseq(from = 1e-4, to = 1e6, length.out = 11),
                labels = parse(text = y_labels)) +
  scale_x_continuous("Number of meadows", breaks = seq(from = -30, to = 30, by = 10), limits = c(-31, 31), labels = abs(seq(from = -30, to = 30, by = 10))) + 
  annotation_logticks(side = "l", colour = "grey70") + 
  theme(axis.text.y = element_text(vjust = 0.25))
ggsave(here("figures/SOM/net_change_all_sites_raw_data_point_size_sqrt.png"), width = 11, height = 7)

# Bioregion facets
# ------------------------------------------------------------------------------
dev.off()
dev.new(width = 16.9 / cm(1), height = 17 / cm(1))
ggplot() + 
  geom_point(
      data = filter(dot_bioregion_df, net_change_sign == "gain"), 
      aes(y = 10^(y), x = row_number, size = abs(net_change), colour = net_change_sign), 
      alpha = 0.5) + 
  geom_point(
      data = filter(dot_bioregion_df, net_change_sign == "loss"), 
      aes(y = 10^(y), x = -1 * row_number, size = (abs(net_change)), colour = net_change_sign), 
      alpha = 0.5) + 
  scale_size(trans = "sqrt", range = c(1, 6)) + 
  scale_y_log10() + 
  guides(size = "none", colour = "none") + 
  scale_colour_manual(values = c("#2c7bb6", "#d7191c")) + 
  geom_vline(xintercept = 0, colour = "grey30") + 
  mytheme(base_size = 12) + 
  scale_x_continuous("Number of meadows", 
                     breaks = c(-10, -5, 0, 5, 10), 
                     labels = c(10, 5, 0, 5, 10),
                     limits = c(-13, 13)) +
  scale_y_log10("Net change in meadow extent", 
                breaks = lseq(from = 1e-4, to = 1e6, length.out = 11),
                labels = parse(text = y_labels), 
                limits = c(1e-4,1.5e6)) +
  annotation_logticks(side = "l", colour = "grey70") + 
  theme(plot.margin = unit(c(3, 1, 5, 1), "mm"),
        axis.text.y = element_text(vjust = 0.25), 
        axis.title.x = element_text(vjust = -1),
        axis.title.y = element_text(vjust = -1)) + 
  facet_wrap(~ bioregion_striptext, ncol = 4)
ggsave(here("figures/SOM/net_change_by_bioregion_raw_data_point_size_sqrt.png"), width = 16.9 / cm(1), height = 17 / cm(1))


# ------------------------------------------------------------------------------
# Survey effort
# ------------------------------------------------------------------------------
y_labels <- c('phantom(0)', '10~m^2', 'phantom(0)', 'phantom(0)', 'phantom(0)~1~ha^phantom(0)', 'phantom(0)', '1~km^2', '10~km^2', '100~km^2', '1000~km^2', '10000~km^2')
dev.off()
dev.new(height = 18 / cm(1), width = 16.9 / cm(1))
sampled_max_area_plot <- 
select(area_ts, bioregion, site_groupings, studyid, study_site, decade, max_area) %>% 
  distinct(study_site, decade, .keep_all = TRUE) %>% 
  group_by(bioregion, decade) %>% 
  summarise(total_area_surveyed = sum(max_area)) %>% 
  mutate(total_area_surveyed_km2 = convertr::convert(total_area_surveyed, "ha", "km2")) %>% 
ggplot(data = ., aes(x = decade, y = total_area_surveyed_km2)) + 
  geom_line() +
  scale_colour_identity() +
  scale_x_continuous(name = "Decade", breaks = pretty_breaks(n = 6)) +
  scale_y_log10(expression(paste("Total area surveyed", " (km"^2, ")", sep = "")), 
                breaks = lseq(from = 1e-4, to = 1e6, length.out = 11),
                labels = trans_format("log10", math_format(10^.x))) +
  mytheme(base_size = 12) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.title.y = element_text(vjust = 2), 
        axis.ticks = element_line()) +
  facet_wrap(~ bioregion)
sampled_max_area_plot
ggsave(here("figures/SOM/sampled_survey_effort_over_time_by_bioregion.png"), height = 18 / cm(1), width = 16.9 / cm(1))
