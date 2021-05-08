# ------------------------------------------------------------------------------
# Plot bioregional trends
# ------------------------------------------------------------------------------

# See: https://stat.ethz.ch/R-manual/R-patched/library/mgcv/html/mgcv-FAQ.html
# Notes: Confidence interval shrinks to nothing when straight line fit from GAM 
# crosses 0 (in our case predicted area = 1 ha: log(1)  = 0) because of 
# sum-to-zero identifiability constraint. Basically a straight line has only 
# one degree of freedom and so there is 'no choice' about where it passes through
# zero. Looking at figure 9 in Marra & Wood 2012, it sems like it is only the 
# confidence around the zero point that is not very informative/precise. At the 
# extremes, the use of bayeisan intervals seem to match the frequentist intervals
# (which I think are what is found using the mgcv package)

#source(here::here('R/03_decadal_data_prep.R'))

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

floor_decade <- function(value) { return(value - value %% 10) }

all_bioregions <-readRDS(file = here::here("data_outputs/GAM_bioregion_fits.RDS"))

initial_values <- 
  all_bioregions %>% 
    group_by(bioregion) %>% 
    arrange(x) %>%
    slice(1) %>% 
    rename(fit1 = fit, lwr1 = lwr, upr1 = upr) %>% 
    select(bioregion, fit1, lwr1, upr1)

all_bioregions_scaled <- left_join(all_bioregions, initial_values)

all_mins <- 
  all_bioregions_scaled %>%
  group_by(bioregion) %>% 
  slice(which.min(lwr)) %>% 
  mutate(min_lwr_exp = exp(lwr), 
         min_lwr_scaled = exp(lwr) / exp(fit1)) %>% 
  ungroup() %>% 
  left_join(., 
    all_bioregions_scaled %>%
      group_by(bioregion) %>% 
      slice(which.min(x)) %>% 
      mutate(min_year = floor_decade(x)) %>% 
      ungroup() %>% 
      select(bioregion, min_year)
  )

sites_per_decade <- 
  decadal_df %>% 
    distinct(study_site, decade, bioregion) %>% 
    group_by(bioregion, decade) %>% 
    summarise(n_sites = n()) %>% 
    ungroup() %>% 
    #unite(n_sites_label, decade, n_sites, sep = "\n", remove = FALSE) %>% 
    mutate(n_sites_label = paste0('(', n_sites, ')')) %>%
    left_join(., select(all_mins, bioregion, min_lwr_exp, min_lwr_scaled, min_year)) %>% 
    left_join(., bioregion_key) %>% 
    mutate(bioregion = factor(bioregion, levels = bioregion_levels)) %>% 
    filter(decade >= min_year)
    #mutate(n_sites_label = paste0(decade, '\n', '(', n_sites, ')'))


# Rescaling to initial values
# ------------------------------------------------------------------------------
n_sites_pos <- tibble(bioregion = bioregion_levels, y = c(0.001, 0.4, 0.05, 0.2, 0.15, 0.5, 0.5))
test <- 
  left_join(sites_per_decade, n_sites_pos) %>% 
  mutate(bioregion = factor(bioregion, levels = bioregion_levels))

dev.off()
dev.new(width = 20 / cm(1), height = 20 / cm(1))
bioregion_gams_plot <- 
  all_bioregions_scaled %>% 
    mutate(bioregion = factor(bioregion, levels = bioregion_levels)) %>% 
  ggplot(data = ., aes(x = x, y = exp(fit) / exp(fit1), fill = bioregion_colours)) + 
    #ggplot(data = ., aes(x = x, y = fit / fit1, fill = bioregion_colours)) + 
    geom_line() + 
    geom_ribbon(aes(ymin = exp(lwr) / exp(fit1), ymax = exp(upr) / exp(fit1)), alpha = 0.6) + 
    #geom_ribbon(aes(ymin = lwr / fit1, ymax = upr / fit1), alpha = 0.6) + 
    mytheme(base_size = 12) + 
    scale_fill_identity() +
    facet_wrap(~ bioregion, ncol = 2) + 
    scale_y_log10(breaks = c(0.001, 0.01, 0.1, 1, 5), labels = c(0.001, 0.01, 0.1, 1, 5)) +
    scale_x_continuous(name = "Year", breaks = scales::pretty_breaks(n = 6, min.n = 3)) + 
    theme(axis.ticks.x = element_line(colour = "grey50")) + 
    annotation_logticks(sides = 'l', colour = 'grey50', short = unit(0.1, 'cm'), long = unit(0.2, 'cm')) + 
    labs(y = "Mean meadow area relative to initial area", xlab = "Year")
  bioregion_gams_plot + 
    geom_text(data = test, aes(x = decade, y = 0.001, label = n_sites), size = 2.75, colour = "grey50")
ggsave(here::here("figures/bioregion_GAM_scaled_to_initial_fixed_axes_area.png"), width = 20 / cm(1), height = 20 / cm(1))

# ------------------------------------------------------------------------------
# Main text Figure 2. 
# ------------------------------------------------------------------------------
dev.off()
dev.new(width = 16.9 / cm(1), height = 20 / cm(1))
n_sites_pos <- tibble(bioregion = bioregion_levels, y = c(0.001, 0.41, 0.0425, 0.195, 0.153, 0.51, 0.51))
test <- 
  left_join(sites_per_decade, n_sites_pos) %>% 
  mutate(bioregion = factor(bioregion, levels = bioregion_levels))
bioregion_gam_plot2 <- 
  bioregion_gams_plot + 
  geom_text(data = test, aes(x = decade, y = y, label = n_sites), size = 3, colour = "grey45") +
  theme(strip.text = element_text(size = 10), 
        axis.title.x = element_text(vjust = -0.5),
        axis.title.y = element_text(vjust = -1)) + 
  facet_wrap_custom(~ bioregion, ncol = 2, scales = "free", scale_overrides = list(
    scale_override(1, scale_y_log10(breaks = c(0.001, 0.01, 0.1, 1), 
                                    labels = c(0.001, 0.01, 0.1, 1), 
                                    lim = c(0.0009, 6))),
    scale_override(2, scale_y_log10(breaks = c(0.5, 1, 2), 
                                    labels = c(0.5, 1, 2), 
                                    lim = c(0.4, 2))),
    scale_override(3, scale_y_log10(breaks = c(0.1, 1, 5), 
                                    labels = c(0.1, 1, 5), 
                                    lim = c(0.04, 5))),
    # scale_override(3, scale_y_log10(breaks = c(0.05, 0.1, 1, 10), 
    #                                 labels = c(0.05, 0.1, 1, 10))),
    scale_override(4, scale_y_log10(breaks = c(0.2, 0.5, 1, 2), 
                                    labels = c(0.2, 0.5, 1, 2), 
                                    lim = c(0.19, 2.3))),
    scale_override(5, scale_y_log10(breaks = c(0.2, 0.5, 1, 2), 
                                    labels = c(0.2, 0.5, 1, 2), 
                                    lim = c(0.15, 4))),
    scale_override(6, scale_y_log10(breaks = c(0.5, 1.0, 2), 
                                    labels = c(0.5, 1.0, 2), 
                                    lim = c(0.5, 2))),
    scale_override(7, scale_y_log10(breaks = c(0.5, 1.0, 2), 
                                    labels = c(0.5, 1.0, 2), 
                                    lim = c(0.5, 2))), 
    scale_override(1, scale_x_continuous(
                        breaks = seq(1950, 2010, by = 20), 
                        labels = seq(1950, 2010, by = 20))), 
    scale_override(2, scale_x_continuous(
                        breaks = seq(1930, 2010, by = 20), 
                        labels = seq(1930, 2010, by = 20))),
    scale_override(3, scale_x_continuous(
                        breaks = seq(1900, 2010, by = 20), 
                        labels = seq(1900, 2010, by = 20))), 
    scale_override(4, scale_x_continuous(
                        lim = c(1880, 2020),
                        breaks = seq(1880, 2010, by = 30), 
                        labels = seq(1880, 2010, by = 30))), 
    scale_override(5, scale_x_continuous(
                        breaks = seq(1940, 2010, by = 20), 
                        labels = seq(1940, 2010, by = 20))),
    scale_override(6, scale_x_continuous(
                        breaks = seq(1930, 2010, by = 20), 
                        labels = seq(1930, 2010, by = 20))), 
    scale_override(7, scale_x_continuous(
                        breaks = seq(1930, 2010, by = 20), 
                        labels = seq(1930, 2010, by = 20)))
    )) +
  theme(plot.margin = unit(c(5, 1, 5, 1), "mm")) +
  xlab("Year")
bioregion_gam_plot2
ggsave(here::here("figures/bioregion_GAM_scaled_to_initial.png"), width = 16.9 / cm(1), height = 20 / cm(1))
ggsave(plot = bioregion_gam_plot2, filename = here::here('figures/main_text/Figure-2_bioregion_GAM_scaled_to_initial.eps'), 
       device = cairo_ps,
       width = 16.9 / cm(1), height = 20 / cm(1))
