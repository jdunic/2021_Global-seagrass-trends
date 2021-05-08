# Status plots 
# ----------------- 
#source(here::here('R/03_decadal_data_prep.R'))

# ------------------------------------------------------------------------------
# Global - all points density plot with blue +/- 2% ribbon (Figure S7)
# ------------------------------------------------------------------------------
dev.off()
dev.new(height = 16.9 / cm(1), width = 16.9 / cm(1))

ggplot(data = decadal_df_ne_1930) + 
  aes(x = frac_of_max, y = trend_ci_50) +
  geom_rect(data = NULL, aes(xmin = -Inf, xmax = Inf, ymin = -2, ymax = 2), fill = "lightblue", alpha = 0.3) +
  stat_density_2d(geom = "polygon", contour = TRUE, size = 0.5,
                  aes(fill = after_stat(level)), colour = "black",
                  bins = 5) +
  scale_fill_distiller(palette = "Oranges", direction = 1) +
  geom_point(aes(size = max_area_km2), colour = "black", alpha = 0.5) +
  geom_hline(yintercept = 0, colour = "grey60") + 
  geom_vline(xintercept = 1, linetype = 1, colour = "grey60") + 
  labs(x = "Fraction of maximum observed area", 
       y = "Annual rate of change (%)") +
  mytheme(base_size = 14) + 
  theme(legend.position = "bottom", 
        #axis.text = element_text(size = 12), 
        panel.grid.major.y = element_line(colour = "grey95"),
        panel.background = element_rect(fill = "grey98", colour = NA)) + 
  scale_x_continuous(breaks = seq(0, 1.3, 0.1), 
                     labels = c(0, '', '', '', '', 0.5, '', '', '', '', 1, '', '', ''), 
                     lim = c(0, 1.2)) + 
  scale_y_continuous(trans = pseudo_log_trans(sigma = 0.5),
                     #lim = c(-400, 400),
                     #breaks = c(-200, -100, -50, -20, -10, -5, -2, 0, 2, 5, 10, 20, 50, 100, 200)) + 
                     breaks = c(-200, -50, -10, -2, 0, 2, 10, 50, 200)) + 
  # scale_size_area(trans = "log10", breaks = 10 * c(0.1, 1, 10, 100, 1000),
  #                 labels = c('\u22640.1', '1', '10', '100', '1000'), 
  #                 max_size = 4) +
  # scale_size_area(trans = "log10", breaks = c(1, 10, 100, 1000),
  #                 labels = c(expression(paste("\u2264", " 0.1", " km"^2)), 
  #                            expression(paste('1-10', " km"^2)), 
  #                            expression(paste('10-100', ' km'^2)), 
  #                            expression(paste()) '10-100 ha', '\u22651000 ha'), 
  #                 max_size = 4) +
  # Sqrt transformation of size did a better job of 
  scale_size(trans = "sqrt", 
             range = c(0.5, 6), 
             breaks = c(0.1, 1, 10, 100, 1000),
             limits = c(1e-6, 6000), 
             labels = c('0.1', '1', '10', '100', '1000')) + 
  theme(legend.key = element_rect(fill = "transparent", colour = "transparent")) +
      theme(axis.title.x = element_text(vjust = -1, margin = margin(t = 0, r = 0, b = 0, l = 0)), 
            legend.title = element_text(size = 12),
            legend.box.margin = margin(t = 10, r = 15, b = -5, l = 0),
            legend.margin = margin(t = -5, r = 0, b = -5, l = 0),
            legend.text = element_text(hjust = 0.5, margin = margin(r = -4, l = -5))) +
  guides(size = guide_legend(override.aes = list(alpha = 1)), 
         fill = "none") +
  labs(size = expression(paste("Maximum observed area", " (km"^2, ")", sep = "")))

# Figure S7
ggsave(here("figures/SOM/global_all-time_density_2percent.png"))

# ------------------------------------------------------------------------------
# Make unmodified status plots for each bioregion
# ------------------------------------------------------------------------------
bioregion_point_colour = NULL
#bioregion_point_colour = quo(studyid)
#bioregion_point_colour = quo(site_groupings)
#bioregion_point_colour = quo(waycott_study)
#event_timeline = quo(waycott_study)
axis_text_size = 8
#
custom_layers <- function(ggplot_obj) {
  ggplot_obj +
  facet_grid(cols = vars(decade)) + 
  guides(colour = "none", size = "none") +
  theme(strip.text = element_blank(),
        axis.text.y = element_text(size = axis_text_size),
        axis.text.x = element_blank(), 
        axis.title.y = element_text(size = 10), 
        axis.title.x = element_blank(),
        plot.title = element_text(size = 8, margin = margin(t = 3, b = 1)), 
        plot.margin = unit(c(0, 0, 1, 3), "mm")) +
  # scale_size_area(trans = "log10", 
  #                 breaks = 10 * c(0.1, 1, 10, 100, 1000),
  #                 limits = c(1, 60000), 
  #                 labels = c('\u22640.1', '1', '10', '100', '1000'), 
  #                 max_size = 3) +
  scale_size(trans = "sqrt", 
             range = c(0.5, 5), 
             breaks = c(0.1, 1, 10, 100, 1000),
             limits = c(1e-6, 6000), 
             labels = c('0.1', '1', '10', '100', '1000')) + 
  #scale_size_binned(breaks = c(0.1, 1, 10, 100, 1000), range = c(0.5, 3)) + 
  theme(axis.title.y = element_blank()) +
  labs(size = expression(paste("Maximum observed area", " (km"^2, ")", sep = "")))
  #ylim(c(-100, 200))  # You lose too much detail when putting everything on the
  # same y-axis scale; variability in the TNAW outweighs everything else
}

status_tnp <- 
  mk_status_plot(df = decadal_df_ne_1930 %>% filter(bioregion == "Temperate North Pacific"), 
                 point_colour = !!bioregion_point_colour, size_bins = (max_area_km2)) %>% 
  custom_layers() + 
  geom_point(data = tibble(waycott_study = FALSE, decade = c(1930, 1940), trend_ci_50 = c(0, 0), frac_of_max = c(1, 1)), alpha = 0) +
  ggtitle("Temperate North Pacific") + 
  theme(strip.text = element_text(size = 8))

status_tnaw <- 
  mk_status_plot(df = decadal_df_ne_1930 %>% filter(bioregion == "Temperate North Atlantic West"), 
                 point_colour = !!bioregion_point_colour, size_bins = max_area_km2) %>%
  custom_layers() +
  ggtitle("Temperate North Atlantic West")

status_tnae <- 
  mk_status_plot(df = decadal_df_ne_1930 %>% filter(bioregion == "Temperate North Atlantic East"), 
                 point_colour = !!bioregion_point_colour, size_bins = max_area_km2) %>% 
  custom_layers() +
  ggtitle("Temperate North Atlantic East")

status_m <- 
  mk_status_plot(df = decadal_df_ne_1930 %>% filter(bioregion == "Mediterranean"), 
                 point_colour = !!bioregion_point_colour, size_bins = max_area_km2) %>%
    custom_layers() +
    ggtitle("Mediterranean") + 
    ylab("Annual rate of change (% per year)") + 
    theme(axis.title.y = element_text(size = 10, angle = 90, vjust = 2))

status_ta <- 
  mk_status_plot(df = decadal_df_ne_1930 %>% filter(bioregion == "Tropical Atlantic"), 
                 point_colour = !!bioregion_point_colour, size_bins = max_area_km2) %>%
  custom_layers() +
  ggtitle("Tropical Atlantic")

status_tip <- 
  mk_status_plot(df = decadal_df_ne_1930 %>% filter(bioregion == "Tropical Indo-Pacific"), 
                 point_colour = !!bioregion_point_colour, size_bins = max_area_km2) %>% 
    custom_layers() +
    ggtitle("Tropical Indo-Pacific")

status_so <- 
  mk_status_plot(df = decadal_df_ne_1930 %>% filter(bioregion == "Temperate Southern Oceans"), 
                 point_colour = !!bioregion_point_colour, size_bins = max_area_km2) %>% 
    custom_layers() +
    guides(colour = guide_legend(title = "Dataset", nrow = 2), size = guide_legend(label.hjust = 0, override.aes = list(alpha = 1))) + 
    theme(axis.title.x = element_text(size = 10, margin = margin(t = 5, r = 0, b = -8, l = 0)), 
          axis.text.x = element_text(size = 8), 
          legend.title = element_text(size = 10, margin = margin(r = 10)),
          legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0),
          legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
          legend.key.width = unit(2, "mm"),
          legend.text = element_text(size = 10, margin = margin(r = 5, l = 0))) +
    ggtitle("Temperate Southern Oceans")
# ------------------------------------------------------------------------------
# Main text Figure 3. 
# ------------------------------------------------------------------------------
dev.off()
dev.new(width = 16.9 / cm(1), height = 20 / cm(1))
status_tnp + status_tnaw + status_tnae + status_m + status_ta + status_tip + status_so + 
  plot_layout(ncol = 1)
ggsave(here("figures/all_bioregions_pseudo_log_y_sigma_0.5.png"), width = 16.9 / cm(1), height = 20 / cm(1))
fig_3 <- status_tnp + status_tnaw + status_tnae + status_m + status_ta + status_tip + status_so + 
  plot_layout(ncol = 1)
ggsave(plot = fig_3, filename = here::here('figures/main_text/Figure-3_all_bioregions_pseudo_log_y_sigma_0.5.eps'), 
       device = cairo_ps,
       width = 16.9 / cm(1), height = 20 / cm(1))

# Status by meadow size bin
# ------------------------------------------------------------------------------
dev.off()
dev.new(width = 16.9 / cm(1), height = 20 / cm(1))
status_size_bins_ne <- 
    mk_status_plot(df = decadal_df_ne_1930, 
                   size_bins = max_area_km2, point_colour = bioregion) +
    facet_grid(rows = vars(size_bins_km2), cols = vars(decade), scale = "free_y") + 
    guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2), title.position = "top", nrow = 4), 
           size = guide_legend( override.aes = list(alpha = 1), title.position = "top")) +
    #ggtitle("All bioregions") + 
    #scale_colour_manual(values = bioregion_colours) + 
    # scale_size_area(trans = "log10", breaks = 10 * c(0.1, 1, 10, 100, 1000), max_size = 3,
    #                 labels = c('\u22640.1', '1', '10', '100', '1000')) +
    scale_size(trans = "sqrt", 
             range = c(0.5, 5), 
             breaks = c(0.1, 1, 10, 100, 1000),
             limits = c(1e-6, 6000), 
             labels = c('0.1', '1', '10', '100', '1000')) + 
    labs(y = "Annual rate of change (%)", colour = "Bioregion",
         size = expression(paste("Maximum observed area", " (km"^2, ")", sep = ""))) + 
    #mytheme(base_size = 12) +
    theme(axis.title.x = element_text(size = 10, margin = margin(t = 5, r = 0, b = 15, l = 0)), 
          axis.text.x = element_text(size = 8), 
          plot.margin = unit(c(3, 5, 3, 5), "mm"),
          legend.title = element_text(size = 10, margin = margin(r = 10)),
          legend.position = "bottom",
          legend.box.just = "left",
          legend.box = "vertical",
          legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0),
          legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
          legend.key.width = unit(2, "mm"),
          legend.text = element_text(size = 10, margin = margin(r = 12, l = 0))) +
    geom_hline(yintercept = c(-2, 2), alpha = 0)
  # scale_y_continuous(trans = pseudo_log_trans(sigma = 0.5),
  #                    lim = c(-400, 400),
  #                    breaks = c(-300, -100, -50, -20, -10, -5, -2, 0, 2, 5, 10, 20, 50, 100, 300))
status_size_bins_ne
ggsave(here("figures/SOM/size-bin_decade_bioregions_pseudo_log_y_sigma_0.5.png"))


# ------------------------------------------------------------------------------
# Compare Waycott studies, pre 2007, full dataset status plots
# ------------------------------------------------------------------------------
bioregion_point_colour = quo(waycott_study)
#event_timeline = quo(waycott_study)
axis_text_size = 8
#
custom_layers <- function(ggplot_obj) {
  ggplot_obj +
  facet_grid(cols = vars(decade)) + 
  guides(colour = "none", size = "none", shape = "none") +
  theme(strip.text = element_blank(),
        axis.text.y = element_text(size = axis_text_size),
        axis.text.x = element_blank(), 
        axis.title.y = element_text(size = 10), 
        axis.title.x = element_blank(),
        plot.title = element_text(size = 8, margin = margin(t = 3, b = 1)), 
        plot.margin = unit(c(0, 0, 1, 3), "mm")) +
  # scale_size_area(trans = "log10", breaks = c(1, 10, 100, 1000),
  #                 limits = c(0.000001, 6000), 
  #                 labels = c('1', '10', '100', '1000'), 
  #                 max_size = 3) +
  scale_size(trans = "sqrt", 
             range = c(0.75, 4), 
             breaks = c(0.1, 1, 10, 100, 1000),
             limits = c(1e-6, 6000), 
             labels = c('0.1', '1', '10', '100', '1000')) + 
  theme(axis.title.y = element_blank()) +
  scale_colour_manual(labels = c('Dunic et al.', 'Waycott et al. 2009'), values = c("#1b9e77", "#d95f02")) + 
  #viridis::scale_colour_viridis(discrete = TRUE, labels = c('Dunic et al.', 'Waycott et al. 2009')) + 
  #scale_shape_manual(labels = c('Dunic et al.', 'Waycott et al. 2009'), values = c(21, 21)) + 
  labs(size = expression(paste("Maximum observed area", " (km"^2, ")", sep = ""))) #+ 
  #ylim(c(-100, 200))  # You lose too much detail when putting everything on the
  # same y-axis scale; variability in the TNAW outweighs everything else
}

status_tnp <- 
  mk_status_plot(df = decadal_df_ne_1930 %>% filter(bioregion == "Temperate North Pacific"), 
                 point_colour = !!bioregion_point_colour, size_bins = max_area_km2) %>% 
  custom_layers() + 
  geom_point(data = tibble(waycott_study = FALSE, decade = c(1930, 1940), trend_ci_50 = c(0, 0), frac_of_max = c(1, 1)), alpha = 0) +
  ggtitle("Temperate North Pacific") + 
  theme(strip.text = element_text(size = 8))

status_tnaw <- 
  mk_status_plot(df = decadal_df_ne_1930 %>% filter(bioregion == "Temperate North Atlantic West"), 
                 point_colour = !!bioregion_point_colour, size_bins = max_area_km2) %>%
  custom_layers() +
  ggtitle("Temperate North Atlantic West")

status_tnae <- 
  mk_status_plot(df = decadal_df_ne_1930 %>% filter(bioregion == "Temperate North Atlantic East"), 
                 point_colour = !!bioregion_point_colour, size_bins = max_area_km2) %>% 
  custom_layers() +
  ggtitle("Temperate North Atlantic East")

status_m <- 
  mk_status_plot(df = decadal_df_ne_1930 %>% filter(bioregion == "Mediterranean"), 
                 point_colour = !!bioregion_point_colour, size_bins = max_area_km2) %>%
    custom_layers() +
    ggtitle("Mediterranean") + 
    ylab("Annual rate of change (% per year)") + 
    theme(axis.title.y = element_text(size = 10, angle = 90, vjust = 2))

status_ta <- 
  mk_status_plot(df = decadal_df_ne_1930 %>% filter(bioregion == "Tropical Atlantic"), 
                 point_colour = !!bioregion_point_colour, size_bins = max_area_km2) %>%
  custom_layers() +
  ggtitle("Tropical Atlantic")

status_tip <- 
  mk_status_plot(df = decadal_df_ne_1930 %>% filter(bioregion == "Tropical Indo-Pacific"), 
                 point_colour = !!bioregion_point_colour, size_bins = max_area_km2) %>% 
    custom_layers() +
    ggtitle("Tropical Indo-Pacific")

status_so <- 
  mk_status_plot(df = decadal_df_ne_1930 %>% filter(bioregion == "Temperate Southern Oceans"), 
                 point_colour = !!bioregion_point_colour, size_bins = max_area_km2) %>% 
    custom_layers() +
    #guides(colour = "none", size = guide_legend(label.hjust = 0, override.aes = list(alpha = 1))) + 
    guides(colour = guide_legend(title = "Dataset", nrow = 2, title.position = 'top'), 
           size = guide_legend(title.position = 'top', label.hjust = 0, nrow = 1,
                               override.aes = list(alpha = 1))) + 
    theme(axis.title.x = element_text(size = 10, margin = margin(t = 5, r = 0, b = 0, l = 0)), 
          axis.text.x = element_text(size = 8), 
          legend.title = element_text(size = 10, vjust = 1, margin = margin(r = 10)),
          legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0),
          legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
          legend.key.width = unit(2, "mm"),
          legend.text = element_text(size = 10, margin = margin(r = 12, l = 0))) +
    ggtitle("Temperate Southern Oceans")

dev.off()
dev.new(width = 16.9 / cm(1), height = 24 / cm(1))
status_tnp + status_tnaw + status_tnae + status_m + status_ta + status_tip + status_so + 
  plot_layout(ncol = 1)
ggsave(here("figures/SOM/waycott_comparison_orange-green.png"), width = 16.9 / cm(1), height = 24 / cm(1))
