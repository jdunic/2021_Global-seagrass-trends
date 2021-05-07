# Viewing
# ------------------------------------------------------------------------------
# My alternative for View() using table in browser
view <- function(object, theme = "scientific", digits = 3) {
  if((length(digits) == 1) & inherits(object, "data.frame")) {
    object <- dplyr::mutate_if(object, is.numeric, round, digits)
  }
  tableHTML::add_theme(tableHTML::tableHTML(object), theme)
}

# Interpolation functions
# -----------------------
do_na_interpolation <- function(site_df, option = "linear", col = "area", new_col = "inter_area") {
  site_df[[new_col]] <- na_interpolation(site_df[[col]], option = option)
  return(site_df)
}

get_na_stats <- function(site_df, col, print = FALSE) {
  statsNA(site_df[[col]], print = print)
}

# Generate MVN for use in lp %*% br
rmvn <- function(n, mu, sig) { ## MVN random deviates

  L <- mroot(sig)
  m <- ncol(L)
  t(mu + L %*%  matrix(rnorm(m * n), m, n))
}

# Get predictions and summary stats for meadows for any year
get_yearly_preds <- function(model_obj = m1, newdata, n_draws = 1000,
                             ci_probs = c(0.025, 0.5, 0.975)) {
  # linear predictors
  lp1 <- predict(model_obj, newdata = newdata, type = "lpmatrix")
  # n_draws of the parameters from covariance matrix
  br <- rmvn(n_draws, coef(model_obj), model_obj$Vp)
#
  pred_area <- matrix(NA, nrow = nrow(newdata), ncol = n_draws)
  for (i in 1:n_draws) {
    pred_area[, i] <- lp1 %*% br[i, ]
  }
#
  area_ci <- as_tibble(t(apply(pred_area, 1, quantile, probs = ci_probs)))
  names(area_ci) <- paste("log_area", names(area_ci), sep = "_")
  dout <- bind_cols(newdata, area_ci) %>% as_tibble
#
  return(dout)
}

floor_decade <- function(value) { return(value - value %% 10) }

# Returns quantile df
get_bioregion_quantiles <- function(gam_obj) {
  newdata <- expand.grid(year = seq(min(gam_obj$model$year), max(gam_obj$model$year), by = 1),
                      study_site = gam_obj$model$study_site[1],
                      # We just use any study_site name here,
                      # as below we predict ignoring this.
                      # We need to include them when we do decadal areas
                     bioregion = unique(gam_obj$model$bioregion))
  lp <- predict(gam_obj, exclude = "study_site", newdata = newdata, type = "lpmatrix")
  lpnames <- colnames(lp)
  ibioregion <- grep("bioregion", lpnames)
  coefnames <- names(coef(gam_obj))
  ibioregion_coef <- grep("bioregion", coefnames)
  br <- rmvn(n = 1000, mu = coef(gam_obj)[ibioregion_coef], 
             sig = gam_obj$Vp[ibioregion_coef, ibioregion_coef])
#
  post <- matrix(NA, nrow = nrow(newdata), ncol = 1000)
#  
  for (i in 1:1000){
    post[, i] <- (lp[, ibioregion] %*% br[i, ])
  }
  CIdat <- t(apply(post, 1, quantile, probs = c(0.025, 0.5, 0.975))) %>%
    data.frame(newdata, .) %>% 
    as_tibble()
  return(CIdat)
}


# Make status plots
#------------------------------------------------------------------------------

mk_status_plot <- function(df, point_colour = studyid, size_bins = max_area_km2, event_timeline = NULL) {
  point_colour <- enquo(point_colour)
  size_bins <- enquo(size_bins)
  event_timeline <- enquo(event_timeline)
  ggplot(data = df) + 
  aes(x = frac_of_max, y = trend_ci_50, colour = !!point_colour) +
  geom_point(aes(size = !!size_bins, shape = !!event_timeline), alpha = 0.4) +
  geom_hline(yintercept = 0, colour = "grey60") + 
  geom_vline(xintercept = 1, linetype = 1, colour = "grey60") + 
  labs(x = "Fraction of maximum observed area", 
       y = "Annual rate of change (%)") +
  mytheme(base_size = 12) + 
  theme(legend.position = "bottom", 
        axis.text.x = element_text(size = 12), 
        panel.grid.major.y = element_line(colour = "grey95"),
        panel.background = element_rect(fill = "grey98", colour = NA)) + 
  scale_x_continuous(breaks = seq(0, 1.3, 0.1), 
                     labels = c(0, '', '', '', '', 0.5, '', '', '', '', 1, '', '', ''), 
                     lim = c(0, 1.3)) + 
  scale_y_continuous(trans = pseudo_log_trans(sigma = 0.5),
                     breaks = c(-200, -50, -10, -2, 0, 2, 10, 50, 200)) + 
  guides(size = guide_legend(override.aes = list(alpha = 1)))
}


# ------------------------------------------------------------------------------
# PLOTTING FUNCTIONS 
# ------------------------------------------------------------------------------

# Themes
# ------------------------------------------------------------------------------
# See: https://rpubs.com/Koundy/71792 for inspiration
mytheme <- function(base_size = 18, base_family = "helvetica") {
  theme_minimal(base_size = base_size) + 
  theme(axis.line = element_line(colour = "black"), 
        panel.grid = element_blank(), 
        axis.title.y = element_text(angle = 90, vjust = 2),
        axis.title.x = element_text(vjust = -5), 
        plot.margin = unit(c(10, 5, 10, 5), "mm")
  )
}

map_theme <- function(base_size = 18, background_fill = "#f5f5f2") {
  theme_minimal(base_size = base_size) + 
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        plot.background = element_rect(fill = "#f5f5f2", color = NA), 
        panel.grid = element_blank(), 
        panel.background = element_rect(fill = background_fill, color = NA), 
        legend.background = element_rect(fill = "#f5f5f2", color = NA))
}


# https://stackoverflow.com/questions/23901907/create-a-log-sequence-across-multiple-orders-of-magnitude
lseq <- function(from=1, to=100000, length.out=10) {
  # logarithmic spaced sequence
  # blatantly stolen from library("emdbook"), because need only this
  exp(seq(log(from), log(to), length.out = length.out))
}

# https://fishandwhistle.net/post/2018/modifying-facet-scales-in-ggplot2/
# set unique scales for each axis
scale_override <- function(which, scale) {
  if(!is.numeric(which) || (length(which) != 1) || (which %% 1 != 0)) {
    stop("which must be an integer of length 1")
  }
  
  if(is.null(scale$aesthetics) || !any(c("x", "y") %in% scale$aesthetics)) {
    stop("scale must be an x or y position scale")
  }
  
  structure(list(which = which, scale = scale), class = "scale_override")
}

CustomFacetWrap <- ggplot2::ggproto(
  "CustomFacetWrap", ggplot2::FacetWrap,
  init_scales = function(self, layout, x_scale = NULL, y_scale = NULL, params) {
    # make the initial x, y scales list
    scales <- ggproto_parent(FacetWrap, self)$init_scales(layout, x_scale, y_scale, params)
    
    if(is.null(params$scale_overrides)) return(scales)
    
    max_scale_x <- length(scales$x)
    max_scale_y <- length(scales$y)
    
    # ... do some modification of the scales$x and scales$y here based on params$scale_overrides
    for(scale_override in params$scale_overrides) {
      which <- scale_override$which
      scale <- scale_override$scale
      
      if("x" %in% scale$aesthetics) {
        if(!is.null(scales$x)) {
          if(which < 0 || which > max_scale_x) stop("Invalid index of x scale: ", which)
          scales$x[[which]] <- scale$clone()
        }
      } else if("y" %in% scale$aesthetics) {
        if(!is.null(scales$y)) {
          if(which < 0 || which > max_scale_y) stop("Invalid index of y scale: ", which)
          scales$y[[which]] <- scale$clone()
        }
      } else {
        stop("Invalid scale")
      }
    }
    
    # return scales
    scales
  }
)

facet_wrap_custom <- function(..., scale_overrides = NULL) {
  # take advantage of the sanitizing that happens in facet_wrap
  facet_super <- facet_wrap(...)
  
  # sanitize scale overrides
  if(inherits(scale_overrides, "scale_override")) {
    scale_overrides <- list(scale_overrides)
  } else if(!is.list(scale_overrides) || 
            !all(vapply(scale_overrides, inherits, "scale_override", FUN.VALUE = logical(1)))) {
    stop("scale_overrides must be a scale_override object or a list of scale_override objects")
  }
  
  facet_super$params$scale_overrides <- scale_overrides
  
  ggplot2::ggproto(NULL, CustomFacetWrap,
    shrink = facet_super$shrink,
    params = facet_super$params
  )
}