library(sf)
library(viridis)
library(maptools)
#library(plotly)
library(here)

# SEAGRASS & SITES -------------------------------------------------------------
seagrass_points <- rgdal::readOGR(here::here('data/spatial-data-layers/WCMC013-014_SeagrassPtPy2020_v7/01_Data/WCMC_013_014_SeagrassesPt_v7.shp'))

points_df <- as_tibble(seagrass_points)

source(here('R/03_decadal_data_prep.R'))

map_df <- distinct(area_ts, study_site, .keep_all = TRUE)

# Get world map that has polygons that can be filled
# ------------------------------------------------------------------------------
# https://stackoverflow.com/questions/62059585/remove-countries-political-borders-from-ggplot2-map
# Defining a general CRS
mycrs <- "+proj=longlat +datum=WGS84 +no_defs"

# Using the original maps package, then converting map into SpatialPolygons object
world <- maps::map("world", fill=TRUE) %>% 
# maps::map("world", wrap = c(-100, 270), fill = TRUE, resolution = 0) %>% 
  maptools::map2SpatialPolygons(., IDs = .$names, proj4string = CRS(mycrs))

# The resulting map has self intersection problems so any further operation reports errors; using buffers of width 0 is a fast fix
while(rgeos::gIsValid(world)==FALSE){
  world <- rgeos::gBuffer(world, byid = TRUE, width = 0, quadsegs = 5, capStyle = "ROUND")
}

# Dissolving polygon's limits
world <- raster::aggregate(world)

# ------------------------------------------------------------------------------
# Add bioregions
# ------------------------------------------------------------------------------
short_bioregions <- read_sf(here::here('data/spatial-data-layers/SeagrassBioregions_Finalised/SeagrassBioregions_Finalised.shp'))

# Split TNAW & TNAE
# ------------------------------------------------------------------------------
# Make Atlantic dividing line at lon = -30:
line_df <- tibble(x = c(-30, -30), y = c(34, 90))
line_sf <- st_sf(id = 'L1', st_sfc(st_linestring(as.matrix(line_df), dim = "XY")))
st_crs(line_sf) <- st_crs(4326) # assign crs

# Split Atlantic and separate polygons
atlantics <- 
  filter(short_bioregions, Bioregion == "TemperateNorthAtlantic") %>% 
  lwgeom::st_split(., line_sf) %>% 
  st_collection_extract(., "POLYGON") %>%
  mutate(side = c("West", "East")) %>% 
  mutate(Bioregion = paste0(Bioregion, side)) %>% 
  distinct(Bioregion)

mod_bioregions <- 
  bind_rows(short_bioregions %>% filter(Bioregion != "TemperateNorthAtlantic"), 
            atlantics
  ) %>%
  left_join(.,
    bioregion_key %>% 
    mutate(bioregion = gsub("Indo-Pacific", "Indo Pacific", bioregion)) %>% 
    mutate(Bioregion = gsub(" ", "", bioregion)) %>% 
    mutate(Bioregion = gsub("Oceans", "Ocean", Bioregion))
  ) %>% 
  mutate(label = as.numeric(factor(bioregion, levels = bioregion_key$bioregion))) %>% 
  arrange(label) %>% 
  mutate(bioregion_colours = fct_reorder(bioregion_colours, bioregion))

lab_xy <- 
  textConnection(
    "bioregion, lat, lon,
    Temperate North Pacific, 45, -149
    Temperate North Atlantic West, 50, -40
    Temperate North Atlantic East, 50.5, -20
    Tropical Atlantic, 8, -32
    Mediterranean, 28, 1
    Tropical Indo-Pacific, -9, 83
    Temperate Southern Oceans, -40, -14"
  ) %>% 
  read.csv(stringsAsFactors = FALSE) %>% 
  mutate(bioregion = trimws(bioregion)) %>% 
  mutate(label = 1:7) %>% 
  select(-X) %>% 
  as_tibble

dev.off()
dev.new(width = 16.9 / cm(1), height = 10 / cm(1)); plot(1)
  ggplot() + 
    geom_sf(data = mod_bioregions, aes(fill = bioregion_colours), colour = NA, size = 0.1, alpha = 0.3) +
    scale_fill_identity() +
    geom_polygon(data = world, aes(x = long, y = lat, group = group), fill = 'grey95', colour = "grey50", size = 0.2) +
    geom_point(data = points_df, aes(x = coords.x1, y = coords.x2), colour = '#75D054FF', size = 1) + 
    geom_point(data = filter(map_df, waycott_study == TRUE), aes(x = lon, y = lat), 
               shape = 1, size = 2, stroke = 1, colour = '#440154FF') + 
    geom_point(data = filter(map_df, waycott_study == FALSE), aes(x = lon, y = lat), 
               shape = 21, size = 2, stroke = 0.25, colour = "grey30", fill = "#FDE725FF", 
               alpha = 0.8) + 
    geom_point(data = area_ts, aes(x = lon, y = lat, colour = waycott_study, shape = waycott_study), alpha = 0, size = 2, stroke = 1) + 
    geom_text(data = lab_xy, aes(x = lon, y = lat, label = label), colour = "grey50", size = 3.5) + 
    scale_colour_manual(values = c('#FDE725FF', '#440154FF'), labels = c("This study", "Waycott et al. 2009")) + 
    scale_shape_manual(values = c(19, 1)) + 
    map_theme(base_size = 12, background_fill = 'white') + 
    theme(plot.background = element_blank(),
          plot.margin = margin(0),
          legend.justification = c(0, 0), 
          #legend.position = c(0.04, 0.035), 
          legend.position = c(0.03, 0.03), 
          legend.direction = "horizontal",
          legend.title = element_blank(),
          legend.text = element_text(margin = margin(t = 0, b = 0, l = -2, unit = 'mm')),
          legend.margin = margin(0, 0, 0, 0),
          legend.background = element_rect(fill = "white", colour = NA)) +
    guides(fill = "none", colour = guide_legend(override.aes = list(shape = c(21, 1), colour = c('grey30', '#440154FF'), fill = c('#FDE725FF', "white"), size = c(4, 2.5), stroke = c(0.5, 1.5), alpha = 1)), 
           shape = "none") + 
    ylim(-70, 90)

ggsave(here::here('figures/site_map_with_bioregions.png'), width = 16.9 / cm(1), height = 10 / cm(1))
