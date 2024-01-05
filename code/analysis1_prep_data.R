



library(tidyverse)
library(here)
library(birdnames)
library(ebirdst)
library(geosphere)
library(sf)
library(rnaturalearth)

# The main data prep is to assign each HEP colony to the appropriate eBird cell

# boundaries of states in the united states
states <- ne_states(iso_a2 = "US", returnclass = "sf") %>%
  filter(iso_a2 == "US", postal %in% c("CA")) %>%
  transmute(state = iso_3166_2)

# 1 filter eBird estimates to just our study area ----
# load eBirdst estimates
# bay area bounding box to filter trend estimates
study_area_bbox <- data.frame(lat.max = 39, lat.min = 37.5, lon.max = -121.6, lon.min = -123.5)
# eBirdst products use 6-7 character species codes. this helper file adds 4 letter codes
species_codes <- data.frame(species_code = c( "greegr", "grbher3", "snoegr"),
                            alpha.code = c("GREG", "GBHE", "SNEG"))

# load data, the location these are saved to is determined by ebirdst; load_trends() knows where that location is
sfba_trends <- map_df(c("Great Egret", "Great Blue Heron", "Snowy Egret"), load_trends) %>% 
  filter(between(latitude, study_area_bbox$lat.min, study_area_bbox$lat.max) & between(longitude, study_area_bbox$lon.min, study_area_bbox$lon.max)) %>% 
  full_join(species_codes)

# plot to double check
ggplot(states) +
  geom_sf() +
  geom_point(data = sfba_trends, aes(x = longitude, y = latitude, color = abd_trend)) +
  facet_wrap(~alpha.code) +
  ylim(37.5, 39) +
  xlim(-123.5, -121.6)

saveRDS(sfba_trends, here("data/ebird_sfba_trends"))

# 2 create a unique id for each eBird cell ----
# ebirdst products are calculated for each 27x27 km cell in N America
# need to assign monitoring data to those cells, and thus need a helper ID for each cell
# also create the 27km squares around each cell center; easiest to do this in UTM

ebird_cells <- distinct(sfba_trends, longitude, latitude) %>% 
  arrange(longitude, latitude) %>% 
  mutate(ebird.cell.num = row_number()) %>% 
  ungroup()

# create lat/long sf object
ebird_p1 <- st_as_sf(ebird_cells, coords = c("longitude", "latitude"), crs = "EPSG:4326")
# convert to utm
ebird_p2 <- st_transform(ebird_p1, crs= "EPSG:32610")
# create 27km boxes
ebird_cells_utm <- ebird_p2 %>%
  distinct(ebird.cell.num, geometry) %>% 
  dplyr::mutate(ebird.utmeast = sf::st_coordinates(.)[,1],
                ebird.utmnorth = sf::st_coordinates(.)[,2],
                max.utmeast = ebird.utmeast + 27000/2,
                max.utmnorth = ebird.utmnorth + 27000/2,
                min.utmeast = ebird.utmeast - 27000/2,
                min.utmnorth = ebird.utmnorth - 27000/2) %>% 
  data.frame() %>% 
  select(-geometry)

# bind lat/long with utm
ebird_ll_utm <- ebird_cells_utm %>% 
  select(contains("ebird")) %>% 
  full_join(ebird_cells)

# then output cell id with utm and lat/lon
saveRDS(ebird_ll_utm, here("data/ebird_ll_utm"))



# 3. assign HEP colonies to appropriate eBird cells ----

# procedure is to draw a 27x27 km box around each cell center, then id hep sites inside each box
# but need to convert back to lat/long for plotting over rnaturalearth basemap

max_ebird_cell_utm <- ebird_cells_utm %>% 
  select(ebird.cell.num, max.utmeast, max.utmnorth) %>% 
  st_as_sf(coords = c("max.utmeast", "max.utmnorth"), crs = "EPSG:32610")

max_ebird_cell_ll <- st_transform(max_ebird_cell_utm, crs= "EPSG:4326") %>% 
  dplyr::mutate(max.longitude = sf::st_coordinates(.)[,1],
                max.latitude = sf::st_coordinates(.)[,2]) %>% 
  data.frame() %>% 
  select(-geometry, -contains("utm"))


min_ebird_cell_utm <- ebird_cells_utm %>% 
  select(ebird.cell.num, min.utmeast, min.utmnorth) %>% 
  st_as_sf(coords = c("min.utmeast", "min.utmnorth"), crs = "EPSG:32610")

min_ebird_cell_ll <- st_transform(min_ebird_cell_utm, crs= "EPSG:4326") %>% 
  dplyr::mutate(min.longitude = sf::st_coordinates(.)[,1],
                min.latitude = sf::st_coordinates(.)[,2]) %>% 
  data.frame() %>% 
  select(-geometry, -contains("utm"))


max_min_ebird_cell_ll <- full_join(min_ebird_cell_ll, max_ebird_cell_ll)

# get list of HEP sites and codes
source("https://raw.githubusercontent.com/scottfjennings/HEP_data_work/master/HEP_code/HEP_utility_functions.R")
hep_sites <- hep_sites_from_access(here("C:/Users/scott.jennings/OneDrive - Audubon Canyon Ranch/Projects/core_monitoring_research/HEP/HEP_data_work/HEP_data/HEPDATA.accdb")) %>% 
  select(code, utmnorth, utmeast) %>% 
  filter(!is.na(utmnorth))

#' hep_site_ebirder
#' 
#' check whether a HEP colony is in an eBird cell
#'
#' @param zebird.cell.num 
#'
#' @return data frame with ebird cell id number if the colony is within that ebird cell
#'
#' @examples
hep_site_ebirder <- function(zebird.cell.num) {
  zebird_cells <- filter(ebird_cells_utm, ebird.cell.num == zebird.cell.num)
  # there's probably a spatial masking function to do this, but paired between() calls works fine
ebird_hep_sites <- hep_sites %>% 
  select(code, utmnorth, utmeast) %>% 
  bind_cols(zebird_cells) %>% 
  filter(between(utmnorth, min.utmnorth, max.utmnorth) &
                                   between(utmeast, min.utmeast, max.utmeast)) %>% 
  select(-contains("min"), -contains("max"))
}


all_hep_ebird <- map_df(ebird_cells$ebird.cell.num, hep_site_ebirder) %>% 
  full_join(hep_sites)


# hep sites are in utm, convert to lat/long for plotting over rnaturalearth base map
hep_p1 <- st_as_sf(all_hep_ebird, coords = c("utmeast", "utmnorth"), crs = "EPSG:32610")
hep_p2 <- st_transform(hep_p1, crs= "EPSG:4326") %>% 
  dplyr::mutate(longitude = sf::st_coordinates(.)[,1],
                latitude = sf::st_coordinates(.)[,2]) %>% 
  data.frame() %>% 
  select(-geometry)
  

plotting_df <- hep_p2 %>% 
  select(code, longitude, latitude, ebird.cell.num) %>% 
  left_join(max_min_ebird_cell_ll)

# plot to double check cell assignments
ggplot(states) +
  geom_sf() +
  geom_point(data = plotting_df, aes(x = longitude, y = latitude, color = as.factor(ebird.cell.num))) +
  #  geom_point(aes(x = ebird.utmeast, y = ebird.utmnorth, color = as.factor(ebird.cell.num)), size = 3, shape = 2) +
  geom_rect(data = plotting_df, aes(xmin = min.longitude, xmax = max.longitude, ymin = min.latitude, ymax = max.latitude, color = as.factor(ebird.cell.num)), fill = NA) +
  ylim(37.5, 39) +
  xlim(-123.5, -121.6)

# note Bolinas colonies don't have an ebird cell

# and save if everything looks right
saveRDS(plotting_df, here("data/plotting_df"))





