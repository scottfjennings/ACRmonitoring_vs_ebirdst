

# fit set of negative binomial models to the HEP data with colonies grouped by ebirdst cells (IDed in analysis1_prep_data)


library(tidyverse)
library(here)
library(MASS)
library(AICcmodavg)
library(birdnames)

ebird_cell_annual <- readRDS(here("data/ebird_cell_annual")) %>% 
  filter(between(year, 2012, 2022)) %>% 
  arrange(ebird.cell.num, species, year)

cell_mean_rain_lag <- ebird_cell_annual %>% 
  group_by(ebird.cell.num, species) %>% 
  summarise(mean.cell.rain = mean(wt.lag.cell.rain)) 

saveRDS(cell_mean_rain_lag, here("data/cell_mean_rain_lag"))

zyears <- ebird_cell_annual %>% 
  data.frame() %>% 
  filter(tot.nests > 0, !is.na(ebird.cell.num)) %>% 
  group_by(ebird.cell.num, species) %>% 
  summarise(num.years = n())

spp_cells <- readRDS(here("data/spp_cells")) %>% 
  filter(species %in% c("GREG", "GBHE", "SNEG", "BCNH")) %>% 
  mutate(spp.cell = paste(species, ebird.cell.num, sep = "_")) %>%  
  left_join(zyears) %>% 
  arrange(species, ebird.cell.num) %>% 
  filter(num.years > 4, num.colonies > 2) %>% 
  left_join(cell_mean_rain_lag) %>% 
  filter(spp.cell != "SNEG_16")

saveRDS(spp_cells, here("data/modeled_spp_cells"))



#' Fit a set of negative binomial glm models on untransformed nest abundance for a species in a cellion
#'
#' @param zspp 4 letter code for species
#' @param zcell cellion code (use codes defined in "C:/Users/scott.jennings/OneDrive - Audubon Canyon Ranch/Projects/core_monitoring_research/HEP/HEP_data_work/HEP_data/cellion_key.csv")
#'
#' @return a list with an element for each fitted model object
#'
#' @examples
#' spp_cell_mods_glmnb <- map2(spp_cell$species, spp_cell$cellion, fit_multi_mods_glmbn)
#' names(spp_cell_mods_glmnb) <- spp_cell$spp.cell
fit_multi_mods_glmbn <- function(zspp, zebird.cell) {
  
#  zspp = spp_cells$species[1]
#  zebird.cell = spp_cells$ebird.cell.num[1]
  
  zdat <- ebird_cell_annual %>% 
    filter(species == zspp, ebird.cell.num == zebird.cell)
  
  
  zmods = list(
    "year2_rain" = MASS::glm.nb(data = zdat, formula = tot.nests ~ poly(year, 2) + wt.lag.cell.rain),
    "year_rain" = MASS::glm.nb(data = zdat, formula = tot.nests ~ year + wt.lag.cell.rain),
    #"lnyear_rain" = MASS::glm.nb(data = zdat, formula = tot.nests ~ log(year - 1994) + wt.lag.cell.rain),
    
    "year2_rain2" = MASS::glm.nb(data = zdat, formula = tot.nests ~ poly(year, 2) + poly(wt.lag.cell.rain, 2)),
    "year_rain2" = MASS::glm.nb(data = zdat, formula = tot.nests ~ year + poly(wt.lag.cell.rain, 2)),
    #"lnyear_rain2" = MASS::glm.nb(data = zdat, formula = tot.nests ~ log(year - 1994) + poly(wt.lag.cell.rain, 2)),
    
    "year2" = MASS::glm.nb(data = zdat, formula = tot.nests ~ poly(year, 2)),
    "year" = MASS::glm.nb(data = zdat, formula = tot.nests ~ year),
    #"lnyear" = MASS::glm.nb(data = zdat, formula = tot.nests ~ log(year - 1994)),
    "rain" = MASS::glm.nb(data = zdat, formula = tot.nests ~ wt.lag.cell.rain),
    "rain2" = MASS::glm.nb(data = zdat, formula = tot.nests ~ poly(wt.lag.cell.rain, 2)),
    
    "intercept" = MASS::glm.nb(data = zdat, formula = tot.nests ~ 1)
  )
}




spp_cell_mods <- map2(spp_cells$species, spp_cells$ebird.cell.num, fit_multi_mods_glmbn)
names(spp_cell_mods) <- spp_cells$spp.cell

saveRDS(spp_cell_mods, here("model_objects/spp_cell_mods"))


