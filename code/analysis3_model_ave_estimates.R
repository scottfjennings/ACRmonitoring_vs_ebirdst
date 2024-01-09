




library(tidyverse)
library(MASS)
library(AICcmodavg)
library(marginaleffects)
library(birdnames)
library(here)

options(scipen = 999)



#' mod_predictions_link
#' 
#' calculate model estimates (predictions) for a species and cell on the scale of the link function
#'
#' @param zspp_cell 4 letter species code pasted to cell code, e.g. GREG_All
#'
#' @return data frame
#' @export
#'
#' @examples
mod_avg_predictions <- function(zspp_cell) {
  
  spp_cell <- data.frame(spp.cell = zspp_cell)%>%
    separate(spp.cell, c("spp", "ebird.cell.num"), sep = "_") %>% 
    mutate(species = spp,
           ebird.cell.num = as.numeric(ebird.cell.num),
           spp = translate_bird_names(species, "alpha.code", "common.name")) %>% 
    left_join(readRDS(here("data/cell_mean_rain_lag")))
  
  znewdat <- readRDS(here("data/ebird_cell_annual")) %>% 
    filter(between(year, 2012, 2022)) %>% 
    arrange(ebird.cell.num, species, year) %>% 
    dplyr::select(year, ebird.cell.num, species, tot.nests) %>% 
    right_join(spp_cell %>% dplyr::select(species, ebird.cell.num, "wt.lag.cell.rain" = mean.cell.rain))
    
    #data.frame(year = seq(1995, 2019), wt.lag.cell.rain = spp_cell$mean.cell.rain)
  
  zmods <- readRDS(here("model_objects/spp_cell_mods"))[[zspp_cell]]
  
  zspp_pred <- modavgPred(zmods, modnames = names(zmods), newdata = znewdat)$matrix.output %>% 
    data.frame() %>% 
    bind_cols(znewdat)  %>%
    mutate(alpha.code = zspp_cell)
  
}


modeled_spp_cells <- readRDS(here("data/modeled_spp_cells")) %>%
  mutate(species = factor(species, levels = c("GREG", "GBHE", "SNEG", "BCNH"))) %>% 
  arrange(species, ebird.cell.num)



all_mod_av_preds <- map_df(modeled_spp_cells$spp.cell, mod_avg_predictions)

saveRDS(all_mod_av_preds, here("model_objects/all_mod_av_preds"))






# --


spp_cell <- data.frame(spp.cell = zspp_cell)%>%
  separate(spp.cell, c("spp", "ebird.cell.num"), sep = "_") %>% 
  mutate(species = spp,
         ebird.cell.num = as.numeric(ebird.cell.num),
         spp = translate_bird_names(species, "alpha.code", "common.name")) %>% 
  left_join(readRDS(here("data/cell_mean_rain_lag")))

znewdat <- readRDS(here("data/ebird_cell_annual")) %>% 
  filter(between(year, 2012, 2022)) %>% 
  arrange(ebird.cell.num, species, year) %>% 
  dplyr::select(year, ebird.cell.num, species, tot.nests) %>% 
  right_join(spp_cell %>% dplyr::select(species, ebird.cell.num, "wt.lag.cell.rain" = mean.cell.rain))

#data.frame(year = seq(1995, 2019), wt.lag.cell.rain = spp_cell$mean.cell.rain)

zmods <- readRDS(here("model_objects/spp_cell_mods"))[[zspp_cell]]

zmod = zmods$year2

zchange <- comparisons(zmod, newdata = znewdat2, comparison = "difference")





zspp_pred <- modavgPred(zmods, modnames = names(zmods), newdata = znewdat2)$matrix.output %>% 
  data.frame() %>% 
  bind_cols(znewdat)  %>%
  mutate(alpha.code = zspp_cell)


