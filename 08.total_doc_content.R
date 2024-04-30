#--------------------------------------------------------------------------#
# Project: DOClimato
# Script purpose: Compute total DOC content
# Date: 29/04/2024
# Author: Thelma Panaiotis
#--------------------------------------------------------------------------#

source("utils.R")

## Load data ----
#--------------------------------------------------------------------------#

# List files of predictions
files <- list.files("data", pattern = "pred.Rdata", full.names = TRUE)

# Load predictod DOC
res <- sapply(files, function(x) mget(load(x)), simplify = TRUE) %>% 
  bind_rows() %>% 
  mutate(season = as.character(season))

# Keep only seasonal and stratified
res <- res %>% filter(season == "0" & cv_type == "stratified")


# Nice layer names
resp_to_layer <- tribble(
  ~resp, ~layer,
  "doc_surf", "Surf.",
  "doc_epi", "Epi.",
  "doc_meso", "Meso.",
  "doc_bathy", "Bathy."
) %>% 
  mutate(layer = fct_inorder(layer))


# DOC in µmol kg⁻¹
# conversion of volume of seawater to mass?
# constant: 1 027.7 kg m⁻³ 


## TODO
# - area & volume of pixels depends on depth
# - conversion from volume to kg of seawater: constant or compute it per pixel?
# - conversion from µmol to PgC of DOC
# - what about pixels with no predictions?
# overall this is an underestimate as we do not have all pixels.


## Prepare all Rsquares ----
#--------------------------------------------------------------------------#
# Rsquares
rsquares <- res %>% 
  select(resp, season, cv_type, fold, preds) %>% 
  unnest(preds) %>% 
  group_by(resp, season, cv_type, fold) %>% 
  mutate(truth = ifelse(
    !is.na(log_doc_surf), 
    log_doc_surf, 
    ifelse(
      !is.na(log_doc_epi),
      log_doc_epi,
      ifelse(
        !is.na(log_doc_meso),
        log_doc_meso,
        log_doc_bathy
      )
    )
  )) %>% 
  select(-contains("log_doc")) %>% 
  rsq(truth = truth, estimate = .pred) %>% 
  rename(rsq = .estimate)


## Read bathymetry data ----
#--------------------------------------------------------------------------#
# Read the csv file containing the bathymetry data
df_bat <- read_delim(
  "data/raw/bathy/marmap_coord_-180;-90;180;90_res_10.csv", 
  delim = ",",
  show_col_types = FALSE
) %>% 
  rename(lon = V1, lat = V2, bathy = V3) %>% 
  filter(bathy <= 0) %>% 
  mutate(bathy = abs(bathy))

# Round to match DOC & env pixels
df_bat <- df_bat %>% 
  mutate(
    lon = roundp(lon, precision = 1, f = floor) + 0.5,
    lat = roundp(lat, precision = 1, f = floor) + 0.5
  ) %>% 
  filter(between(lat, -90, 90) & between(lon, -180, 180)) %>% 
  group_by(lon, lat) %>% 
  summarise(bottom = mean(bathy, na.rm = TRUE), .groups = "drop")


## Compute pixel volume ----
#--------------------------------------------------------------------------#
# Earth radius
R <- 6378.137 # km # https://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.html

# From bathymetry, compute bottom of each layer, and set it to NA if layer does not exist
df_lay <- df_bat %>% 
  rowwise() %>% 
  mutate(
    bathy_bottom = ifelse(bottom <= 1000, NA, bottom),
    meso_bottom = ifelse(bottom <= 200, NA, min(1000, bottom)),
    epi_bottom = ifelse(bottom <= 10, NA, min(200, bottom)),
    surf_bottom = min(10, bottom),
  ) %>% 
  ungroup()

# Finally, compute range of each layer
df_lay <- df_lay %>% 
  mutate(
    surf_range = surf_bottom,
    epi_range = epi_bottom - surf_bottom,
    meso_range = meso_bottom - epi_bottom,
    bathy_range = bathy_bottom - meso_bottom
  )

# Compute pixel area and volume per layer
pixels <- df_lay %>% 
  select(lon, lat, bottom, contains("range")) %>% 
  pivot_longer(contains("range"), names_to = "layer", values_to = "range") %>% 
  mutate(
    layer = str_remove(layer, "_range"),
    lat1 = (lat - 0.5) * pi / 180, # lower border
    lat2 = (lat + 0.5) * pi / 180, # upper border
    area_m2 = ((pi/180)*(R)^2 * abs(sin(lat1) - sin(lat2)))*10^6, # area in m²
    area_km2 = area_m2 / 1e6, # area in km²
    vol_m3 = area_m2 * range, # vol in m³
    vol_km3 = area_km2 * (range/1000) # vol in km³
  ) 

ggplot(pixels) + geom_histogram(aes(x = range, fill = layer), alpha = 0.5, bins = 50) + scale_fill_brewer(palette = "Blues")
ggplot(pixels) + geom_histogram(aes(x = vol_km3, fill = layer), alpha = 0.5, bins = 50) + scale_fill_brewer(palette = "Blues")


## Compute total DOC a minima ----
#--------------------------------------------------------------------------#
# Get predictions per fold
doc_preds <- res %>% 
  select(resp, fold, new_preds) %>% # keep only relevant columns
  unnest(new_preds) %>% 
  mutate(pred_doc = exp(pred_doc_log), .before = pred_doc_log) %>%  # compute doc from log(doc)
  select(resp, fold, lon, lat, pred_doc) %>% 
  rename(layer = resp) %>% 
  mutate(layer = str_remove(layer, "doc_"))

# Plot all layers for Fold01 (N.B. folds are not related across layers)
doc_preds %>% 
  filter(fold == "Fold01") %>% 
  ggplot() + 
  geom_polygon(data = world, aes(x = lon, y = lat, group = group), fill = "gray") +
  geom_raster(aes(x = lon, y = lat, fill = pred_doc)) +
  scale_fill_cmocean(name = "matter", direction = -1, trans = "log1p") +
  coord_quickmap(expand = 0) +
  facet_wrap(~layer)


# Join with pixel volumes
tot_doc <- doc_preds %>% 
  left_join(pixels, by = join_by(layer, lon, lat)) %>% 
  # Compute total DOC
  mutate(
    mass_kg = vol_m3 * 1027.7, # mass of seawater from volume
    tot_doc_mmol = pred_doc * mass_kg, # total DOC in one pixel, in µmol
    tot_doc_gc = (tot_doc_mmol / 1e6) * 12, # g of C
    tot_doc_Pgc = tot_doc_gc / 1e15 # T of C
  ) %>% 
  # Sum DOC content per fold
  group_by(fold, layer) %>% 
  summarise(tot_doc_Pgc = sum(tot_doc_Pgc, na.rm = TRUE), .groups = "drop")

tot_doc %>% 
  mutate(layer = factor(layer, levels = c("surf", "epi", "meso", "bathy"))) %>% 
  ggplot() + 
  geom_boxplot(aes(x = layer, y = tot_doc_Pgc, group = layer, colour = layer), show.legend = FALSE) +
  scale_colour_manual(values = c("#bdd7e7", "#6baed6", "#3182bd", "#08519c")) +
  facet_wrap(~layer, scales = "free", nrow = 1)


# Compute avg and sd of POC per layer
l_tot_doc <- tot_doc %>% 
  group_by(layer) %>% 
  summarise(
    tot_doc_avg = mean(tot_doc_Pgc),
    tot_doc_sd = sd(tot_doc_Pgc)
  )
l_tot_doc

# Compute avg and sd of total DOC
g_tot_doc <- tot_doc %>% 
  group_by(fold) %>% 
  summarise(tot_doc_Pgc = sum(tot_doc_Pgc, na.rm = TRUE), .groups = "drop") %>% 
  summarise(
    tot_doc_avg = mean(tot_doc_Pgc),
    tot_doc_sd = sd(tot_doc_Pgc)
  )

g_tot_doc$tot_doc_avg
g_tot_doc$tot_doc_sd
#  678.63 ± 0.36 PgC


## Compute total DOC with empty pixels replaced by mean ----
#--------------------------------------------------------------------------#
# Load pixels to include
load("data/00.bathymetry_inc.Rdata")

# Existing pixels in each layer (i.e. should be predicted if predictors were available)
to_pred <- df_inc %>% 
  pivot_longer(contains("include"), names_to = "layer", values_to = "include") %>% 
  mutate(layer = str_replace(layer, "_include", "")) %>% 
  drop_na()

# Fill NA with average of each fold in each layer
doc_preds_fill <- to_pred %>% 
  crossing(fold = unique(doc_preds$fold)) %>% 
  # Join with doc predictions
  left_join(doc_preds, by = join_by(lon, lat, layer, fold)) %>% 
  # In each layer and for each fold, replace NA by average value
  group_by(layer, fold) %>% 
  mutate(across(pred_doc, ~replace_na(., mean(., na.rm = TRUE)))) %>% 
  ungroup()

# Join with pixel volumes
tot_doc_fill <- doc_preds_fill %>% 
  left_join(pixels, by = join_by(layer, lon, lat)) %>% 
  # Compute total DOC
  mutate(
    mass_kg = vol_m3 * 1027.7, # mass of seawater from volume
    tot_doc_mmol = pred_doc * mass_kg, # total DOC in one pixel, in µmol
    tot_doc_gc = (tot_doc_mmol / 1e6) * 12, # g of C
    tot_doc_Pgc = tot_doc_gc / 1e15 # T of C
  ) %>% 
  # Sum DOC content per fold
  group_by(fold, layer) %>% 
  summarise(tot_doc_Pgc = sum(tot_doc_Pgc, na.rm = TRUE), .groups = "drop")

# Compute avg and sd of POC per layer
l_tot_doc_fill <- tot_doc_fill %>% 
  group_by(layer) %>% 
  summarise(
    tot_doc_avg = mean(tot_doc_Pgc),
    tot_doc_sd = sd(tot_doc_Pgc)
  )
l_tot_doc_fill

# Compute avg and sd of total DOC
g_tot_doc_fill <- tot_doc_fill %>% 
  group_by(fold) %>% 
  summarise(tot_doc_Pgc = sum(tot_doc_Pgc, na.rm = TRUE), .groups = "drop") %>% 
  summarise(
    tot_doc_avg = mean(tot_doc_Pgc),
    tot_doc_sd = sd(tot_doc_Pgc)
  )

g_tot_doc_fill$tot_doc_avg
g_tot_doc_fill$tot_doc_sd
#  690.71 ± 0.36 PgC

