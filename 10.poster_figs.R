#--------------------------------------------------------------------------#
# Project: DOClimato
# Script purpose: Prepare poster figures
# Date: 12/06/2023
# Author: Thelma Panaïotis
#--------------------------------------------------------------------------#


source("utils.R")
library(ggtext)
library(patchwork)

## Load data ----
#--------------------------------------------------------------------------#

# List files of predictions
files <- list.files("data", pattern = "pred.Rdata", full.names = TRUE)

# load data
res <- sapply(files, function(x) mget(load(x)), simplify = TRUE) %>% 
  bind_rows() %>% 
  mutate(season = as.character(season))

# Nice layer names
resp_to_layer <- tribble(
  ~resp, ~layer,
  "doc_surf", "Surf.",
  "doc_epi", "Epi.",
  "doc_meso", "Meso.",
  "doc_bathy", "Bathy."
) %>% 
  mutate(layer = fct_inorder(layer))


# Theme set-up
p_theme <- theme(
  legend.title = element_markdown(size = 8), 
  text = element_text(size = 8),
  #axis.title.x = element_text(size = 6),
  #axis.title.y = element_text(size = 6), 
  plot.margin = margin(0, 0, 0, 0, "pt"),
  legend.key.height = unit(0.4, "cm"),
  legend.key.width = unit(0.3, "cm")
) 

# World map centered on Pacific
# work around for polygons crossing the 0 latitude
polygon <- st_polygon(x = list(rbind(
  c(-0.0001, 90),
  c(0, 90),
  c(0, -90),
  c(-0.0001, -90),
  c(-0.0001, 90)
))) %>%
  st_sfc() %>%
  st_set_crs(4326)

# perform transformation on modified version of world dataset
world_rob <- st_transform(world_sf %>% st_difference(polygon), crs = '+proj=robin +lon_0=180 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs')

# Graticules
grat <- list(
  st_sf(geometry = st_sfc(
    st_multilinestring(x = lapply(c(0, 60, 120, 180, 240, 300, 359.9), 
                                  function(x) cbind(x, seq(-90, 90, 1)))), crs = 'WGS84')),
  st_sf(geometry = st_sfc(
    st_multilinestring(x = lapply(c(-90, -60, -30, 0, 30, 60, 90), function(x) {
      cbind(seq(0, 360, 1), x)
    })), crs = 'WGS84'))) %>% 
  bind_rows() %>% 
  st_transform(crs = '+proj=robin +lon_0=180 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs')


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


# Nice table of rsquares
rsquares %>% 
  filter(season == "0") %>% 
  filter(cv_type == "stratified") %>% 
  #filter(cv_type == "spatial") %>% 
  group_by(resp) %>% 
  summarise(
    mean = mean(rsq),
    sd = sd(rsq)
  )


## Map of annual observations ----
#--------------------------------------------------------------------------#
load("data/02.ann_surf.Rdata")
load("data/02.ann_epi.Rdata")
load("data/02.ann_meso.Rdata")
load("data/02.ann_bathy.Rdata")

p_obs_theme <- theme(
  axis.title = element_blank(),
  text = element_text(size = 8),
  plot.margin = unit(c(0, 0, 0, 0), "mm")
)

# All layers together
df_obs <- df_ann_surf_fit %>% 
  bind_rows(df_ann_epi_fit) %>% 
  bind_rows(df_ann_meso_fit) %>% 
  bind_rows(df_ann_bathy_fit)

p_obs <- ggplot(df_obs) +
  geom_sf(data = grat, alpha = 0.8, color = "gray80", linewidth = 0.2) +
  geom_point(aes(x = lon, y = lat), size = 0.05, alpha = 0.3) +
  geom_sf(data = world_rob, fill = "gray80", colour = NA) +
  coord_sf(crs = '+proj=robin +lon_0=180 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs', default_crs = sf::st_crs(4326), datum = NA) +
  p_obs_theme

## Save
ggsave(file = "figures/poster/obs.png", plot = p_obs, width = 120, height = 64, unit = "mm", bg = "white")


## Surface climatology ----
#--------------------------------------------------------------------------#
# Surface climatology:
# a - annual climatology
# b - uncertainty

p_surf_theme <- theme(
  axis.title = element_blank(),
  legend.position = "bottom",
  legend.title = element_markdown(size = 10),
  legend.key.width = unit(1.5, "cm"),
  legend.key.height = unit(2, "mm"),
  legend.direction = "horizontal",
  legend.title.position = "top",
  legend.margin = margin(0, 0, 0, 0),
  legend.box.margin = margin(-10, -10, 0, -10),
  plot.margin = unit(c(-10, -10, -10, -10), "mm"),
  text = element_text(size = 8)
  )

## a + b
# Get projections, average and sd by pixel
df_surf <- res %>% 
  filter(resp == "doc_surf" & season == "0") %>% 
  filter(cv_type == "stratified") %>% 
  select(fold, new_preds) %>% 
  unnest(new_preds) %>% 
  # Apply exp to predictions as we predicted log(doc)
  mutate(pred_doc = exp(pred_doc_log), .after = pred_doc_log) %>% 
  select(fold, lon, lat, contains("doc")) %>% 
  left_join(rsquares %>% filter(resp == "doc_surf" & season == "0" & cv_type == "stratified"), by = join_by(fold)) %>% 
  group_by(lon, lat) %>% 
  summarise(
    doc_avg = wtd.mean(pred_doc, weights = rsq, na.rm = TRUE), 
    doc_sd = sqrt(wtd.var(pred_doc, weights = rsq, na.rm = TRUE)), 
    .groups = "drop"
  )

# Plot a
# Colour bar limit
doc_lims <- c(min(df_surf$doc_avg), 200)
p_surf_a <- ggplot(df_surf %>% filter(lon != -0.5)) + # Need to remove lon = -0.5
  geom_sf(data = grat, alpha = 0.8, color = "gray80", linewidth = 0.2) +
  geom_tile(aes(x = lon, y = lat, fill = doc_avg, colour = doc_avg)) +
  geom_sf(data = world_rob, fill = "gray80", colour = NA) +
  ggplot2::scale_fill_viridis_c(option = "F", trans = "log1p", limits = doc_lims) +
  ggplot2::scale_colour_viridis_c(option = "F", trans = "log1p", guide = "none", limits = doc_lims) +
  coord_sf(datum = NA, default_crs = sf::st_crs(4326)) +
  labs(fill = "Predicted DOC (μmol kg<sup>-1</sup>)") +
  p_surf_theme

# Plot b
p_surf_b <- ggplot(df_surf  %>% filter(lon != -0.5))  +
  geom_sf(data = grat, alpha = 0.8, color = "gray80", linewidth = 0.2) +
  geom_tile(aes(x = lon, y = lat, fill = doc_sd, colour = doc_sd)) +
  geom_sf(data = world_rob, fill = "gray80", colour = NA) +
  ggplot2::scale_fill_viridis_c(option = "E", trans = "log1p", breaks = c(1, 5, 10, 20, 30)) +
  ggplot2::scale_colour_viridis_c(option = "E", trans = "log1p", guide = "none") +
  coord_sf(crs = '+proj=robin +lon_0=180 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs', default_crs = sf::st_crs(4326), datum = NA) +
  labs(fill = "Prediction uncertainty (μmol kg<sup>-1</sup>)") +
  p_surf_theme 


## Assemble
p_surf <- p_surf_a / p_surf_b + plot_layout(axis_titles = "collect")
#p_surf

# Save
ggsave(file = "figures/poster/surf_doc.png", plot = p_surf, width = 120, height = 150, unit = "mm", bg = "white")


## Deeper climatologies ----
#--------------------------------------------------------------------------#
# Deeper climatologies:
# a - epi climatology
# b - meso climatology
# c - bathy climatology

p_deep_theme <- theme(
  axis.title = element_blank(),
  legend.position = "bottom",
  legend.title = element_markdown(size = 10),
  legend.key.width = unit(1.2, "cm"),
  legend.key.height = unit(2, "mm"),
  legend.direction = "horizontal",
  legend.title.position = "top",
  legend.margin = margin(0, 0, 0, 0),
  legend.box.margin = margin(-10, -10, 0, -10),
  plot.margin = unit(c(0, 0, 0, 0), "mm"),
  text = element_text(size = 8)
)

## a, epipelagic predictions
df_deep_a <- res %>% 
  filter(resp == "doc_epi" & season == "0") %>% 
  filter(cv_type == "stratified") %>% 
  select(fold, new_preds) %>% 
  unnest(new_preds) %>% 
  # Apply exp to predictions as we predicted log(doc)
  mutate(pred_doc = exp(pred_doc_log), .after = pred_doc_log) %>% 
  select(fold, lon, lat, contains("doc")) %>% 
  left_join(rsquares %>% filter(resp == "doc_epi" & season == "0" & cv_type == "stratified"), by = join_by(fold)) %>% 
  group_by(lon, lat) %>% 
  summarise(
    doc_avg = wtd.mean(pred_doc, weights = rsq, na.rm = TRUE), 
    doc_sd = sqrt(wtd.var(pred_doc, weights = rsq, na.rm = TRUE)), 
    .groups = "drop"
  )

## b, mesopelagic predictions
df_deep_b <- res %>% 
  filter(resp == "doc_meso" & season == "0") %>% 
  filter(cv_type == "stratified") %>% 
  select(fold, new_preds) %>% 
  unnest(new_preds) %>% 
  # Apply exp to predictions as we predicted log(doc)
  mutate(pred_doc = exp(pred_doc_log), .after = pred_doc_log) %>% 
  select(fold, lon, lat, contains("doc")) %>% 
  left_join(rsquares %>% filter(resp == "doc_meso" & season == "0" & cv_type == "stratified"), by = join_by(fold)) %>% 
  group_by(lon, lat) %>% 
  summarise(
    doc_avg = wtd.mean(pred_doc, weights = rsq, na.rm = TRUE), 
    doc_sd = sqrt(wtd.var(pred_doc, weights = rsq, na.rm = TRUE)), 
    .groups = "drop"
  )

## c, bathypelagic predictions
df_deep_c <- res %>% 
  filter(resp == "doc_bathy" & season == "0") %>% 
  filter(cv_type == "stratified") %>% 
  select(fold, new_preds) %>% 
  unnest(new_preds) %>% 
  # Apply exp to predictions as we predicted log(doc)
  mutate(pred_doc = exp(pred_doc_log), .after = pred_doc_log) %>% 
  select(fold, lon, lat, contains("doc")) %>% 
  left_join(rsquares %>% filter(resp == "doc_bathy" & season == "0" & cv_type == "stratified"), by = join_by(fold)) %>% 
  group_by(lon, lat) %>% 
  summarise(
    doc_avg = wtd.mean(pred_doc, weights = rsq, na.rm = TRUE), 
    doc_sd = sqrt(wtd.var(pred_doc, weights = rsq, na.rm = TRUE)), 
    .groups = "drop"
  )

# Plot a
p_deep_a <- ggplot(df_deep_a %>% filter(lon != -0.5))  +
  geom_sf(data = grat, alpha = 0.8, color = "gray80", linewidth = 0.2) +
  geom_tile(aes(x = lon, y = lat, fill = doc_avg, colour = doc_avg)) +
  geom_sf(data = world_rob, fill = "gray80", colour = NA) +
  ggplot2::scale_fill_viridis_c(option = "F", trans = "log1p") +
  ggplot2::scale_colour_viridis_c(option = "F", trans = "log1p", guide = "none") +
  coord_sf(crs = '+proj=robin +lon_0=180 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs', default_crs = sf::st_crs(4326), datum = NA) +
  labs(fill = "Predicted DOC (μmol kg<sup>-1</sup>)") +
  p_deep_theme

# Plot b
p_deep_b <- ggplot(df_deep_b %>% filter(lon != -0.5))  +
  geom_sf(data = grat, alpha = 0.8, color = "gray80", linewidth = 0.2) +
  geom_tile(aes(x = lon, y = lat, fill = doc_avg, colour = doc_avg)) +
  geom_sf(data = world_rob, fill = "gray80", colour = NA) +
  ggplot2::scale_fill_viridis_c(option = "F", trans = "log1p") +
  ggplot2::scale_colour_viridis_c(option = "F", trans = "log1p", guide = "none") +
  coord_sf(crs = '+proj=robin +lon_0=180 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs', default_crs = sf::st_crs(4326), datum = NA) +
  labs(fill = "Predicted DOC (μmol kg<sup>-1</sup>)") +
  p_deep_theme

# Plot c
p_deep_c <- ggplot(df_deep_c %>% filter(lon != -0.5))  +
  geom_sf(data = grat, alpha = 0.8, color = "gray80", linewidth = 0.2) +
  geom_tile(aes(x = lon, y = lat, fill = doc_avg, colour = doc_avg)) +
  geom_sf(data = world_rob, fill = "gray80", colour = NA) +
  ggplot2::scale_fill_viridis_c(option = "F", trans = "log1p") +
  ggplot2::scale_colour_viridis_c(option = "F", trans = "log1p", guide = "none") +
  coord_sf(crs = '+proj=robin +lon_0=180 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs', default_crs = sf::st_crs(4326), datum = NA) +
  labs(fill = "Predicted DOC (μmol kg<sup>-1</sup>)") +
  p_deep_theme


## Assemble
p_deep <- p_deep_a / p_deep_b / p_deep_c + plot_layout(axis_titles = "collect")
p_deep

## Save
ggsave(file = "figures/poster/deep_doc.png", plot = p_deep, width = 120, height = 225, unit = "mm", bg = "white")


## Predictors ----
#--------------------------------------------------------------------------#
p_obs_theme <- theme(
  axis.title = element_blank(),
  text = element_text(size = 8),
  plot.margin = unit(c(0, 0, 0, 0), "mm")
)

# White globe background
bg <- crossing(
  lon = seq(from = -179.5, to = 179.5, by = 1),
  lat = seq(from = -89.5, to = 89.5, by = 1)
)



# Temperature
p_temp <- ggplot(df_ann_epi_pred %>% filter(lon != -0.5)) + # Need to remove lon = -0.5
  geom_tile(data = bg, aes(x = lon, y = lat), fill = "white", colour = "white") +
  geom_sf(data = grat, alpha = 0.8, color = "gray80", linewidth = 0.2) +
  geom_tile(aes(x = lon, y = lat, fill = temperature_surf, colour = temperature_surf), show.legend = FALSE) +
  geom_sf(data = world_rob, fill = "gray80", colour = NA) +
  scale_fill_cmocean(name = "thermal") +
  scale_colour_cmocean(name = "thermal") +
  coord_sf(datum = NA, default_crs = sf::st_crs(4326)) +
  p_obs_theme
ggsave(file = "figures/poster/temp.png", plot = p_temp, width = 120, height = 64, unit = "mm", bg = "transparent")


# Salinity
p_sal <- ggplot(df_ann_epi_pred %>% filter(lon != -0.5)) + # Need to remove lon = -0.5
  geom_tile(data = bg, aes(x = lon, y = lat), fill = "white", colour = "white") +
  geom_sf(data = grat, alpha = 0.8, color = "gray80", linewidth = 0.2) +
  geom_tile(aes(x = lon, y = lat, fill = salinity_surf, colour = salinity_surf), show.legend = FALSE) +
  geom_sf(data = world_rob, fill = "gray80", colour = NA) +
  scale_fill_cmocean(name = "haline") +
  scale_colour_cmocean(name = "haline") +
  coord_sf(datum = NA, default_crs = sf::st_crs(4326)) +
  p_obs_theme
ggsave(file = "figures/poster/sal.png", plot = p_sal, width = 120, height = 64, unit = "mm", bg = "transparent")

# Nitrates
p_nit <- ggplot(df_ann_epi_pred %>% filter(lon != -0.5)) + # Need to remove lon = -0.5
  geom_tile(data = bg, aes(x = lon, y = lat), fill = "white", colour = "white") +
  geom_sf(data = grat, alpha = 0.8, color = "gray80", linewidth = 0.2) +
  geom_tile(aes(x = lon, y = lat, fill = nitrate_surf, colour = nitrate_surf), show.legend = FALSE) +
  geom_sf(data = world_rob, fill = "gray80", colour = NA) +
  scale_fill_cmocean(name = "tempo") +
  scale_colour_cmocean(name = "tempo") +
  coord_sf(datum = NA, default_crs = sf::st_crs(4326)) +
  p_obs_theme
ggsave(file = "figures/poster/nit.png", plot = p_nit, width = 120, height = 64, unit = "mm", bg = "transparent")
