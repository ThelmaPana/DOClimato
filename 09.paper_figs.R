#--------------------------------------------------------------------------#
# Project: DOClimato
# Script purpose: Prepare paper figures
# Date: 27/03/2023
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

## Axis labels
#ylabs <- lapply(c(-90, -60, -30, 0, 30, 60, 90), function(x) {
#  st_sf(label = paste0(abs(x), '\u00b0', 
#                       ifelse(x == 0, '', ifelse(x < 0, 'S', 'N'))),
#        geometry = st_sfc(st_point(c(0, x)), crs = 'WGS84'))
#}) %>% 
#  bind_rows() %>% 
#  st_transform(crs = '+proj=robin +lon_0=180 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs')
#xlabs <- lapply(c(0, 60, 120, 180, 240, 300, 359.9), function(x) {
#  st_sf(label = paste0(round(x), '\u00b0', 'E'),
#        geometry = st_sfc(st_point(c(x, 90)), crs = 'WGS84'))
#}) %>% 
#  bind_rows() %>% 
#  st_transform(crs = '+proj=robin +lon_0=180 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs')


## Prepare all Rsquares and RMSE ----
#--------------------------------------------------------------------------#
# Compute Rsquares and RMSE
rsquares <- res %>% 
  select(resp, season, cv_type, fold, preds) %>% 
  unnest(preds) %>% 
  pivot_longer(log_doc_surf:log_doc_bathy, names_to = "layer", values_to = "truth") %>% 
  drop_na(truth) %>% 
  select(-layer) %>% 
  group_by(resp, season, cv_type, fold) %>% 
  summarise(
    rsq = rsq_vec(truth, .pred),
    rmse = rmse_vec(truth, .pred)
    )


## 1: Surface climatology ----
#--------------------------------------------------------------------------#
# Surface climatology:
# a - annual climatology
# b - uncertainty
# c - amplitude of seasonal cycle

p1_theme <- theme(
  axis.title = element_blank(),
  legend.position = "bottom",
  legend.title = element_markdown(size = 10),
  legend.key.width = unit(1.5, "cm"),
  legend.key.height = unit(2, "mm"),
  legend.direction = "horizontal",
  legend.title.position = "top",
  legend.margin = margin(0, 0, 0, 0),
  legend.box.margin = margin(-10, -10, 0, -10),
  plot.margin = unit(c(0, 0, 0, 0), "mm"),
  text = element_text(size = 8)
  )

## a + b
# Get projections, average and sd by pixel
df_1ab <- res %>% 
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
# Get surface climatology with nice names
surf_clim <- df_1ab %>% rename(surf_doc_avg = doc_avg, surf_doc_sd = doc_sd)


# Plot a
# Colour bar limit
doc_lims <- c(min(df_1ab$doc_avg), 200)
p1a <- ggplot(df_1ab %>% filter(lon != -0.5)) + # Need to remove lon = -0.5
  geom_sf(data = grat, alpha = 0.8, color = "gray80", linewidth = 0.2) +
  geom_tile(aes(x = lon, y = lat, fill = doc_avg, colour = doc_avg)) +
  geom_sf(data = world_rob, fill = "gray80", colour = NA) +
  ggplot2::scale_fill_viridis_c(option = "F", trans = "log1p", limits = doc_lims) +
  ggplot2::scale_colour_viridis_c(option = "F", trans = "log1p", guide = "none", limits = doc_lims) +
  coord_sf(datum = NA, default_crs = sf::st_crs(4326)) +
  labs(fill = "Predicted DOC (μmol kg<sup>-1</sup>)") +
  p1_theme

# Plot b
p1b <- ggplot(df_1ab  %>% filter(lon != -0.5))  +
  geom_sf(data = grat, alpha = 0.8, color = "gray80", linewidth = 0.2) +
  geom_tile(aes(x = lon, y = lat, fill = doc_sd, colour = doc_sd)) +
  geom_sf(data = world_rob, fill = "gray80", colour = NA) +
  ggplot2::scale_fill_viridis_c(option = "E", trans = "log1p", breaks = c(1, 5, 10, 20, 30)) +
  ggplot2::scale_colour_viridis_c(option = "E", trans = "log1p", guide = "none") +
  coord_sf(crs = '+proj=robin +lon_0=180 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs', default_crs = sf::st_crs(4326), datum = NA) +
  labs(fill = "Prediction uncertainty (μmol kg<sup>-1</sup>)") +
  p1_theme 

## c
# Seasonal surface climatology
seas_clim <- res %>% 
  filter(resp == "doc_surf" & season != "0") %>% 
  select(fold, new_preds) %>% 
  unnest(new_preds) %>% 
  mutate(season = as.character(season)) %>% 
  # Apply exp to predictions as we predicted log(doc)
  mutate(pred_doc = exp(pred_doc_log), .after = pred_doc_log) %>% 
  select(season, fold, lon, lat, contains("doc")) %>% 
  left_join(rsquares %>% filter(resp == "doc_surf" & season != "0" & cv_type == "stratified"), by = join_by(season, fold)) %>% 
  # Compute avg prediction per pixel
  group_by(lon, lat, season) %>% 
  summarise(
    doc_avg = wtd.mean(pred_doc, weights = rsq, na.rm = TRUE), 
    doc_sd = sqrt(wtd.var(pred_doc, weights = rsq, na.rm = TRUE)),
    .groups = "drop"
  )

# Compute seasonal amplitude
df_1c <- seas_clim %>% 
  # compute amplitude of seasonal cycle
  group_by(lon, lat) %>% 
  # Keep only pixels where DOC is predicted for all seasons
  filter(n() == 4) %>% 
  summarise(
    seas_amp = max(doc_avg, na.rm = TRUE) - min(doc_avg, na.rm = TRUE), 
    seas_var = var(doc_avg, na.rm = TRUE),
    .groups = "drop"
  )

p1c <- ggplot(df_1c %>% filter(lon != -0.5))  +
  geom_sf(data = grat, alpha = 0.8, color = "gray80", linewidth = 0.2) +
  geom_tile(aes(x = lon, y = lat, fill = seas_amp, colour = seas_amp)) +
  geom_sf(data = world_rob, fill = "gray80", colour = NA) +
  ggplot2::scale_fill_viridis_c(trans = "log1p", breaks = c(1, 5, 10, 20, 30)) +
  ggplot2::scale_colour_viridis_c(trans = "log1p", guide = "none") +
  coord_sf(crs = '+proj=robin +lon_0=180 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs', default_crs = sf::st_crs(4326), datum = NA) +
  labs(fill = "Seasonal amplitude (μmol kg<sup>-1</sup>)") +
  p1_theme 

## Assemble
p1 <- p1a / p1b / p1c + plot_layout(axis_titles = "collect") + plot_annotation(tag_levels = "a")

# Save
ggsave(file = "figures/figure_1.png", plot = p1, width = 110, height = 188, unit = "mm", bg = "white")


## 2: Deeper climatologies ----
#--------------------------------------------------------------------------#
# Deeper climatologies:
# a - epi climatology
# b - meso climatology
# c - bathy climatology

p2_theme <- theme(
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
df_2a <- res %>% 
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
# Get epipelagic climatology with nice names
epi_clim <- df_2a %>% rename(epi_doc_avg = doc_avg, epi_doc_sd = doc_sd)

## b, mesopelagic predictions
df_2b <- res %>% 
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
# Get mesopelagic climatology with nice names
meso_clim <- df_2b %>% rename(meso_doc_avg = doc_avg, meso_doc_sd = doc_sd)

## c, bathypelagic predictions
df_2c <- res %>% 
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
# Get bathypelagic climatology with nice names
bathy_clim <- df_2c %>% rename(bathy_doc_avg = doc_avg, bathy_doc_sd = doc_sd)

# Plot a
p2a <- ggplot(df_2a %>% filter(lon != -0.5))  +
  geom_sf(data = grat, alpha = 0.8, color = "gray80", linewidth = 0.2) +
  geom_tile(aes(x = lon, y = lat, fill = doc_avg, colour = doc_avg)) +
  geom_sf(data = world_rob, fill = "gray80", colour = NA) +
  ggplot2::scale_fill_viridis_c(option = "F", trans = "log1p") +
  ggplot2::scale_colour_viridis_c(option = "F", trans = "log1p", guide = "none") +
  coord_sf(crs = '+proj=robin +lon_0=180 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs', default_crs = sf::st_crs(4326), datum = NA) +
  labs(fill = "Predicted DOC (μmol kg<sup>-1</sup>)") +
  p2_theme

# Plot b
p2b <- ggplot(df_2b %>% filter(lon != -0.5))  +
  geom_sf(data = grat, alpha = 0.8, color = "gray80", linewidth = 0.2) +
  geom_tile(aes(x = lon, y = lat, fill = doc_avg, colour = doc_avg)) +
  geom_sf(data = world_rob, fill = "gray80", colour = NA) +
  ggplot2::scale_fill_viridis_c(option = "F", trans = "log1p") +
  ggplot2::scale_colour_viridis_c(option = "F", trans = "log1p", guide = "none") +
  coord_sf(crs = '+proj=robin +lon_0=180 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs', default_crs = sf::st_crs(4326), datum = NA) +
  labs(fill = "Predicted DOC (μmol kg<sup>-1</sup>)") +
  p2_theme

# Plot c
p2c <- ggplot(df_2c %>% filter(lon != -0.5))  +
  geom_sf(data = grat, alpha = 0.8, color = "gray80", linewidth = 0.2) +
  geom_tile(aes(x = lon, y = lat, fill = doc_avg, colour = doc_avg)) +
  geom_sf(data = world_rob, fill = "gray80", colour = NA) +
  ggplot2::scale_fill_viridis_c(option = "F", trans = "log1p") +
  ggplot2::scale_colour_viridis_c(option = "F", trans = "log1p", guide = "none") +
  coord_sf(crs = '+proj=robin +lon_0=180 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs', default_crs = sf::st_crs(4326), datum = NA) +
  labs(fill = "Predicted DOC (μmol kg<sup>-1</sup>)") +
  p2_theme

## Assemble
p2 <- p2a / p2b / p2c + plot_layout(axis_titles = "collect") + plot_annotation(tag_levels = "a")
p2

## Save
ggsave(file = "figures/figure_2.png", plot = p2, width = 110, height = 188, unit = "mm", bg = "white")


## 3: Variable importance ----
#--------------------------------------------------------------------------#
# Variable importance
# Unnest variable importance

full_vip <- res %>% 
  filter(cv_type == "stratified") %>% 
  select(cv_type, resp, season, fold, importance) %>% 
  unnest(importance) %>%
  mutate(variable = forcats::fct_reorder(variable, dropout_loss)) %>%
  group_by(cv_type, resp, season, fold, variable) %>%
  summarise(dropout_loss = mean(dropout_loss), .groups = "drop")

# Surface
vip_surf_ann <- full_vip %>% 
  filter(variable != "_full_model_") %>%
  filter(resp == "doc_surf" & season == "0") %>% 
  group_by(variable) %>% 
  mutate(mean_dl = mean(dropout_loss)) %>% 
  ungroup() %>% 
  mutate(
    variable = forcats::fct_drop(variable),
    variable = forcats::fct_reorder(variable, mean_dl),
    variable = forcats::fct_other(variable, keep = tail(levels(variable), n = 5), other_level = "other"),
    variable = forcats::fct_reorder(variable, mean_dl)
  )

# Epi
vip_epi_ann <- full_vip %>% 
  filter(variable != "_full_model_") %>%
  filter(resp == "doc_epi" & season == "0") %>% 
  group_by(variable) %>% 
  mutate(mean_dl = mean(dropout_loss)) %>% 
  ungroup() %>% 
  mutate(
    variable = forcats::fct_drop(variable),
    variable = forcats::fct_reorder(variable, mean_dl),
    variable = forcats::fct_other(variable, keep = tail(levels(variable), n = 5), other_level = "other"),
    variable = forcats::fct_reorder(variable, mean_dl)
  ) 

# Meso
vip_meso_ann <- full_vip %>% 
  filter(variable != "_full_model_") %>%
  filter(resp == "doc_meso" & season == "0") %>% 
  group_by(variable) %>% 
  mutate(mean_dl = mean(dropout_loss)) %>% 
  ungroup() %>% 
  mutate(
    variable = forcats::fct_drop(variable),
    variable = forcats::fct_reorder(variable, mean_dl),
    variable = forcats::fct_other(variable, keep = tail(levels(variable), n = 5), other_level = "other"),
    variable = forcats::fct_reorder(variable, mean_dl)
  ) 

# Bathy
vip_bathy_ann <- full_vip %>% 
  filter(variable != "_full_model_") %>%
  filter(resp == "doc_bathy" & season == "0") %>% 
  group_by(variable) %>% 
  mutate(mean_dl = mean(dropout_loss)) %>% 
  ungroup() %>% 
  mutate(
    variable = forcats::fct_drop(variable),
    variable = forcats::fct_reorder(variable, mean_dl),
    variable = forcats::fct_other(variable, keep = tail(levels(variable), n = 5), other_level = "other"),
    variable = forcats::fct_reorder(variable, mean_dl)
  )

## Drop-out loss for full model
full_model_vip <- full_vip %>% 
  filter(variable == "_full_model_") %>% 
  mutate(resp = factor(resp, levels = c("doc_surf", "doc_epi", "doc_meso", "doc_bathy"))) %>% 
  group_by(resp) %>% 
  summarise(dropout_loss = mean(dropout_loss))

# Assemble layers
df_3 <- vip_surf_ann %>% 
  bind_rows(vip_epi_ann) %>% 
  bind_rows(vip_meso_ann) %>% 
  bind_rows(vip_bathy_ann) %>% 
  # Set layer as factor for ordered labels
  mutate(resp = fct_inorder(resp)) %>% 
  # nice formatting for variables
  mutate(
    variable = str_replace_all(variable, "_", " "),
    variable = str_to_sentence(variable),
    variable = str_replace_all(variable, "Mld", "MLD"),
    variable = str_replace_all(variable, "Bbp", "b<sub>bp</sub>"),
    variable = str_replace_all(variable, "Fmicro", "F micro"),
    variable = str_replace_all(variable, "surf", "surf."),
    variable = str_replace_all(variable, "epi", "epi."),
    variable = str_replace_all(variable, "meso", "meso."),
    variable = str_replace_all(variable, "bathy", "bathy."),
    variable = str_replace_all(variable, "Log chl", "log chl a"),
    variable = str_replace_all(variable, "Aou", "AOU")
  ) %>% 
  # Generate interaction between layers and variables to order variables in each layer
  mutate(variable_full = paste0(variable, "-" ,resp)) %>% 
  arrange(mean_dl) %>% 
  mutate(variable_full = fct_inorder(variable_full))

#($\log([\mathrm{chl}_a])$)


# Extract labels for y axis
df_3_names <- df_3 %>% 
  select(variable, variable_full) %>% 
  unique()
labeller <- setNames(str_c(df_3_names$variable), str_c(df_3_names$variable_full))

# Plot
p3 <- ggplot(df_3) + 
  geom_vline(data = full_model_vip, aes(xintercept = dropout_loss), colour = "grey", linewidth = 0.8) +
  geom_boxplot(aes(x = dropout_loss, y = variable_full, colour = resp), outlier.size = 0.3, linewidth = 0.2) +
  scale_colour_manual(
    values = c("#bdd7e7", "#6baed6", "#3182bd", "#08519c"),
    labels = c(
      doc_surf = "Surf.",
      doc_epi = "Epi.",
      doc_meso = "Meso.",
      doc_bathy = "Bathy."
    )) +
  expand_limits(x = 0) +
  labs(x = "RMSE", y = "Variable", colour = "Layer") +
  facet_wrap(~resp, ncol = 1, scales = "free_y") +
  scale_y_discrete(labels = function(x) str_replace_all(x, labeller)) +
  theme(
    strip.background = element_blank(), strip.text.x = element_blank(), # empty facets
    legend.position = "inside", legend.position.inside = c(0.75, 0.4), # legend position
    legend.background = element_rect(fill = "white", colour = NA),
    axis.text.y = element_markdown(),
    text = element_text(size = 8)
  )
p3

## Save
ggsave(file = "figures/figure_3.png", plot = p3, width = 80, height = 100, unit = "mm", bg = "white")


## Format and save climatologies ----
#--------------------------------------------------------------------------#
# Annual climatologies
ann_clim <- surf_clim %>% 
  left_join(epi_clim, by = join_by(lon, lat)) %>% 
  left_join(meso_clim, by = join_by(lon, lat)) %>% 
  left_join(bathy_clim, by = join_by(lon, lat))

# Seasonal surface climatology
seas_clim <- seas_clim %>% rename(surf_doc_avg = doc_avg, surf_doc_sd = doc_sd)

## Save
write_csv(ann_clim, file = "output/annual_climatologies.csv")
write_csv(seas_clim, file = "output/seasonal_climatologies.csv")


## S1: Map of annual observations ----
#--------------------------------------------------------------------------#
load("data/02.ann_surf.Rdata")
load("data/02.ann_epi.Rdata")
load("data/02.ann_meso.Rdata")
load("data/02.ann_bathy.Rdata")
load("data/02.seas_surf.Rdata")

# Count of observations
df_seas_surf_fit %>% count(season)


ps1_theme <- theme(
  axis.title = element_blank(),
  text = element_text(size = 8)
)

# Surface
ps1a <- ggplot(df_ann_surf_fit) +
  geom_sf(data = grat, alpha = 0.8, color = "gray80", linewidth = 0.2) +
  geom_point(aes(x = lon, y = lat), size = 0.05, alpha = 0.3) +
  geom_sf(data = world_rob, fill = "gray80", colour = NA) +
  coord_sf(crs = '+proj=robin +lon_0=180 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs', default_crs = sf::st_crs(4326), datum = NA) +
  ps1_theme

# Epi
ps1b <- ggplot(df_ann_epi_fit) +
  geom_sf(data = grat, alpha = 0.8, color = "gray80", linewidth = 0.2) +
  geom_point(aes(x = lon, y = lat), size = 0.05, alpha = 0.3) +
  geom_sf(data = world_rob, fill = "gray80", colour = NA) +
  coord_sf(crs = '+proj=robin +lon_0=180 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs', default_crs = sf::st_crs(4326), datum = NA) +
  ps1_theme

# Meso
ps1c <- ggplot(df_ann_meso_fit) +
  geom_sf(data = grat, alpha = 0.8, color = "gray80", linewidth = 0.2) +
  geom_point(aes(x = lon, y = lat), size = 0.05, alpha = 0.3) +
  geom_sf(data = world_rob, fill = "gray80", colour = NA) +
  coord_sf(crs = '+proj=robin +lon_0=180 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs', default_crs = sf::st_crs(4326), datum = NA) +
  ps1_theme

# Bathy
ps1d <- ggplot(df_ann_bathy_fit) +
  geom_sf(data = grat, alpha = 0.8, color = "gray80", linewidth = 0.2) +
  geom_point(aes(x = lon, y = lat), size = 0.05, alpha = 0.3) +
  geom_sf(data = world_rob, fill = "gray80", colour = NA) +
  coord_sf(crs = '+proj=robin +lon_0=180 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs', default_crs = sf::st_crs(4326), datum = NA) +
  ps1_theme

ps1 <- ps1a + ps1b + ps1c + ps1d + plot_layout(ncol = 2, guides = "collect") + plot_annotation(tag_levels = "a")
ps1

## Save
ggsave(file = "figures/figure_s1.png", plot = ps1, width = 180, height = 110, unit = "mm", bg = "white")


## S2: Map of seasonal observations ----
#--------------------------------------------------------------------------#
load("data/02.seas_surf.Rdata")

# Count of observations
df_seas_surf_fit %>% count(season)

ps2_theme <- theme(
  axis.title = element_blank(),
  text = element_text(size = 8)
)

# Winter
ps2a <- df_seas_surf_fit %>% 
  filter(season == 1) %>% 
  ggplot() +
  geom_sf(data = grat, alpha = 0.8, color = "gray80", linewidth = 0.2) +
  geom_point(aes(x = lon, y = lat), size = 0.05, alpha = 0.3) +
  geom_sf(data = world_rob, fill = "gray80", colour = NA) +
  coord_sf(crs = '+proj=robin +lon_0=180 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs', default_crs = sf::st_crs(4326), datum = NA) +
  ps2_theme

# Spring
ps2b <- df_seas_surf_fit %>% 
  filter(season == 2) %>% 
  ggplot() +
  geom_sf(data = grat, alpha = 0.8, color = "gray80", linewidth = 0.2) +
  geom_point(aes(x = lon, y = lat), size = 0.05, alpha = 0.3) +
  geom_sf(data = world_rob, fill = "gray80", colour = NA) +
  coord_sf(crs = '+proj=robin +lon_0=180 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs', default_crs = sf::st_crs(4326), datum = NA) +
  ps1_theme

# Summer
ps2c <- df_seas_surf_fit %>% 
  filter(season == 3) %>% 
  ggplot() +
  geom_sf(data = grat, alpha = 0.8, color = "gray80", linewidth = 0.2) +
  geom_point(aes(x = lon, y = lat), size = 0.05, alpha = 0.3) +
  geom_sf(data = world_rob, fill = "gray80", colour = NA) +
  coord_sf(crs = '+proj=robin +lon_0=180 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs', default_crs = sf::st_crs(4326), datum = NA) +
  ps1_theme

# Autumn
ps2d <- df_seas_surf_fit %>% 
  filter(season == 4) %>% 
  ggplot() +
  geom_sf(data = grat, alpha = 0.8, color = "gray80", linewidth = 0.2) +
  geom_point(aes(x = lon, y = lat), size = 0.05, alpha = 0.3) +
  geom_sf(data = world_rob, fill = "gray80", colour = NA) +
  coord_sf(crs = '+proj=robin +lon_0=180 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs', default_crs = sf::st_crs(4326), datum = NA) +
  ps1_theme

ps2 <- ps2a + ps2b + ps2c + ps2d + plot_layout(ncol = 2, guides = "collect") + plot_annotation(tag_levels = "a")
ps2

## Save
ggsave(file = "figures/figure_s2.png", plot = ps2, width = 180, height = 110, unit = "mm", bg = "white")


 ## S3: Distribution of predictors for annual predictions ----
#--------------------------------------------------------------------------#
# Distribution of predictors in the learning set VS in new data for annual prediction

ps3_theme <- theme(
  axis.text.x = element_text(angle = 45, hjust = 1),
  axis.title.x = element_blank(),
  strip.background = element_blank(),
  strip.placement = "outside",
  legend.margin = margin(0, 0, 0, 0),
  legend.box.margin = margin(-10, -10, 0, -10),
  plot.margin = unit(c(0, 0, 0, 0), "mm"),
  text = element_text(size = 8),
  strip.text.x = element_text(vjust = 5), 
  panel.spacing = unit(0, "lines")
)

ps3_facet <- facet_wrap(~variable, scales = "free", ncol = 5, strip.position = "bottom")

# Surface
pred_surf <- bind_rows(
  df_ann_surf_fit %>% select(temperature_surf:par) %>% pivot_longer(everything(), names_to = "variable") %>% mutate(type = "Annotated"),
  df_ann_surf_pred %>% select(temperature_surf:par) %>% pivot_longer(everything(), names_to = "variable") %>% mutate(type = "New")
)
  
ps3a <- ggplot(pred_surf) + 
  geom_density(aes(x = value, linetype = type), linewidth = 0.3) +
  labs(x = "Value", y = "Density", linetype = "Data type") +
  ps3_facet +
  ps3_theme

# Epi
pred_epi <- bind_rows(
  df_ann_epi_fit %>% select(temperature_surf:nitrate_epi) %>% pivot_longer(everything(), names_to = "variable") %>% mutate(type = "Annotated"),
  df_ann_epi_pred %>% select(temperature_surf:nitrate_epi) %>% pivot_longer(everything(), names_to = "variable") %>% mutate(type = "New")
) %>% 
  filter(!variable %in% pred_surf$variable) # drop predictors already shown for surface layer

ps3b <- ggplot(pred_epi) + 
  geom_density(aes(x = value, linetype = type), linewidth = 0.3) +
  labs(x = "Value", y = "Density", linetype = "Data type") +
  ps3_facet +
  ps3_theme

# Meso
pred_meso <- bind_rows(
  df_ann_meso_fit %>% select(temperature_surf:nitrate_meso) %>% pivot_longer(everything(), names_to = "variable") %>% mutate(type = "Annotated"),
  df_ann_meso_pred %>% select(temperature_surf:nitrate_meso) %>% pivot_longer(everything(), names_to = "variable") %>% mutate(type = "New")
) %>% 
  filter(!variable %in% pred_surf$variable) %>% # drop predictors already shown for surface layer
  filter(!variable %in% pred_epi$variable) # and those already shown for epi layer

ps3c <- ggplot(pred_meso) + 
  geom_density(aes(x = value, linetype = type), linewidth = 0.3) +
  labs(x = "Value", y = "Density", linetype = "Data type") +
  ps3_facet +
  ps3_theme

# Bathy
pred_bathy <- bind_rows(
  df_ann_bathy_fit %>% select(temperature_surf:nitrate_bathy) %>% pivot_longer(everything(), names_to = "variable") %>% mutate(type = "Annotated"),
  df_ann_bathy_pred %>% select(temperature_surf:nitrate_bathy) %>% pivot_longer(everything(), names_to = "variable") %>% mutate(type = "New")
) %>% 
  filter(!variable %in% pred_surf$variable) %>% # drop predictors already shown for surface layer
  filter(!variable %in% pred_epi$variable) %>% # and those already shown for epi layer
  filter(!variable %in% pred_meso$variable) # and those already shown for meso layer

ps3d <- ggplot(pred_bathy) + 
  geom_density(aes(x = value, linetype = type), linewidth = 0.3) +
  labs(x = "Value", y = "Density", linetype = "Data type") +
  ps3_facet +
  ps3_theme

# Assemble
ps3 <- ps3a / ps3b / ps3c / ps3d + plot_annotation(tag_levels = "a") + plot_layout(guides = "collect", heights = c(6, 2, 2, 2))
ps3


# Save
ggsave(file = "figures/figure_s3.png", plot = ps3, width = 188, height = 250, unit = "mm", bg = "white")



## S4: Distribution of predictors for annual predictions ----
#--------------------------------------------------------------------------#
# Distribution of predictors in the learning set VS in new data for seasonal prediction

ps4_facet <- ps3_facet # use same facet as previous figure
ps4_theme <- theme(
  axis.text.x = element_text(angle = 45, hjust = 1),
  axis.title.x = element_blank(),
  strip.background = element_blank(),
  strip.placement = "outside",
  text = element_text(size = 8),
  strip.text.x = element_text(vjust = 5), 
  panel.spacing = unit(0, "lines")
)

df_s4 <- bind_rows(
  df_seas_surf_fit %>% select(season, temperature_surf:par) %>% pivot_longer(temperature_surf:par, names_to = "variable") %>% mutate(type = "Annotated"),
  df_seas_surf_pred %>% select(season, temperature_surf:par) %>% pivot_longer(temperature_surf:par, names_to = "variable") %>% mutate(type = "New")
) %>% 
  mutate(season = as.character(season))

# Here we need to save figure by figure because they are big. 
# Some formatting for each season
seasons <- tibble(
  season = c(1:4) %>% as.character(),
  colour = c("#F29762", "#3F9D86", "#486E9E", "#F2BB62"),
  tag = c("a", "b", "c", "d")
)

for (i in 1:nrow(seasons)) {
  # Generate plot
  p <- df_s4 %>% 
    filter(season == seasons$season[i]) %>% 
    ggplot() + 
    geom_density(aes(x = value, linetype = type), colour = seasons$colour[i], linewidth = 0.3) +
    labs(x = "Value", y = "Density", linetype = "Data type", tag = seasons$tag[i]) +
    ps4_facet +
    ps4_theme
  
  # Save it
  ggsave(file = paste0("figures/figure_s4", seasons$tag[i], ".png"), plot = p, width = 188, height = 120, unit = "mm", bg = "white")
}


## S5: Rsquares and RMSE dist in all layers ----
#--------------------------------------------------------------------------#
ps5_theme <- theme(
  text = element_text(size = 8),
  axis.text.x = element_text(angle = 45, hjust = 1)
)

# Rsquares annual
ps5a <- rsquares %>% 
  left_join(resp_to_layer, by = join_by(resp)) %>% 
  filter(season == "0") %>% 
  mutate(bp_factor = interaction(layer, cv_type)) %>% 
  ggplot() +
  geom_boxplot(aes(x = layer, y = rsq, group = bp_factor, colour = cv_type), position = position_dodge(1), linewidth = 0.2, outlier.size = 0.2) +
  scale_colour_manual(values = c("grey", "black")) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  labs(x = "Layer", y = "R²", colour = "CV type") +
  ps5_theme

# Rsquares seasonnal
ps5b <- rsquares %>% 
  left_join(resp_to_layer, by = join_by(resp)) %>% 
  filter(season != "0") %>% 
  mutate(bp_factor = interaction(layer, season)) %>% 
  ggplot() +
  geom_boxplot(aes(x = layer, y = rsq, group = bp_factor, colour = season), position = position_dodge(1),  linewidth = 0.2, outlier.size = 0.2) +
  scale_colour_manual(
    values = c( "#F29762", "#3F9D86", "#486E9E", "#F2BB62"),
    labels = c(
      `1` = "winter",
      `2` = "spring",
      `3` = "summer",
      `4` = "autumn"
    )
    ) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  labs(x = "Layer", y = "R²", colour = "Season") +
  ps5_theme

# RMSE annual
ps5c <- rsquares %>% 
  left_join(resp_to_layer, by = join_by(resp)) %>% 
  filter(season == "0") %>% 
  mutate(bp_factor = interaction(layer, cv_type)) %>% 
  ggplot() +
  geom_boxplot(aes(x = layer, y = rmse, group = bp_factor, colour = cv_type), position = position_dodge(1), linewidth = 0.2, outlier.size = 0.2) +
  scale_colour_manual(values = c("grey", "black")) +
  scale_y_continuous(limits = c(0, 0.2), expand = c(0, 0)) +
  labs(x = "Layer", y = "RMSE", colour = "CV type") +
  ps5_theme

# RMSE seasonnal
ps5d <- rsquares %>% 
  left_join(resp_to_layer, by = join_by(resp)) %>% 
  filter(season != "0") %>% 
  mutate(bp_factor = interaction(layer, season)) %>% 
  ggplot() +
  geom_boxplot(aes(x = layer, y = rmse, group = bp_factor, colour = season), position = position_dodge(1),  linewidth = 0.2, outlier.size = 0.2) +
  scale_colour_manual(
    values = c("#F29762", "#3F9D86", "#486E9E", "#F2BB62"),
    labels = c(
      `1` = "winter",
      `2` = "spring",
      `3` = "summer",
      `4` = "autumn"
    )
  ) +
  scale_y_continuous(limits = c(0, 0.2), expand = c(0, 0)) +
  labs(x = "Layer", y = "RMSE", colour = "Season") +
  ps5_theme

# Assemble
ps5 <- (
  ((ps5a / ps5c) + plot_layout(guides = "collect", axis_titles = "collect")) | # Left column
  ((ps5b / ps5d) + plot_layout(guides = "collect", axis_titles = "collect")) # Right column
  ) + 
  plot_annotation(tag_levels = list(c("a", "c", "b", "d"))) # Reordered labels

# Save
ggsave(file = "figures/figure_s5.png", plot = ps5, width = 150, height = 80, unit = "mm", bg = "white")


# Nice table of rsquares and RMSE
rsquares %>% 
  filter(season == "0") %>% 
  filter(cv_type == "stratified") %>% 
  #filter(cv_type == "spatial") %>% 
  group_by(resp) %>% 
  summarise(
    mean = mean(rsq),
    sd = sd(rsq)
  )

rsquares %>% 
  filter(season == "0") %>% 
  #filter(cv_type == "stratified") %>% 
  filter(cv_type == "spatial") %>% 
  group_by(resp) %>% 
  summarise(
    mean = mean(rmse),
    sd = sd(rmse)
  )

rsquares %>% 
  filter(season != "0") %>% 
  filter(cv_type == "stratified") %>% 
  #filter(cv_type == "spatial") %>% 
  group_by(season) %>% 
  summarise(
    mean = mean(rmse),
    sd = sd(rmse)
  )

## S6: Model error VS measure error ----
#--------------------------------------------------------------------------#
# Compare RMSE to within pixel sd values
load("data/01.doc_log_data.Rdata")

# RMSE for non log-transformed values
rmse_log <- res %>% 
  select(resp, season, cv_type, fold, preds) %>% 
  # Unnest predictions
  unnest(preds) %>% 
  # Reformat
  pivot_longer(log_doc_surf:log_doc_bathy, names_to = "layer", values_to = "truth") %>% 
  drop_na(truth) %>% 
  select(-layer) %>% 
  group_by(resp, season, cv_type, fold) %>% 
  # Compute RMSE for each layer / season / cv_type and fold
  summarise(
    rmse = rmse_vec(truth, .pred), .groups = "drop"
  ) %>% 
  # Join with nice layer names
  left_join(resp_to_layer, by = join_by(resp)) %>% 
  # Some formatting to join with sd values
  rename(value = rmse) %>% 
  mutate(metric = "rmse") %>% 
  select(layer, season, cv_type, metric, value)

sd_log <- df_doc_log %>% 
  # Convert season to character
  mutate(season = as.character(season)) %>% 
  # Keep only sd columns
  select(season, contains("sd")) %>% 
  # Reformat
  pivot_longer(surf_sd:bathy_sd, names_to = "resp", values_to = "sd") %>% 
  # Nice layer names
  mutate(resp = paste0("doc_", str_remove(resp, "_sd"))) %>% 
  left_join(resp_to_layer, by = join_by(resp)) %>% 
  # Some formatting to join with RMSE values
  rename(value = sd) %>% 
  mutate(metric = "sd") %>% 
  select(layer, season, metric, value)

df_s6 <- bind_rows(rmse_log, sd_log) %>% 
  # Regroup both CV type for RMSE and SD in one metric name
  mutate(metric = ifelse(metric == "rmse", paste0("RMSE (", cv_type, " CV)"), "Within pixel SD\nof log-transformed\nDOC values"))

# Annual
ps6a <- df_s6 %>% 
  filter(season == "0") %>% 
  ggplot() +
  geom_boxplot(aes(x = layer, y = value, colour = metric, linetype = metric), linewidth = 0.2, outlier.size = 0.2) +
  scale_y_log10() +
  scale_linetype_manual(values = c("solid", "solid", "22")) +
  scale_colour_manual(values = c("grey", "black", "#756bb1")) +
  labs(x = "Layer", y = "Metric value", colour = "Metric", linetype = "Metric") +
  ps5_theme

# Seasonal
ps6b <- df_s6 %>% 
  filter(season != "0") %>% 
  filter(layer == "Surf.") %>% 
  ggplot() +
  geom_boxplot(aes(x = layer, y = value, colour = season, linetype = metric), linewidth = 0.2, outlier.size = 0.2) +
  scale_colour_manual(
    values = c( "#F29762", "#3F9D86", "#486E9E", "#F2BB62"),
    labels = c(
      `1` = "winter",
      `2` = "spring",
      `3` = "summer",
      `4` = "autumn"
    )
  ) +
  scale_linetype_manual(values = linetype_pal()(2)) +
  scale_y_log10() +
  labs(x = "Layer", y = "Metric value", linetype = "Metric", colour = "Season") +
  ps5_theme +
  theme(legend.margin = margin(t = 0, b = 0, unit = "mm"))


ps6 <- ps6a + ps6b + plot_annotation(tag_levels = "a")
ps6

# Save
ggsave(file = "figures/figure_s6.png", plot = ps6, width = 180, height = 60, unit = "mm", bg = "white")


## S7: Spatial VS stratified in surface layer ----
#--------------------------------------------------------------------------#
# - rsquares value
# - diff in projections

ps7_theme <- theme(
  axis.title = element_blank(),
  legend.position = "bottom",
  legend.title = element_markdown(size = 10),
  legend.key.width = unit(2.5, "cm"),
  legend.key.height = unit(2, "mm"),
  legend.direction = "horizontal",
  legend.title.position = "top",
  text = element_text(size = 8)
)

# Get projections, average and sd by pixel for stratified
df_7sa <- res %>% 
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
    .groups = "drop"
  )

# Get projections, average and sd by pixel for spatial
df_7sb <- res %>% 
  filter(resp == "doc_surf" & season == "0") %>% 
  filter(cv_type == "spatial") %>% 
  select(fold, new_preds) %>% 
  unnest(new_preds) %>% 
  # Apply exp to predictions as we predicted log(doc)
  mutate(pred_doc = exp(pred_doc_log), .after = pred_doc_log) %>% 
  select(fold, lon, lat, contains("doc")) %>% 
  left_join(rsquares %>% filter(resp == "doc_surf" & season == "0" & cv_type == "spatial"), by = join_by(fold)) %>% 
  group_by(lon, lat) %>% 
  summarise(
    doc_avg = wtd.mean(pred_doc, weights = rsq, na.rm = TRUE), 
    .groups = "drop"
  )

# Compute difference between stratified and spatial
df_7s <- df_7sa %>% 
  rename(strat = doc_avg) %>% 
  left_join(df_7sb %>% rename(spat = doc_avg), by = join_by(lon, lat)) %>% 
  mutate(diff = strat - spat)

# Plot it
ps7 <- ggplot(df_7s %>% filter(lon != -0.5)) +
  geom_sf(data = grat, alpha = 0.8, color = "gray80", linewidth = 0.2) +
  geom_tile(aes(x = lon, y = lat, fill = diff, colour = diff)) + 
  geom_sf(data = world_rob, fill = "gray80", colour = NA) +
  scale_fill_gradient2(low = "#4575b4", mid = "#ffffbf", high = "#d73027") +
  scale_colour_gradient2(low = "#4575b4", mid = "#ffffbf", high = "#d73027", guide = "none") +
  coord_sf(crs = '+proj=robin +lon_0=180 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs', default_crs = sf::st_crs(4326), datum = NA) +
  labs(fill = "Difference between projection from stratified CV and spatial CV (μmol kg<sup>-1</sup>)") +
  ps7_theme

# Save
ggsave(file = "figures/figure_s7.png", plot = ps7, width = 180, height = 90, unit = "mm", bg = "white")


## S8: Seasonal projections in the surface layer ----
#--------------------------------------------------------------------------#
ps8_theme <- theme(
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

df_s8 <- res %>% 
  filter(resp == "doc_surf" & season != "0") %>% 
  select(fold, new_preds) %>% 
  unnest(new_preds) %>% 
  mutate(season = as.character(season)) %>% 
  # Apply exp to predictions as we predicted log(doc)
  mutate(pred_doc = exp(pred_doc_log), .after = pred_doc_log) %>% 
  select(season, fold, lon, lat, contains("doc")) %>% 
  left_join(rsquares %>% filter(resp == "doc_surf" & season != "0" & cv_type == "stratified"), by = join_by(season, fold)) %>% 
  # Compute avg prediction per pixel
  group_by(lon, lat, season) %>% 
  summarise(
    doc_avg = wtd.mean(pred_doc, weights = rsq, na.rm = TRUE),
    doc_sd = sqrt(wtd.var(pred_doc, weights = rsq, na.rm = TRUE)), 
    .groups = "drop"
  ) 

doc_avg_lims <- c(min(df_s8$doc_avg), max(df_s8$doc_avg))

ps8a <- ggplot(df_s8 %>% filter(season == 1) %>% filter(lon != -0.5)) + # Need to remove lon = -0.5
  geom_sf(data = grat, alpha = 0.8, color = "gray80", linewidth = 0.2) +
  geom_tile(aes(x = lon, y = lat, fill = doc_avg, colour = doc_avg)) +
  geom_sf(data = world_rob, fill = "gray80", colour = NA) +
  ggplot2::scale_fill_viridis_c(option = "F", trans = "log1p", limits = doc_avg_lims) +
  ggplot2::scale_colour_viridis_c(option = "F", trans = "log1p", guide = "none", limits = doc_avg_lims) +
  coord_sf(crs = '+proj=robin +lon_0=180 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs', default_crs = sf::st_crs(4326), datum = NA) +
  labs(fill = "Predicted DOC (μmol kg<sup>-1</sup>)") +
  ps8_theme

ps8b <- ggplot(df_s8 %>% filter(season == 2) %>% filter(lon != -0.5)) + # Need to remove lon = -0.5
  geom_sf(data = grat, alpha = 0.8, color = "gray80", linewidth = 0.2) +
  geom_tile(aes(x = lon, y = lat, fill = doc_avg, colour = doc_avg)) +
  geom_sf(data = world_rob, fill = "gray80", colour = NA) +
  ggplot2::scale_fill_viridis_c(option = "F", trans = "log1p", limits = doc_avg_lims) +
  ggplot2::scale_colour_viridis_c(option = "F", trans = "log1p", guide = "none", limits = doc_avg_lims) +
  coord_sf(crs = '+proj=robin +lon_0=180 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs', default_crs = sf::st_crs(4326), datum = NA) +
  labs(fill = "Predicted DOC (μmol kg<sup>-1</sup>)") +
  ps8_theme

ps8c <- ggplot(df_s8 %>% filter(season == 3) %>% filter(lon != -0.5)) + # Need to remove lon = -0.5
  geom_sf(data = grat, alpha = 0.8, color = "gray80", linewidth = 0.2) +
  geom_tile(aes(x = lon, y = lat, fill = doc_avg, colour = doc_avg)) +
  geom_sf(data = world_rob, fill = "gray80", colour = NA) +
  ggplot2::scale_fill_viridis_c(option = "F", trans = "log1p", limits = doc_avg_lims) +
  ggplot2::scale_colour_viridis_c(option = "F", trans = "log1p", guide = "none", limits = doc_avg_lims) +
  coord_sf(crs = '+proj=robin +lon_0=180 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs', default_crs = sf::st_crs(4326), datum = NA) +
  labs(fill = "Predicted DOC (μmol kg<sup>-1</sup>)") +
  ps8_theme

ps8d <- ggplot(df_s8 %>% filter(season == 4) %>% filter(lon != -0.5)) + # Need to remove lon = -0.5
  geom_sf(data = grat, alpha = 0.8, color = "gray80", linewidth = 0.2) +
  geom_tile(aes(x = lon, y = lat, fill = doc_avg, colour = doc_avg)) +
  geom_sf(data = world_rob, fill = "gray80", colour = NA) +
  ggplot2::scale_fill_viridis_c(option = "F", trans = "log1p", limits = doc_avg_lims) +
  ggplot2::scale_colour_viridis_c(option = "F", trans = "log1p", guide = "none", limits = doc_avg_lims) +
  coord_sf(crs = '+proj=robin +lon_0=180 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs', default_crs = sf::st_crs(4326), datum = NA) +
  labs(fill = "Predicted DOC (μmol kg<sup>-1</sup>)") +
  ps8_theme

ps8 <- ps8a + ps8b + ps8c + ps8d + plot_annotation(tag_levels = "a") + plot_layout(guides = "collect") & theme(legend.position = 'bottom')
ps8

## Save
ggsave(file = "figures/figure_s8.png", plot = ps8, width = 180, height = 120, unit = "mm", bg = "white")


## S9: Uncertainty for seasonal surface climatologies ----
#--------------------------------------------------------------------------#
# Use same df as previous figure
doc_sd_lims <- c(min(df_s8$doc_sd), max(df_s8$doc_sd))
ps9_theme <- ps8_theme

ps9a <- ggplot(df_s8 %>% filter(season == 1) %>% filter(lon != -0.5)) + # Need to remove lon = -0.5
  geom_sf(data = grat, alpha = 0.8, color = "gray80", linewidth = 0.2) +
  geom_tile(aes(x = lon, y = lat, fill = doc_sd, colour = doc_sd)) +
  geom_sf(data = world_rob, fill = "gray80", colour = NA) +
  ggplot2::scale_fill_viridis_c(option = "E", trans = "log1p", limits = doc_sd_lims) +
  ggplot2::scale_colour_viridis_c(option = "E", trans = "log1p", guide = "none", limits = doc_sd_lims) +
  coord_sf(crs = '+proj=robin +lon_0=180 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs', default_crs = sf::st_crs(4326), datum = NA) +
  labs(fill = "Prediction uncertainty (μmol kg<sup>-1</sup>)") +
  ps9_theme

ps9b <- ggplot(df_s8 %>% filter(season == 2) %>% filter(lon != -0.5)) + # Need to remove lon = -0.5
  geom_sf(data = grat, alpha = 0.8, color = "gray80", linewidth = 0.2) +
  geom_tile(aes(x = lon, y = lat, fill = doc_sd, colour = doc_sd)) +
  geom_sf(data = world_rob, fill = "gray80", colour = NA) +
  ggplot2::scale_fill_viridis_c(option = "E", trans = "log1p", limits = doc_sd_lims) +
  ggplot2::scale_colour_viridis_c(option = "E", trans = "log1p", guide = "none", limits = doc_sd_lims) +
  coord_sf(crs = '+proj=robin +lon_0=180 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs', default_crs = sf::st_crs(4326), datum = NA) +
  labs(fill = "Prediction uncertainty (μmol kg<sup>-1</sup>)") +
  ps9_theme

ps9c <- ggplot(df_s8 %>% filter(season == 3) %>% filter(lon != -0.5)) + # Need to remove lon = -0.5
  geom_sf(data = grat, alpha = 0.8, color = "gray80", linewidth = 0.2) +
  geom_tile(aes(x = lon, y = lat, fill = doc_sd, colour = doc_sd)) +
  geom_sf(data = world_rob, fill = "gray80", colour = NA) +
  ggplot2::scale_fill_viridis_c(option = "E", trans = "log1p", limits = doc_sd_lims) +
  ggplot2::scale_colour_viridis_c(option = "E", trans = "log1p", guide = "none", limits = doc_sd_lims) +
  coord_sf(crs = '+proj=robin +lon_0=180 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs', default_crs = sf::st_crs(4326), datum = NA) +
  labs(fill = "Prediction uncertainty (μmol kg<sup>-1</sup>)") +
  ps9_theme

ps9d <- ggplot(df_s8 %>% filter(season == 4) %>% filter(lon != -0.5)) + # Need to remove lon = -0.5
  geom_sf(data = grat, alpha = 0.8, color = "gray80", linewidth = 0.2) +
  geom_tile(aes(x = lon, y = lat, fill = doc_sd, colour = doc_sd)) +
  geom_sf(data = world_rob, fill = "gray80", colour = NA) +
  ggplot2::scale_fill_viridis_c(option = "E", trans = "log1p", limits = doc_sd_lims) +
  ggplot2::scale_colour_viridis_c(option = "E", trans = "log1p", guide = "none", limits = doc_sd_lims) +
  coord_sf(crs = '+proj=robin +lon_0=180 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs', default_crs = sf::st_crs(4326), datum = NA) +
  labs(fill = "Prediction uncertainty (μmol kg<sup>-1</sup>)") +
  ps9_theme

ps9 <- ps9a + ps9b + ps9c + ps9d + plot_annotation(tag_levels = "a") + plot_layout(guides = "collect") & theme(legend.position = 'bottom')
ps9
## Save
ggsave(file = "figures/figure_s9.png", plot = ps9, width = 180, height = 120, unit = "mm", bg = "white")


## S10: Uncertainty for deeper climatologies ----
#--------------------------------------------------------------------------#
# Uncertainties for deeper climatologies:
# a - epi uncertainty
# b - meso uncertainty
# c - bathy uncertainty

ps10_theme <- p2_theme

# Get common colourbar limits for meso and bathy
doc_sd_lims <- c(
  min(df_2a$doc_sd, df_2b$doc_sd, df_2c$doc_sd), 
  max(df_2a$doc_sd, df_2b$doc_sd, df_2c$doc_sd)  
)

# Plot a
ps10a <- ggplot(df_2a %>% filter(lon != -0.5)) + # Need to remove lon = -0.5
  geom_sf(data = grat, alpha = 0.8, color = "gray80", linewidth = 0.2) +
  geom_tile(aes(x = lon, y = lat, fill = doc_sd, colour = doc_sd)) +
  geom_sf(data = world_rob, fill = "gray80", colour = NA) +
  ggplot2::scale_fill_viridis_c(option = "E", trans = "log1p", limits = doc_sd_lims, breaks = c(1, 5, 10)) +
  ggplot2::scale_colour_viridis_c(option = "E", trans = "log1p", limits = doc_sd_lims, guide = "none") +
  coord_sf(crs = '+proj=robin +lon_0=180 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs', default_crs = sf::st_crs(4326), datum = NA) +
  labs(fill = "Prediction uncertainty (μmol kg<sup>-1</sup>)") +
  ps10_theme

# Plot b
ps10b <- ggplot(df_2b %>% filter(lon != -0.5)) + # Need to remove lon = -0.5
  geom_sf(data = grat, alpha = 0.8, color = "gray80", linewidth = 0.2) +
  geom_tile(aes(x = lon, y = lat, fill = doc_sd, colour = doc_sd)) +
  geom_sf(data = world_rob, fill = "gray80", colour = NA) +
  ggplot2::scale_fill_viridis_c(option = "E", trans = "log1p", limits = doc_sd_lims, breaks = c(1, 5, 10)) +
  ggplot2::scale_colour_viridis_c(option = "E", trans = "log1p", limits = doc_sd_lims, guide = "none") +
  coord_sf(crs = '+proj=robin +lon_0=180 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs', default_crs = sf::st_crs(4326), datum = NA) +
  labs(fill = "Prediction uncertainty (μmol kg<sup>-1</sup>)") +
  ps10_theme

# Plot c
ps10c <- ggplot(df_2c %>% filter(lon != -0.5)) + # Need to remove lon = -0.5
  geom_sf(data = grat, alpha = 0.8, color = "gray80", linewidth = 0.2) +
  geom_tile(aes(x = lon, y = lat, fill = doc_sd, colour = doc_sd)) +
  geom_sf(data = world_rob, fill = "gray80", colour = NA) +
  ggplot2::scale_fill_viridis_c(option = "E", trans = "log1p", limits = doc_sd_lims, breaks = c(1, 5, 10)) +
  ggplot2::scale_colour_viridis_c(option = "E", trans = "log1p", limits = doc_sd_lims, guide = "none") +
  coord_sf(crs = '+proj=robin +lon_0=180 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs', default_crs = sf::st_crs(4326), datum = NA) +
  labs(fill = "Prediction uncertainty (μmol kg<sup>-1</sup>)") +
  ps10_theme

## Assemble
ps10 <- ps10a / ps10b / ps10c + plot_layout(axis_titles = "collect", guides = "collect") + plot_annotation(tag_levels = "a") & theme(legend.position = 'bottom')
ps10

## Save
ggsave(file = "figures/figure_s10.png", plot = ps10, width = 110, height = 188, unit = "mm", bg = "white")


## S11: Partial dependance plots ----
#--------------------------------------------------------------------------#
n_pdp <- 3 # number of pdp for each layer


ps11_theme <- theme(
  axis.title.x = element_blank(), 
  strip.placement = "outside",
  text = element_text(size = 8)
)
ps11_facet <- facet_wrap(~var_name, scales = "free_x", strip.position = "bottom")


## Surface
# Variables to use
vars_surf <- vip_surf_ann %>% 
  filter(variable != "other") %>% 
  group_by(cv_type, variable) %>% 
  summarise(dropout_loss = mean(dropout_loss), .groups = "drop") %>% 
  arrange(desc(dropout_loss)) %>% 
  slice(1:n_pdp)

# CP profiles
cp_surf <- res %>% 
  filter(cv_type == "stratified") %>% 
  filter(resp == "doc_surf" & season == "0") %>% 
  select(cv_type, fold, cp_profiles) %>% 
  unnest(cp_profiles)

# Average across folds
pdp_surf <- prep_pdp(cp = cp_surf, vars = vars_surf)

ps11a <- pdp_surf %>% 
  left_join(vars_surf %>% rename(var_name = variable), by = join_by(var_name)) %>% 
  mutate(var_name = fct_reorder(var_name, dropout_loss, .desc = TRUE)) %>% 
  ggplot() +
  geom_path(aes(x = x, y = yhat_loc)) +
  geom_ribbon(aes(x = x, ymin = yhat_loc - yhat_spr, ymax = yhat_loc + yhat_spr), alpha = 0.2) +
  labs(y = "Predicted log(DOC)") +
  ps11_facet +
  ps11_theme


## Epi
# Variables to use
vars_epi <- vip_epi_ann %>% 
  filter(variable != "other") %>% 
  group_by(cv_type, variable) %>% 
  summarise(dropout_loss = mean(dropout_loss), .groups = "drop") %>% 
  arrange(desc(dropout_loss)) %>% 
  slice(1:n_pdp)

# CP profiles
cp_epi <- res %>% 
  filter(cv_type == "stratified") %>% 
  filter(resp == "doc_epi" & season == "0") %>% 
  select(cv_type, fold, cp_profiles) %>% 
  unnest(cp_profiles)

# Average across folds
pdp_epi <- prep_pdp(cp = cp_epi, vars = vars_epi)

ps11b <- pdp_epi %>% 
  left_join(vars_epi %>% rename(var_name = variable), by = join_by(var_name)) %>% 
  mutate(var_name = fct_reorder(var_name, dropout_loss, .desc = TRUE)) %>% 
  ggplot() +
  geom_path(aes(x = x, y = yhat_loc)) +
  geom_ribbon(aes(x = x, ymin = yhat_loc - yhat_spr, ymax = yhat_loc + yhat_spr), alpha = 0.2) +
  labs(y = "Predicted log(DOC)") +
  ps11_facet +
  ps11_theme


## Meso
# Variables to use
vars_meso <- vip_meso_ann %>% 
  filter(variable != "other") %>% 
  group_by(cv_type, variable) %>% 
  summarise(dropout_loss = mean(dropout_loss), .groups = "drop") %>% 
  arrange(desc(dropout_loss)) %>% 
  slice(1:n_pdp)

# CP profiles
cp_meso <- res %>% 
  filter(cv_type == "stratified") %>% 
  filter(resp == "doc_meso" & season == "0") %>% 
  select(cv_type, fold, cp_profiles) %>% 
  unnest(cp_profiles)

# Average across folds
pdp_meso <- prep_pdp(cp = cp_meso, vars = vars_meso)

ps11c <- pdp_meso %>% 
  left_join(vars_meso %>% rename(var_name = variable), by = join_by(var_name)) %>% 
  mutate(var_name = fct_reorder(var_name, dropout_loss, .desc = TRUE)) %>% 
  ggplot() +
  geom_path(aes(x = x, y = yhat_loc)) +
  geom_ribbon(aes(x = x, ymin = yhat_loc - yhat_spr, ymax = yhat_loc + yhat_spr), alpha = 0.2) +
  labs(y = "Predicted log(DOC)") +
  ps11_facet +
  ps11_theme


## Bathy
# Variables to use
vars_bathy <- vip_bathy_ann %>% 
  filter(variable != "other") %>% 
  group_by(cv_type, variable) %>% 
  summarise(dropout_loss = mean(dropout_loss), .groups = "drop") %>% 
  arrange(desc(dropout_loss)) %>% 
  slice(1:n_pdp)

# CP profiles
cp_bathy <- res %>% 
  filter(cv_type == "stratified") %>% 
  filter(resp == "doc_bathy" & season == "0") %>% 
  select(cv_type, fold, cp_profiles) %>% 
  unnest(cp_profiles)

# Average across folds
pdp_bathy <- prep_pdp(cp = cp_bathy, vars = vars_bathy)

ps11d <- pdp_bathy %>% 
  left_join(vars_bathy %>% rename(var_name = variable), by = join_by(var_name)) %>% 
  mutate(var_name = fct_reorder(var_name, dropout_loss, .desc = TRUE)) %>% 
  ggplot() +
  geom_path(aes(x = x, y = yhat_loc)) +
  geom_ribbon(aes(x = x, ymin = yhat_loc - yhat_spr, ymax = yhat_loc + yhat_spr), alpha = 0.2) +
  labs(y = "Predicted log(DOC)") +
  ps11_facet +
  ps11_theme


ps11 <- ps11a / ps11b / ps11c / ps11d + plot_annotation(tag_levels = "a")
ps11

## Save
ggsave(file = "figures/figure_s11.png", plot = ps11, width = 188, height = 200, unit = "mm", bg = "white")


## Presentation figure: all 4 layers together ----
#--------------------------------------------------------------------------#
# Plot a
pfa <- ggplot(df_1ab %>% filter(lon != -0.5))  +
  geom_sf(data = grat, alpha = 0.8, color = "gray80", linewidth = 0.2) +
  geom_tile(aes(x = lon, y = lat, fill = doc_avg, colour = doc_avg)) +
  geom_sf(data = world_rob, fill = "gray80", colour = NA) +
  ggplot2::scale_fill_viridis_c(option = "F", trans = "log1p") +
  ggplot2::scale_colour_viridis_c(option = "F", trans = "log1p", guide = "none") +
  coord_sf(crs = '+proj=robin +lon_0=180 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs', default_crs = sf::st_crs(4326), datum = NA) +
  labs(x = "Longitude", y = "Latitude", fill = "Predicted DOC (μmol kg<sup>-1</sup>)") +
  p2_theme

# Plot b
pfb <- ggplot(df_2a %>% filter(lon != -0.5))  +
  geom_sf(data = grat, alpha = 0.8, color = "gray80", linewidth = 0.2) +
  geom_tile(aes(x = lon, y = lat, fill = doc_avg, colour = doc_avg)) +
  geom_sf(data = world_rob, fill = "gray80", colour = NA) +
  ggplot2::scale_fill_viridis_c(option = "F", trans = "log1p") +
  ggplot2::scale_colour_viridis_c(option = "F", trans = "log1p", guide = "none") +
  coord_sf(crs = '+proj=robin +lon_0=180 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs', default_crs = sf::st_crs(4326), datum = NA) +
  labs(x = "Longitude", y = "Latitude", fill = "Predicted DOC (μmol kg<sup>-1</sup>)") +
  p2_theme

# Plot c
pfc <- ggplot(df_2b %>% filter(lon != -0.5))  +
  geom_sf(data = grat, alpha = 0.8, color = "gray80", linewidth = 0.2) +
  geom_tile(aes(x = lon, y = lat, fill = doc_avg, colour = doc_avg)) +
  geom_sf(data = world_rob, fill = "gray80", colour = NA) +
  ggplot2::scale_fill_viridis_c(option = "F", trans = "log1p") +
  ggplot2::scale_colour_viridis_c(option = "F", trans = "log1p", guide = "none") +
  coord_sf(crs = '+proj=robin +lon_0=180 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs', default_crs = sf::st_crs(4326), datum = NA) +
  labs(x = "Longitude", y = "Latitude", fill = "Predicted DOC (μmol kg<sup>-1</sup>)") +
  p2_theme

# Plot d
pfd <- ggplot(df_2c %>% filter(lon != -0.5))  +
  geom_sf(data = grat, alpha = 0.8, color = "gray80", linewidth = 0.2) +
  geom_tile(aes(x = lon, y = lat, fill = doc_avg, colour = doc_avg)) +
  geom_sf(data = world_rob, fill = "gray80", colour = NA) +
  ggplot2::scale_fill_viridis_c(option = "F", trans = "log1p") +
  ggplot2::scale_colour_viridis_c(option = "F", trans = "log1p", guide = "none") +
  coord_sf(crs = '+proj=robin +lon_0=180 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs', default_crs = sf::st_crs(4326), datum = NA) +
  labs(x = "Longitude", y = "Latitude", fill = "Predicted DOC (μmol kg<sup>-1</sup>)") +
  p2_theme

pf <- (pfa + pfb) / (pfc + pfd) + plot_layout(axis_titles = "collect") + plot_annotation(tag_levels = "a") & theme(legend.position = 'bottom')
pf
