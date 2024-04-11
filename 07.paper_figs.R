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
      !is.na(log_doc_meso),
      log_doc_meso,
      log_doc_bathy
    )
  )) %>% 
  select(-contains("log_doc")) %>% 
  rsq(truth = truth, estimate = .pred) %>% 
  rename(rsq = .estimate)




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
  legend.title.position = "top"
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

# Plot a
p1a <- ggplot(df_1ab) + 
  geom_polygon(data = world, aes(x = lon, y = lat, group = group), fill = "grey") +
  geom_tile(aes(x = lon, y = lat, fill = doc_avg, colour = doc_avg)) + 
  ggplot2::scale_fill_viridis_c(option = "F", trans = "log1p") +
  ggplot2::scale_colour_viridis_c(option = "F", trans = "log1p", guide = "none") +
  labs(x = "Longitude", y = "Latitude", fill = "Pred. DOC<br>(μmol kg<sup>-1</sup>)") +
  coord_map(projection = "mollweide") +
  p_theme

p1a <- ggplot(df_1ab)  +
  geom_sf(data = world_sf, fill = "gray80", colour = NA) +
  geom_tile(aes(x = lon, y = lat, fill = doc_avg, colour = doc_avg)) +
  ggplot2::scale_fill_viridis_c(option = "F", trans = "log1p") +
  ggplot2::scale_colour_viridis_c(option = "F", trans = "log1p", guide = "none") +
  coord_sf(crs = "+proj=moll +lat_0=20 +lon_0=0", default_crs = 4326) +
  labs(x = "Longitude", y = "Latitude", fill = "Predicted DOC (μmol kg<sup>-1</sup>)") +
  p1_theme




# Plot b
p1b <- ggplot(df_1ab) + 
  geom_polygon(data = world, aes(x = lon, y = lat, group = group), fill = "grey") +
  geom_tile(aes(x = lon, y = lat, fill = doc_sd, colour = doc_sd)) + 
  ggplot2::scale_fill_viridis_c(option = "E", trans = "log1p") +
  ggplot2::scale_colour_viridis_c(option = "E", trans = "log1p", guide = "none") +
  labs(x = "Longitude", y = "Latitude", fill = "Pred. uncert.<br>(μmol kg<sup>-1</sup>)") +
  coord_map(projection = "mollweide") +
  p_theme

p1b <- ggplot(df_1ab)  +
  geom_sf(data = world_sf, fill = "gray80", colour = NA) +
  geom_tile(aes(x = lon, y = lat, fill = doc_sd, colour = doc_sd)) +
  ggplot2::scale_fill_viridis_c(option = "E", trans = "log1p", breaks = c(1, 5, 10, 20, 30)) +
  ggplot2::scale_colour_viridis_c(option = "E", trans = "log1p", guide = "none") +
  coord_sf(crs = "+proj=moll +lat_0=20 +lon_0=0", default_crs = 4326) +
  labs(x = "Longitude", y = "Latitude", fill = "Prediction uncertainty (μmol kg<sup>-1</sup>)") +
  p1_theme



## c
df_1c_other <- res %>% 
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
    .groups = "drop"
  ) %>% 
  # compute amplitude of seasonal cycle
  group_by(lon, lat) %>% 
  # Keep only pixels where DOC is predicted for all seasons
  filter(n() == 4) %>% 
  summarise(
    seas_amp = max(doc_avg, na.rm = TRUE) - min(doc_avg, na.rm = TRUE), 
    seas_var = var(doc_avg, na.rm = TRUE),
    .groups = "drop"
  )

p1c <- ggplot(df_1c_other) + 
  geom_polygon(data = world, aes(x = lon, y = lat, group = group), fill = "grey") +
  geom_tile(aes(x = lon, y = lat, colour = seas_amp, fill = seas_amp)) + 
  ggplot2::scale_fill_viridis_c(trans = "log1p") +
  ggplot2::scale_colour_viridis_c(trans = "log1p", guide = "none") +
  labs(x = "Longitude", y = "Latitude", fill = "Seasonal<br>amplitude<br>(μmol kg<sup>-1</sup>)") +
  coord_map(projection = "mollweide") +
  p_theme

p1c <- ggplot(df_1c_other)  +
  geom_sf(data = world_sf, fill = "gray80", colour = NA) +
  geom_tile(aes(x = lon, y = lat, fill = seas_amp, colour = seas_amp)) +
  ggplot2::scale_fill_viridis_c(trans = "log1p", breaks = c(1, 5, 10, 20, 30)) +
  ggplot2::scale_colour_viridis_c(trans = "log1p", guide = "none") +
  coord_sf(crs = "+proj=moll +lat_0=20 +lon_0=0", default_crs = 4326) +
  labs(x = "Longitude", y = "Latitude", fill = "Seasonal amplitude (μmol kg<sup>-1</sup>)") +
  p1_theme

#df_1c <- res %>% 
#  filter(resp == "doc_surf" & season != "0") %>% 
#  select(fold, new_preds) %>% 
#  unnest(new_preds) %>% 
#  mutate(season = as.character(season)) %>% 
#  # Apply exp to predictions as we predicted log(doc)
#  mutate(pred_doc = exp(pred_doc_log), .after = pred_doc_log) %>% 
#  select(season, fold, lon, lat, contains("doc")) %>% 
#  left_join(rsquares %>% filter(resp == "doc_surf" & season != "0" & cv_type == "stratified"), by = join_by(season, fold)) %>% 
#  # Compute avg prediction per pixel
#  group_by(lon, lat, season) %>% 
#  summarise(
#    doc_avg = wtd.mean(pred_doc, weights = rsq, na.rm = TRUE), 
#    .groups = "drop"
#  ) %>% 
#  # compute amplitude of seasonal cycle
#  group_by(lon, lat) %>% 
#  summarise(
#    seas_up = max(doc_avg, na.rm = TRUE) - mean(doc_avg, na.rm = TRUE), 
#    seas_down = min(doc_avg, na.rm = TRUE) - mean(doc_avg, na.rm = TRUE), 
#    .groups = "drop"
#  )
#
#p1c <- ggplot(df_1c) + 
#  geom_polygon(data = world, aes(x = lon, y = lat, group = group), fill = "grey") +
#  geom_raster(aes(x = lon, y = lat, fill = seas_up)) + 
#  scale_fill_distiller(palette = "Reds", direction = 1) +
#  labs(x = "Longitude", y = "Latitude", fill = "Positive<br>seas. amp.<br>(μmol kg<sup>-1</sup>)") +
#  coord_map(projection = "moll") +
#  scale_xy_map() +
#  p_theme
#
#p1d <- ggplot(df_1c) + 
#  geom_polygon(data = world, aes(x = lon, y = lat, group = group), fill = "grey") +
#  geom_raster(aes(x = lon, y = lat, fill = seas_down)) + 
#  #ggplot2::scale_fill_viridis_c(trans = "log1p") +
#  scale_fill_distiller(palette = "Blues", direction = -1) +
#  labs(x = "Longitude", y = "Latitude", fill = "Negative<br>seas. amp.<br>(μmol kg<sup>-1</sup>)") +
#  coord_map(projection = "moll") +
#  scale_xy_map() +
#  p_theme


## Assemble
#p1 <- p1a + p1b + p1c + p1d + plot_layout(axis_titles = "collect") + plot_annotation(tag_levels = "a")
p1 <- p1a / p1b / p1c + plot_layout(axis_titles = "collect") + plot_annotation(tag_levels = "a")
p1

# Save
ggsave(file = "figures/figure_1.png", plot = p1, width = 180, height = 250, unit = "mm", bg = "white")


## 2: Deep climatologies ----
#--------------------------------------------------------------------------#
# Deep climatologies:
# a - meso climatology
# b - meso uncertainty
# c - bathy climatology
# d - bathy uncertainty

p2_theme <- theme(
  axis.title = element_blank(),
  legend.position = "bottom",
  legend.title = element_markdown(size = 10),
  legend.key.width = unit(1.2, "cm"),
  legend.key.height = unit(2, "mm"),
  legend.direction = "horizontal",
  legend.title.position = "top"
  )

## a + b
# Get projections, average and sd by pixel
df_2ab <- res %>% 
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

## c + d
# Get projections, average and sd by pixel
df_2cd <- res %>% 
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

# Get common colourbar limits for meso and bathy
doc_avg_lims <- c(
  min(df_2ab$doc_avg, df_2cd$doc_avg), 
  max(df_2ab$doc_avg, df_2cd$doc_avg)  
)
doc_sd_lims <- c(
  min(df_2ab$doc_sd, df_2cd$doc_sd), 
  max(df_2ab$doc_sd, df_2cd$doc_sd)  
)

# Plot a
p2a <- ggplot(df_2ab) + 
  geom_polygon(data = world, aes(x = lon, y = lat, group = group), fill = "grey") +
  geom_tile(aes(x = lon, y = lat, fill = doc_avg, colour = doc_avg)) + 
  ggplot2::scale_fill_viridis_c(option = "F", trans = "log1p", limits = doc_avg_lims) +
  ggplot2::scale_colour_viridis_c(option = "F", trans = "log1p", limits = doc_avg_lims, guide = "none") +
  labs(x = "Longitude", y = "Latitude", fill = "Predicted<br>DOC<br>(μmol kg<sup>-1</sup>)") +
  coord_map(projection = "moll") +
  p_theme

p2a <- ggplot(df_2ab)  +
  geom_sf(data = world_sf, fill = "gray80", colour = NA) +
  geom_tile(aes(x = lon, y = lat, fill = doc_avg, colour = doc_avg)) +
  ggplot2::scale_fill_viridis_c(option = "F", trans = "log1p", limits = doc_avg_lims) +
  ggplot2::scale_colour_viridis_c(option = "F", trans = "log1p", guide = "none", limits = doc_avg_lims) +
  coord_sf(crs = "+proj=moll +lat_0=20 +lon_0=0", default_crs = 4326) +
  labs(x = "Longitude", y = "Latitude", fill = "Predicted DOC (μmol kg<sup>-1</sup>)") +
  p2_theme


# Plot b
p2b <- ggplot(df_2ab) + 
  geom_polygon(data = world, aes(x = lon, y = lat, group = group), fill = "grey") +
  geom_tile(aes(x = lon, y = lat, fill = doc_sd, colour = doc_sd)) + 
  ggplot2::scale_fill_viridis_c(option = "E", trans = "log1p", limits = doc_sd_lims) +
  ggplot2::scale_colour_viridis_c(option = "E", trans = "log1p", limits = doc_sd_lims, guide = "none") +
  labs(x = "Longitude", y = "Latitude", fill = "Prediction<br>uncertainty<br>(μmol kg<sup>-1</sup>)") +
  coord_map(projection = "moll") +
  p_theme

p2b <- ggplot(df_2ab)  +
  geom_sf(data = world_sf, fill = "gray80", colour = NA) +
  geom_tile(aes(x = lon, y = lat, fill = doc_sd, colour = doc_sd)) +
  ggplot2::scale_fill_viridis_c(option = "E", trans = "log1p", limits = doc_sd_lims, breaks = c(1, 2, 4, 6)) +
  ggplot2::scale_colour_viridis_c(option = "E", trans = "log1p", limits = doc_sd_lims, guide = "none") +
  coord_sf(crs = "+proj=moll +lat_0=20 +lon_0=0", default_crs = 4326) +
  labs(x = "Longitude", y = "Latitude", fill = "Prediction uncertainty (μmol kg<sup>-1</sup>)") +
  p2_theme


# Plot c
p2c <- ggplot(df_2cd) + 
  geom_polygon(data = world, aes(x = lon, y = lat, group = group), fill = "grey") +
  geom_tile(aes(x = lon, y = lat, fill = doc_avg, colour = doc_avg)) + 
  ggplot2::scale_fill_viridis_c(option = "F", trans = "log1p", limits = doc_avg_lims) +
  ggplot2::scale_colour_viridis_c(option = "F", trans = "log1p", limits = doc_avg_lims, guide = "none") +
  labs(x = "Longitude", y = "Latitude", fill = "Predicted<br>DOC<br>(μmol kg<sup>-1</sup>)") +
  coord_map(projection = "moll") +
  p_theme

p2c <- ggplot(df_2cd)  +
  geom_sf(data = world_sf, fill = "gray80", colour = NA) +
  geom_tile(aes(x = lon, y = lat, fill = doc_avg, colour = doc_avg)) +
  ggplot2::scale_fill_viridis_c(option = "F", trans = "log1p", limits = doc_avg_lims) +
  ggplot2::scale_colour_viridis_c(option = "F", trans = "log1p", guide = "none", limits = doc_avg_lims) +
  coord_sf(crs = "+proj=moll +lat_0=20 +lon_0=0", default_crs = 4326) +
  labs(x = "Longitude", y = "Latitude", fill = "Predicted DOC (μmol kg<sup>-1</sup>)") +
  p2_theme


# Plot d
p2d <- ggplot(df_2cd) + 
  geom_polygon(data = world, aes(x = lon, y = lat, group = group), fill = "grey") +
  geom_tile(aes(x = lon, y = lat, fill = doc_sd, colour = doc_sd)) + 
  ggplot2::scale_fill_viridis_c(option = "E", trans = "log1p", limits = doc_sd_lims) +
  ggplot2::scale_colour_viridis_c(option = "E", trans = "log1p", limits = doc_sd_lims, guide = "none") +
  labs(x = "Longitude", y = "Latitude", fill = "Prediction<br>uncertainty<br>(μmol kg<sup>-1</sup>)") +
  coord_map(projection = "moll") +
  p_theme

p2d <- ggplot(df_2cd)  +
  geom_sf(data = world_sf, fill = "gray80", colour = NA) +
  geom_tile(aes(x = lon, y = lat, fill = doc_sd, colour = doc_sd)) +
  ggplot2::scale_fill_viridis_c(option = "E", trans = "log1p", limits = doc_sd_lims, breaks = c(1, 2, 4, 6)) +
  ggplot2::scale_colour_viridis_c(option = "E", trans = "log1p", limits = doc_sd_lims, guide = "none") +
  coord_sf(crs = "+proj=moll +lat_0=20 +lon_0=0", default_crs = 4326) +
  labs(x = "Longitude", y = "Latitude", fill = "Prediction uncertainty (μmol kg<sup>-1</sup>)") +
  p2_theme

## Assemble
p2 <- p2a + p2b + p2c + p2d + plot_annotation(tag_levels = "a") + plot_layout(guides = "collect") & theme(legend.position = 'bottom', legend.spacing.x = unit(2.7, "cm"))
#p2 <- p2a + p2b + p2c + p2d + plot_annotation(tag_levels = "a") + plot_layout(guides = "collect") & theme(legend.position = 'bottom', legend.justification = "center")
p2

## Save
ggsave(file = "figures/figure_2.png", plot = p2, width = 180, height = 130, unit = "mm", bg = "white")



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

## 1st version with patchwork
# Surface
vip_surf_ann <- full_vip %>% 
  filter(variable != "_full_model_") %>%
  filter(resp == "doc_surf" & season == "0") %>% 
  mutate(
    variable = forcats::fct_drop(variable),
    variable = forcats::fct_reorder(variable, dropout_loss),
    variable = forcats::fct_other(variable, keep = tail(levels(variable), n = 5), other_level = "other"),
    variable = forcats::fct_reorder(variable, dropout_loss)
    )

p3a <- vip_surf_ann %>% 
  ggplot() + 
  geom_vline(data = full_vip %>% filter(variable == "_full_model_"), aes(xintercept = mean(dropout_loss)), colour = "grey", linewidth = 2) +
  geom_boxplot(aes(x = dropout_loss, y = variable)) +
  expand_limits(x = 0)

# Meso
vip_meso_ann <- full_vip %>% 
  filter(variable != "_full_model_") %>%
  filter(resp == "doc_meso" & season == "0") %>% 
  mutate(
    variable = forcats::fct_drop(variable),
    variable = forcats::fct_reorder(variable, dropout_loss),
    variable = forcats::fct_other(variable, keep = tail(levels(variable), n = 5), other_level = "other"),
    variable = forcats::fct_reorder(variable, dropout_loss)
  )

p3b <- vip_meso_ann %>% 
  ggplot() + 
  geom_vline(data = full_vip %>% filter(variable == "_full_model_"), aes(xintercept = mean(dropout_loss)), colour = "grey", linewidth = 2) +
  geom_boxplot(aes(x = dropout_loss, y = variable)) +
  expand_limits(x = 0)

# Bathy
vip_bathy_ann <- full_vip %>% 
  filter(variable != "_full_model_") %>%
  filter(resp == "doc_bathy" & season == "0") %>% 
  mutate(
    variable = forcats::fct_drop(variable),
    variable = forcats::fct_reorder(variable, dropout_loss),
    variable = forcats::fct_other(variable, keep = tail(levels(variable), n = 5), other_level = "other"),
    variable = forcats::fct_reorder(variable, dropout_loss)
  )

p3c <- vip_bathy_ann %>% 
  ggplot() + 
  geom_vline(data = full_vip %>% filter(variable == "_full_model_"), aes(xintercept = mean(dropout_loss)), colour = "grey", linewidth = 2) +
  geom_boxplot(aes(x = dropout_loss, y = variable)) +
  expand_limits(x = 0)

# Assemble
p3 <- p3a / p3b / p3c + plot_layout(axis_titles = "collect")
p3

## 2nd version with facet_wrap
full_model_vip <- full_vip %>% 
  filter(variable == "_full_model_") %>% 
  mutate(resp = factor(resp, levels = c("doc_surf", "doc_meso", "doc_bathy"))) %>% 
  group_by(resp) %>% 
  summarise(dropout_loss = mean(dropout_loss))


df_3 <- vip_surf_ann %>% 
  bind_rows(vip_meso_ann) %>% 
  bind_rows(vip_bathy_ann) %>% 
  mutate(resp = fct_inorder(resp))

# nice formatting for variables
df_3 <- df_3 %>% 
  mutate(
    variable = str_replace_all(variable, "_", " "),
    variable = str_to_sentence(variable),
    variable = str_replace_all(variable, "Mld", "MLD"),
    variable = str_replace_all(variable, "Bpp", "BBP"),
    variable = str_replace_all(variable, "Fmicro", "F micro"),
    variable = str_replace_all(variable, "surf", "surf."),
    variable = str_replace_all(variable, "meso", "meso."),
    variable = str_replace_all(variable, "bathy", "bathy.")
  ) %>% 
  arrange(dropout_loss) %>% 
  mutate(variable = fct_inorder(variable))

p3 <- ggplot(df_3) + 
  geom_vline(data = full_model_vip, aes(xintercept = dropout_loss), colour = "grey", linewidth = 2) +
  geom_boxplot(aes(x = dropout_loss, y = variable, colour = resp), outlier.size = 0.3) +
  scale_colour_manual(
    values = c("#bdd7e7", "#6baed6", "#2171b5"),
    labels = c(
      doc_surf = "Surf.",
      doc_meso = "Meso.",
      doc_bathy = "Bathy."
    )) +
  expand_limits(x = 0) +
  labs(x = "Dropout loss", y = "Variable", colour = "Layer") +
  facet_wrap(~resp, ncol = 1, scales = "free_y") +
  theme(
    strip.background = element_blank(), strip.text.x = element_blank(), # empty facets
    legend.position = "inside", legend.position.inside = c(0.75, 0.5), # legend position
    legend.background = element_rect(fill = "white", colour = NA)
  )

p3


## S1: Map of observations ----
#--------------------------------------------------------------------------#
load("data/02.ann_surf.Rdata")
load("data/02.ann_meso.Rdata")
load("data/02.ann_bathy.Rdata")

ps1_theme <- theme(
  axis.title = element_blank()
)

doc_lims <- c(
  min(df_ann_surf_fit$doc_surf, df_ann_meso_fit$doc_meso, df_ann_bathy_fit$doc_bathy),
  max(df_ann_surf_fit$doc_surf, df_ann_meso_fit$doc_meso, df_ann_bathy_fit$doc_bathy)
)

# Surface
ps1a <-   df_ann_surf_fit %>% 
  ggplot() + 
  geom_polygon(data = world, aes(x = lon, y = lat, group = group), fill = "grey") +
  geom_point(aes(x = lon, y = lat), size = 0.5, shape = 1, alpha = 0.5) +
  coord_map(projection = "mollweide") +
  theme(
    legend.title = element_markdown()
  )

ps1a <- ggplot(df_ann_surf_fit)  +
  geom_sf(data = world_sf, fill = "gray80", colour = NA) +
  geom_point(aes(x = lon, y = lat), size = 0.5, alpha = 0.5) +
  coord_sf(crs = "+proj=moll +lat_0=20 +lon_0=0", default_crs = 4326) +
  labs(x = "Longitude", y = "Latitude") +
  ps1_theme


# Meso
ps1b <- df_ann_meso_fit %>% 
  ggplot() + 
  geom_polygon(data = world, aes(x = lon, y = lat, group = group), fill = "grey") +
  geom_point(aes(x = lon, y = lat), size = 0.5, shape = 1, alpha = 0.5) +
  coord_map(projection = "mollweide") +
  theme(
    legend.title = element_markdown()
  )

ps1b <- ggplot(df_ann_meso_fit)  +
  geom_sf(data = world_sf, fill = "gray80", colour = NA) +
  geom_point(aes(x = lon, y = lat), size = 0.5, alpha = 0.5) +
  coord_sf(crs = "+proj=moll +lat_0=20 +lon_0=0", default_crs = 4326) +
  labs(x = "Longitude", y = "Latitude") +
  ps1_theme

# Bathy
ps1c <- df_ann_bathy_fit %>% 
  ggplot() + 
  geom_polygon(data = world, aes(x = lon, y = lat, group = group), fill = "grey") +
  geom_point(aes(x = lon, y = lat), size = 0.5, shape = 1, alpha = 0.5) +
  coord_map(projection = "mollweide") +
  theme(
    legend.title = element_markdown()
  )

ps1c <- ggplot(df_ann_bathy_fit)  +
  geom_sf(data = world_sf, fill = "gray80", colour = NA) +
  geom_point(aes(x = lon, y = lat), size = 0.5, alpha = 0.5) +
  coord_sf(crs = "+proj=moll +lat_0=20 +lon_0=0", default_crs = 4326) +
  labs(x = "Longitude", y = "Latitude") +
  ps1_theme

ps1 <- ps1a + ps1b + ps1c + plot_layout(ncol = 2, guides = "collect") + plot_annotation(tag_levels = "a")
ps1


## S2: Distribution of predictors ----
#--------------------------------------------------------------------------#
# Surface
pred_surf <- bind_rows(
  df_ann_surf_fit %>% select(temperature_surf:par) %>% pivot_longer(everything(), names_to = "variable") %>% mutate(type = "Learning"),
  df_ann_surf_pred %>% select(temperature_surf:par) %>% pivot_longer(everything(), names_to = "variable") %>% mutate(type = "New")
)
  
ps2a <- ggplot(pred_surf) + 
  geom_density(aes(x = value, linetype = type)) +
  labs(x = "Value", y = "Density", linetype = "Data type") +
  #scale_colour_manual(values = c("grey30", "black")) +
  facet_wrap(~variable, scales = "free", ncol = 5) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Meso
pred_meso <- bind_rows(
  df_ann_meso_fit %>% select(temperature_surf:nitrate_meso) %>% pivot_longer(everything(), names_to = "variable") %>% mutate(type = "Learning"),
  df_ann_meso_pred %>% select(temperature_surf:nitrate_meso) %>% pivot_longer(everything(), names_to = "variable") %>% mutate(type = "New")
) %>% 
  filter(!variable %in% pred_surf$variable) # drop predictors already shown for surface layer

ps2b <- ggplot(pred_meso) + 
  geom_density(aes(x = value, linetype = type)) +
  labs(x = "Value", y = "Density", linetype = "Data type") +
  facet_wrap(~variable, scales = "free", ncol = 5) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Bathy
pred_bathy <- bind_rows(
  df_ann_bathy_fit %>% select(temperature_surf:nitrate_meso) %>% pivot_longer(everything(), names_to = "variable") %>% mutate(type = "Learning"),
  df_ann_bathy_pred %>% select(temperature_surf:nitrate_meso) %>% pivot_longer(everything(), names_to = "variable") %>% mutate(type = "New")
) %>% 
  filter(!variable %in% pred_surf$variable) %>% # drop predictors already shown for surface layer
  filter(!variable %in% pred_meso$variable) # and those already shown for meso layer

ps2c <- ggplot(pred_bathy) + 
  geom_density(aes(x = value, linetype = type)) +
  labs(x = "Value", y = "Density", linetype = "Data type") +
  facet_wrap(~variable, scales = "free", ncol = 5) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Assemble
ps2 <- ps2a / ps2b / ps2c + plot_annotation(tag_levels = "a") + plot_layout(guides = "collect", heights = c(5, 2, 1))
ps2


## S3: Rsquares dist in all layers ----
#--------------------------------------------------------------------------#
ps3a <- rsquares %>% 
  left_join(resp_to_layer, by = join_by(resp)) %>% 
  filter(season == "0") %>% 
  mutate(bp_factor = interaction(layer, cv_type)) %>% 
  ggplot() +
  geom_boxplot(aes(x = layer, y = rsq, group = bp_factor, colour = cv_type), position = position_dodge(1)) +
  scale_colour_manual(values = c("grey", "black")) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  labs(x = "Layer", y = "R²", colour = "CV type")


ps3b <- rsquares %>% 
  left_join(resp_to_layer, by = join_by(resp)) %>% 
  filter(season != "0") %>% 
  mutate(bp_factor = interaction(layer, season)) %>% 
  ggplot() +
  geom_boxplot(aes(x = layer, y = rsq, group = bp_factor, colour = season), position = position_dodge(1)) +
  scale_colour_manual(
    #values = c("#3F9D86", "#486E9E", "#F2BB62", "#F29762"),
    values = c( "#F29762", "#3F9D86", "#486E9E", "#F2BB62"),
    labels = c(
      `1` = "winter",
      `2` = "spring",
      `3` = "summer",
      `4` = "autumn"
    )
    ) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  labs(x = "Layer", y = "R²", colour = "Season")

ps3 <- ps3a + ps3b + plot_annotation(tag_levels = "a")
ps3


## S4: Spatial VS stratified in surface layer ----
#--------------------------------------------------------------------------#
# - rsquares value
# - diff in projections

ps4_theme <- theme(
  axis.title = element_blank(),
  legend.position = "bottom",
  legend.title = element_markdown(size = 10),
  legend.key.width = unit(2.5, "cm"),
  legend.key.height = unit(2, "mm"),
  legend.direction = "horizontal",
  legend.title.position = "top"
)


# Get projections, average and sd by pixel for stratified
df_4sa <- res %>% 
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
df_4sb <- res %>% 
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
df_4s <- df_4sa %>% 
  rename(strat = doc_avg) %>% 
  left_join(df_4sb %>% rename(spat = doc_avg), by = join_by(lon, lat)) %>% 
  mutate(diff = strat - spat)

ps4 <- ggplot(df_4s) + 
  geom_polygon(data = world, aes(x = lon, y = lat, group = group), fill = "grey") +
  geom_tile(aes(x = lon, y = lat, fill = diff, colour = diff)) + 
  scale_fill_gradient2(low = "#4575b4", mid = "#ffffbf", high = "#d73027") +
  scale_colour_gradient2(low = "#4575b4", mid = "#ffffbf", high = "#d73027", guide = "none") +
  labs(x = "Longitude", y = "Latitude", fill = "strat - spat<br>(μmol kg<sup>-1</sup>)") +
  coord_map(projection = "moll") +
  theme(legend.title = element_markdown())


ps4 <- ggplot(df_4s)  +
  geom_sf(data = world_sf, fill = "gray80", colour = NA) +
  geom_tile(aes(x = lon, y = lat, fill = diff, colour = diff)) + 
  scale_fill_gradient2(low = "#4575b4", mid = "#ffffbf", high = "#d73027") +
  scale_colour_gradient2(low = "#4575b4", mid = "#ffffbf", high = "#d73027", guide = "none") +
  coord_sf(crs = "+proj=moll +lat_0=20 +lon_0=0", default_crs = 4326) +
  labs(x = "Longitude", y = "Latitude", fill = "Difference between projection from stratified CV and spatial CV (μmol kg<sup>-1</sup>)") +
  ps4_theme

# Save
ggsave(file = "figures/figure_s4.png", plot = ps4, width = 180, height = 90, unit = "mm", bg = "white")


## S5: Seasonal projections in the surface layer ----
#--------------------------------------------------------------------------#
df_s5 <- res %>% 
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
    .groups = "drop"
  ) 

doc_avg_lims <- c(min(df_s5$doc_avg), max(df_s5$doc_avg))

ps5a <- ggplot(df_s5 %>% filter(season == 1)) + 
  geom_sf(data = world_sf, fill = "gray80", colour = NA) +
  geom_tile(aes(x = lon, y = lat, fill = doc_avg, colour = doc_avg)) +
  ggplot2::scale_fill_viridis_c(option = "F", trans = "log1p", na.value = NA, limits = doc_avg_lims) +
  ggplot2::scale_colour_viridis_c(option = "F", trans = "log1p", guide = "none", na.value = NA, limits = doc_avg_lims) +
  coord_sf(crs = "+proj=moll +lat_0=20 +lon_0=0", default_crs = 4326) +
  labs(x = "Longitude", y = "Latitude", fill = "Predicted DOC (μmol kg<sup>-1</sup>)") +
  p2_theme

ps5b <- ggplot(df_s5 %>% filter(season == 2)) + 
  geom_sf(data = world_sf, fill = "gray80", colour = NA) +
  geom_tile(aes(x = lon, y = lat, fill = doc_avg, colour = doc_avg)) +
  ggplot2::scale_fill_viridis_c(option = "F", trans = "log1p", na.value = NA, limits = doc_avg_lims) +
  ggplot2::scale_colour_viridis_c(option = "F", trans = "log1p", guide = "none", na.value = NA, limits = doc_avg_lims) +
  coord_sf(crs = "+proj=moll +lat_0=20 +lon_0=0", default_crs = 4326) +
  labs(x = "Longitude", y = "Latitude", fill = "Predicted DOC (μmol kg<sup>-1</sup>)") +
  p2_theme

ps5c <- ggplot(df_s5 %>% filter(season == 3)) + 
  geom_sf(data = world_sf, fill = "gray80", colour = NA) +
  geom_tile(aes(x = lon, y = lat, fill = doc_avg, colour = doc_avg)) +
  ggplot2::scale_fill_viridis_c(option = "F", trans = "log1p", na.value = NA, limits = doc_avg_lims) +
  ggplot2::scale_colour_viridis_c(option = "F", trans = "log1p", guide = "none", na.value = NA, limits = doc_avg_lims) +
  coord_sf(crs = "+proj=moll +lat_0=20 +lon_0=0", default_crs = 4326) +
  labs(x = "Longitude", y = "Latitude", fill = "Predicted DOC (μmol kg<sup>-1</sup>)") +
  p2_theme

ps5d <- ggplot(df_s5 %>% filter(season == 4)) + 
  geom_sf(data = world_sf, fill = "gray80", colour = NA) +
  geom_tile(aes(x = lon, y = lat, fill = doc_avg, colour = doc_avg)) +
  ggplot2::scale_fill_viridis_c(option = "F", trans = "log1p", na.value = NA, limits = doc_avg_lims) +
  ggplot2::scale_colour_viridis_c(option = "F", trans = "log1p", guide = "none", na.value = NA, limits = doc_avg_lims) +
  coord_sf(crs = "+proj=moll +lat_0=20 +lon_0=0", default_crs = 4326) +
  labs(x = "Longitude", y = "Latitude", fill = "Predicted DOC (μmol kg<sup>-1</sup>)") +
  p2_theme

ps5 <- ps5a + ps5b + ps5c + ps5d + plot_annotation(tag_levels = "a") + plot_layout(guides = "collect") & theme(legend.position = 'bottom')

## Save
ggsave(file = "figures/figure_s5.png", plot = ps5, width = 180, height = 130, unit = "mm", bg = "white")


## S6: Partial dependance plots ----
#--------------------------------------------------------------------------#
## Surface
# Variables to use
vars_surf <- vip_surf_ann %>% 
  filter(variable != "other") %>% 
  group_by(cv_type, variable) %>% 
  summarise(dropout_loss = mean(dropout_loss))

# CP profiles
cp_surf <- res %>% 
  filter(cv_type == "stratified") %>% 
  filter(resp == "doc_surf" & season == "0") %>% 
  select(cv_type, fold, cp_profiles) %>% 
  unnest(cp_profiles)

# Average across folds
pdp_surf <- prep_pdp(cp = cp_surf, vars = vars_surf)

ps6a <- pdp_surf %>% 
  left_join(vars_surf %>% rename(var_name = variable), by = join_by(var_name)) %>% 
  mutate(var_name = fct_reorder(var_name, dropout_loss, .desc = TRUE)) %>% 
  ggplot() +
  geom_path(aes(x = x, y = yhat_loc)) +
  geom_ribbon(aes(x = x, ymin = yhat_loc - yhat_spr, ymax = yhat_loc + yhat_spr), alpha = 0.2) +
  labs(y = "Predicted log(DOC)") +
  facet_wrap(~var_name, scales = "free_x", strip.position = "bottom") +
  theme(axis.title.x = element_blank(), strip.placement = "outside")


## Meso
# Variables to use
vars_meso <- vip_meso_ann %>% 
  filter(variable != "other") %>% 
  group_by(cv_type, variable) %>% 
  summarise(dropout_loss = mean(dropout_loss))

# CP profiles
cp_meso <- res %>% 
  filter(cv_type == "stratified") %>% 
  filter(resp == "doc_meso" & season == "0") %>% 
  select(cv_type, fold, cp_profiles) %>% 
  unnest(cp_profiles)

# Average across folds
pdp_meso <- prep_pdp(cp = cp_meso, vars = vars_meso)

ps6b <- pdp_meso %>% 
  left_join(vars_meso %>% rename(var_name = variable), by = join_by(var_name)) %>% 
  mutate(var_name = fct_reorder(var_name, dropout_loss, .desc = TRUE)) %>% 
  ggplot() +
  geom_path(aes(x = x, y = yhat_loc)) +
  geom_ribbon(aes(x = x, ymin = yhat_loc - yhat_spr, ymax = yhat_loc + yhat_spr), alpha = 0.2) +
  labs(y = "Predicted log(DOC)") +
  facet_wrap(~var_name, scales = "free_x", strip.position = "bottom") +
  theme(axis.title.x = element_blank(), strip.placement = "outside")


## Bathy
# Variables to use
vars_bathy <- vip_bathy_ann %>% 
  filter(variable != "other") %>% 
  group_by(cv_type, variable) %>% 
  summarise(dropout_loss = mean(dropout_loss))

# CP profiles
cp_bathy <- res %>% 
  filter(cv_type == "stratified") %>% 
  filter(resp == "doc_bathy" & season == "0") %>% 
  select(cv_type, fold, cp_profiles) %>% 
  unnest(cp_profiles)

# Average across folds
pdp_bathy <- prep_pdp(cp = cp_bathy, vars = vars_bathy)

ps6c <- pdp_bathy %>% 
  left_join(vars_bathy %>% rename(var_name = variable), by = join_by(var_name)) %>% 
  mutate(var_name = fct_reorder(var_name, dropout_loss, .desc = TRUE)) %>% 
  ggplot() +
  geom_path(aes(x = x, y = yhat_loc)) +
  geom_ribbon(aes(x = x, ymin = yhat_loc - yhat_spr, ymax = yhat_loc + yhat_spr), alpha = 0.2) +
  labs(y = "Predicted log(DOC)") +
  facet_wrap(~var_name, scales = "free_x", strip.position = "bottom") +
  theme(axis.title.x = element_blank(), strip.placement = "outside")


ps6 <- ps6a / ps6b / ps6c + plot_annotation(tag_levels = "a")
ps6