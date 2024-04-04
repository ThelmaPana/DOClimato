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

## Prepare all projections ----
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


## Figure 1 ----
#--------------------------------------------------------------------------#
# Surface climatology:
# a - annual climatology
# b - uncertainty
# c - amplitude of seasonal cycle

p1_theme <- theme(
    legend.title = element_markdown(size = 8), 
    text = element_text(size = 8),
    #axis.title.x = element_text(size = 6),
    #axis.title.y = element_text(size = 6), 
    plot.margin = margin(0, 0, 0, 0, "pt"),
    legend.key.height = unit(0.4, "cm"),
    legend.key.width = unit(0.3, "cm")
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
  geom_raster(aes(x = lon, y = lat, fill = doc_avg)) + 
  ggplot2::scale_fill_viridis_c(option = "F", trans = "log1p") +
  labs(x = "Longitude", y = "Latitude", fill = "Pred. DOC<br>(μmol kg<sup>-1</sup>)") +
  coord_quickmap(expand = 0) +
  scale_xy_map() +
  p1_theme


# Plot b
p1b <- ggplot(df_1ab) + 
  geom_polygon(data = world, aes(x = lon, y = lat, group = group), fill = "grey") +
  geom_raster(aes(x = lon, y = lat, fill = doc_sd)) + 
  ggplot2::scale_fill_viridis_c(option = "E", trans = "log1p") +
  labs(x = "Longitude", y = "Latitude", fill = "Pred. uncert.<br>(μmol kg<sup>-1</sup>)") +
  coord_quickmap(expand = 0) +
  scale_xy_map() +
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
  summarise(
    seas_amp = max(doc_avg, na.rm = TRUE) - min(doc_avg, na.rm = TRUE), 
    seas_var = var(doc_avg, na.rm = TRUE),
    .groups = "drop"
  )

p1c_other <- ggplot(df_1c) + 
  geom_polygon(data = world, aes(x = lon, y = lat, group = group), fill = "grey") +
  geom_raster(aes(x = lon, y = lat, fill = seas_amp)) + 
  ggplot2::scale_fill_viridis_c(trans = "log1p") +
  labs(x = "Longitude", y = "Latitude", fill = "Seasonal<br>amplitude<br>(μmol kg<sup>-1</sup>)") +
  coord_quickmap(expand = 0) +
  scale_xy_map() +
  theme(legend.title = element_markdown(), text = element_text(size = 10, family = "Helvetica"))


df_1c <- res %>% 
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
  summarise(
    seas_up = max(doc_avg, na.rm = TRUE) - mean(doc_avg, na.rm = TRUE), 
    seas_down = min(doc_avg, na.rm = TRUE) - mean(doc_avg, na.rm = TRUE), 
    .groups = "drop"
  )

p1c <- ggplot(df_1c) + 
  geom_polygon(data = world, aes(x = lon, y = lat, group = group), fill = "grey") +
  geom_raster(aes(x = lon, y = lat, fill = seas_up)) + 
  scale_fill_distiller(palette = "Reds", direction = 1) +
  labs(x = "Longitude", y = "Latitude", fill = "Positive<br>seas. amp.<br>(μmol kg<sup>-1</sup>)") +
  coord_quickmap(expand = 0) +
  scale_xy_map() +
  p1_theme

p1d <- ggplot(df_1c) + 
  geom_polygon(data = world, aes(x = lon, y = lat, group = group), fill = "grey") +
  geom_raster(aes(x = lon, y = lat, fill = seas_down)) + 
  #ggplot2::scale_fill_viridis_c(trans = "log1p") +
  scale_fill_distiller(palette = "Blues", direction = -1) +
  labs(x = "Longitude", y = "Latitude", fill = "Negative<br>seas. amp.<br>(μmol kg<sup>-1</sup>)") +
  coord_quickmap(expand = 0) +
  scale_xy_map() +
  p1_theme


## Assemble
p1 <- p1a + p1b + p1c + p1d + plot_layout(axis_titles = "collect") + plot_annotation(tag_levels = "a")
p1

# Save
ggsave(file = "figures/figure_1.png", plot = p1, width = 180, height = 90, unit = "mm", bg = "white")


## Figure 2 ----
#--------------------------------------------------------------------------#
# Deep climatologies:
# a - meso climatology
# b - meso uncertainty
# c - bathy climatology
# d - bathy uncertainty


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
  geom_raster(aes(x = lon, y = lat, fill = doc_avg)) + 
  ggplot2::scale_fill_viridis_c(option = "F", trans = "log1p", limits = doc_avg_lims) +
  labs(x = "Longitude", y = "Latitude", fill = "Predicted<br>DOC<br>(μmol kg<sup>-1</sup>)") +
  coord_quickmap(expand = 0) +
  scale_xy_map() +
  p1_theme

# Plot b
p2b <- ggplot(df_2ab) + 
  geom_polygon(data = world, aes(x = lon, y = lat, group = group), fill = "grey") +
  geom_raster(aes(x = lon, y = lat, fill = doc_sd)) + 
  ggplot2::scale_fill_viridis_c(option = "E", trans = "log1p", limits = doc_sd_lims) +
  labs(x = "Longitude", y = "Latitude", fill = "Prediction<br>uncertainty<br>(μmol kg<sup>-1</sup>)") +
  coord_quickmap(expand = 0) +
  scale_xy_map() +
  p1_theme


# Plot c
p2c <- ggplot(df_2cd) + 
  geom_polygon(data = world, aes(x = lon, y = lat, group = group), fill = "grey") +
  geom_raster(aes(x = lon, y = lat, fill = doc_avg)) + 
  ggplot2::scale_fill_viridis_c(option = "F", trans = "log1p", limits = doc_avg_lims) +
  labs(x = "Longitude", y = "Latitude", fill = "Predicted<br>DOC<br>(μmol kg<sup>-1</sup>)") +
  coord_quickmap(expand = 0) +
  scale_xy_map() +
  p1_theme

# Plot d
p2d <- ggplot(df_2cd) + 
  geom_polygon(data = world, aes(x = lon, y = lat, group = group), fill = "grey") +
  geom_raster(aes(x = lon, y = lat, fill = doc_sd)) + 
  ggplot2::scale_fill_viridis_c(option = "E", trans = "log1p", limits = doc_sd_lims) +
  labs(x = "Longitude", y = "Latitude", fill = "Prediction<br>uncertainty<br>(μmol kg<sup>-1</sup>)") +
  coord_quickmap(expand = 0) +
  scale_xy_map() +
  p1_theme

## Assemble
p2 <- p2a + p2b + p2c + p2d + plot_annotation(tag_levels = "a")
p2

## Save
ggsave(file = "figures/figure_2.png", plot = p2, width = 180, height = 90, unit = "mm", bg = "white")


## Figure 3 ----
#--------------------------------------------------------------------------#
# Variable importance
# Unnest variable importance

full_vip_surf_ann <- res %>% 
  filter(cv_type == "stratified") %>% 
  select(resp, season, fold, importance) %>% 
  unnest(importance) %>%
  mutate(variable = forcats::fct_reorder(variable, dropout_loss)) %>%
  group_by(resp, season, fold, variable) %>%
  summarise(dropout_loss = mean(dropout_loss), .groups = "drop")

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
  geom_boxplot(aes(x = dropout_loss, y = variable))


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
  geom_boxplot(aes(x = dropout_loss, y = variable))



p3a / p3b
