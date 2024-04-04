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
  bind_rows()

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

## a
# Get rsquares of folds for weighted average and sd
rsquares <- res %>% 
  filter(resp == "doc_surf" & season == "0") %>% 
  filter(cv_type == "stratified") %>% 
  select(fold, preds) %>% 
  unnest(preds) %>% 
  group_by(fold) %>%
  rsq(truth = log_doc_surf, estimate = .pred) %>% 
  rename(rsq = .estimate)

# Get projections, average and sd by pixel
df_1ab <- res %>% 
  filter(resp == "doc_surf" & season == "0") %>% 
  filter(cv_type == "stratified") %>% 
  select(fold, new_preds) %>% 
  unnest(new_preds) %>% 
  # Apply exp to predictions as we predicted log(doc)
  mutate(pred_doc = exp(pred_doc_log), .after = pred_doc_log) %>% 
  select(fold, lon, lat, contains("doc")) %>% 
  left_join(rsquares, by = join_by(fold)) %>% 
  group_by(lon, lat) %>% 
  summarise(
    doc_avg = wtd.mean(pred_doc, weights = rsq, na.rm = TRUE), 
    doc_sd = sqrt(wtd.var(pred_doc, weights = rsq, na.rm = TRUE)), 
    .groups = "drop"
  )

# Plot
p1a <- ggplot(df_1ab) + 
  geom_polygon(data = world, aes(x = lon, y = lat, group = group), fill = "grey") +
  geom_raster(aes(x = lon, y = lat, fill = doc_avg)) + 
  ggplot2::scale_fill_viridis_c(option = "F", trans = "log1p") +
  labs(x = "Longitude", y = "Latitude", fill = "Predicted<br>DOC<br>(μmol kg<sup>-1</sup>)") +
  coord_quickmap(expand = 0) +
  scale_xy_map() +
  theme(legend.title = element_markdown())

## b
# Plot
p1b <- ggplot(df_1ab) + 
  geom_polygon(data = world, aes(x = lon, y = lat, group = group), fill = "grey") +
  geom_raster(aes(x = lon, y = lat, fill = doc_sd)) + 
  ggplot2::scale_fill_viridis_c(option = "E", trans = "log1p") +
  labs(x = "Longitude", y = "Latitude", fill = "Prediction<br>uncertainty<br>(μmol kg<sup>-1</sup>)") +
  coord_quickmap(expand = 0) +
  scale_xy_map() +
  theme(legend.title = element_markdown())

## c
res %>% 
  filter(resp == "doc_surf" & season != "0") %>% 
  select(fold, new_preds) %>% 
  unnest(new_preds) %>% 
  # Apply exp to predictions as we predicted log(doc)
  mutate(pred_doc = exp(pred_doc_log), .after = pred_doc_log) %>% 
  select(season, fold, lon, lat, contains("doc"))


## Assemble
p1 <- p1a + p1b
p1
# Save