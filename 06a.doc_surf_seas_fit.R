## ----set_up--------------------------------------------------------------------------------------------------------------------------------------------------------------------
#|output: false
#|cache: false
source("utils.R")
load("data/02.seas_surf.Rdata")
output_filename <- "data/06.doc_surf_seas_pred.Rdata"

df_fit <- df_seas_surf_fit
df_pred <- df_seas_surf_pred


## ----roles---------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Response variable
resp_var <- c("log_doc_surf", "doc_surf") # keep both untransformed and transformed predictor
# Explanatory variables
exp_vars <- df_fit %>% select(temperature_surf:par) %>% colnames()
# Metadata
meta_vars <- c("lon", "lat", "season")


## ----resp_dist-----------------------------------------------------------------------------------------------------------------------------------------------------------------
ggplot(df_fit) + geom_histogram(aes(x = log_doc_surf), bins = 100) + facet_wrap(~season)


## ----resp_map------------------------------------------------------------------------------------------------------------------------------------------------------------------
#| fig-column: body-outset
#| out-width: 100%
ggplot(df_fit) + 
  geom_polygon(data = world, aes(x = lon, y = lat, group = group), fill = "grey") +
  geom_point(aes(x = lon, y = lat, colour = log_doc_surf), size = 0.5) +
  ggplot2::scale_colour_viridis_c(option = "A") +
  coord_quickmap(expand = 0) +
  facet_wrap(~season)


## ----exp_pca-------------------------------------------------------------------------------------------------------------------------------------------------------------------
#| fig-column: body-outset
#| out-width: 100%
# Need to remove lon and lat and to scale because units differ between variables
df_pca <- df_fit %>% select(c(lon, lat, season, doc_surf, log_doc_surf, everything())) %>% mutate(season = as.character(season))
pca_all <- FactoMineR::PCA(df_pca, quanti.sup = 1:2, quali.sup = 3, scale.unit = TRUE, graph = FALSE)

# Plot eigenvalues
plot_eig(pca_all)

# Plot variables
plot(pca_all, choix = "var", axes = c(1, 2), cex = 0.7)

## Get coordinates of individuals
inds <- pca_all$ind$coord %>% as_tibble()
# Set nice names for columns
colnames(inds) <- str_c("dim", paste(c(1:ncol(inds))))
# And join with initial dataframe of objects
inds <- df_fit %>% bind_cols(inds)

# Plot maps of PCA projections
ggmap(inds, "dim1", type = "point", palette = div_pal)
ggmap(inds, "dim2", type = "point", palette = div_pal)


## ----split---------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Split by season
df_fit_seas <- df_fit %>% group_by(season) %>% group_split()

# For each season, generate a stratified nested CV
folds <- lapply(1:length(df_fit_seas), function(i) {
  # Get data for this season
  df_seas <- df_fit_seas[[i]]
  
  set.seed(seed)
  # Stratified CV on deciles of response variable
  folds <- nested_cv(
    df_seas,
    outside = vfold_cv(v = 10, strata = log_doc_surf, breaks = 9),
    inside = vfold_cv(v = 10, strata = log_doc_surf, breaks = 9)) %>%
    mutate(cv_type = "stratified") %>% 
    mutate(season = as.character(i))
}) %>% 
  bind_rows()


## ----def_mod-------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Define a xgboost model with hyperparameters to tune
xgb_spec <- boost_tree(
  trees = tune(),
  tree_depth = tune(),
  min_n = tune(),
  learn_rate = tune()
) %>%
  set_mode("regression") %>%
  set_engine("xgboost")


## ----def_form------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Generate formula from list of explanatory variables
xgb_form <- as.formula(paste("log_doc_surf ~ ", paste(c("doc_surf", exp_vars), collapse = " + "), sep = ""))


## ----def_grid------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Define one grid for all folds
set.seed(seed)
xgb_grid <- grid_latin_hypercube(
  trees(),
  learn_rate(),
  tree_depth(),
  min_n(),
  size = 30
)


## ----gridsearch----------------------------------------------------------------------------------------------------------------------------------------------------------------
message("Gridsearch")
#|cache.lazy: false
res <- pbmclapply(1:nrow(folds), function(i){
  
  ## Get fold
  x <- folds[i,]
  
  message(paste0("Processing ", x$id, " of season ", x$season))

  ## Train and test sets
  df_train <- analysis(x$splits[[1]]) %>% as_tibble()
  df_test <- assessment(x$splits[[1]]) %>% as_tibble()

  ## Recipe
  xgb_rec <- recipe(xgb_form, data = df_train) %>%
    update_role(doc_surf, new_role = "untransformed outcome")

  ## Workflow
  xgb_wflow <- workflow() %>%
    add_recipe(xgb_rec) %>%
    add_model(xgb_spec)

  ## Gridsearch
  set.seed(seed)
  xgb_res <- tune_grid(
    xgb_wflow,
    resamples = x$inner_resamples[[1]],
    grid = xgb_grid,
    metrics = metric_set(rmse),
    control = control_grid(save_pred = TRUE)
  )
  best_params <- select_best(xgb_res)

  ## Final fit
  final_xgb <- finalize_workflow(
    xgb_wflow,
    best_params
  )
  final_res <- fit(final_xgb, df_train)

  ## Prediction on outer folds
  preds <- predict(final_res, new_data = df_test) %>%
    bind_cols(df_test %>% select(log_doc_surf))

  ## Model explainer
  # Select only predictors
  vip_train <- xgb_rec %>% prep() %>% bake(new_data = NULL, all_predictors())

  # Explainer
  xgb_explain <- explain_tidymodels(
      model = extract_fit_parsnip(final_res),
      data = vip_train,
      y = df_train %>%  pull(log_doc_surf),
      verbose = FALSE
    )

  # Variable importance
  full_vip <- model_parts(xgb_explain) %>%
    bind_rows() %>%
    filter(variable != "_baseline_")

  # CP profiles for all variables
  #selected_points <- ingredients::select_sample(df_train, n = 100, seed = seed)
  #cp_profiles <- ingredients::ceteris_paribus(xgb_explain, selected_points) %>% as_tibble()
  # CP profiles
  cp_profiles <- lapply(exp_vars, function(my_var){
    model_profile(explainer = xgb_explain, variables = my_var)$cp_profiles %>% as_tibble()
  }) %>%
    bind_rows()
  
  ## Predict new data
  new_preds <- augment(final_res, new_data = df_pred %>% mutate(doc_surf = NA)) %>% 
    select(-doc_surf) %>% 
    rename(pred_doc_log = .pred) %>% 
    relocate(pred_doc_log, .after = lat)
  
  ## Return results
  return(tibble(
      resp = resp_var[2],
      season = x$season,
      cv_type = x$cv_type,
      fold = x$id,
      preds = list(preds),
      importance = list(full_vip),
      cp_profiles = list(cp_profiles),
      new_preds = list(new_preds)
    ))
}, mc.cores = 10) %>%
  bind_rows()


## ----save----------------------------------------------------------------------------------------------------------------------------------------------------------------------
#|cache.lazy: false
save(res, file = output_filename)

