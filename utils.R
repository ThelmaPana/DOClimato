## Packages ----
#--------------------------------------------------------------------------#
library(tidyverse)
library(parallel)
library(here)
library(pbmcapply)

## Reading & downloading
library(ncdf4)
library(R.matlab)
library(readxl)

## Processing
library(castr)
library(fields)
library(oce)
library(vegan)
library(plotpca)
library(Hmisc)
library(sf)
library(spatialsample)
library(abind)
library(marmap)

## Plots
library(cmocean)
library(chroma)
library(rnaturalearth)

## Modeling
library(tidymodels)
library(DALEXtra)
library(bonsai)



## Default values ----
#--------------------------------------------------------------------------#
# Seed for reproducibility
seed <- 11
set.seed(seed)


# Number of cores for parallel computing
n_cores <- 10

# GGplot theme
theme_set(theme_minimal())

# Path to data folder
data_dir <- "data"

## WOA
# Path to local WOA links
woa_loc <- here(file.path(data_dir, "raw/woa"))
dir.create(woa_loc, showWarnings=FALSE, recursive=TRUE)

# Path to WOA dataset
woa_dir <- "~/Documents/work/datasets/woa"


# Whether to perform download for WOA data
download_woa <- FALSE

# Whether to perform download for bathymetry data
download_bathy <- FALSE

# Max depth to compute clines from WOA
max_depth_woa <- 500 # compromise between coverage and computation cost

# Max depth of the layers we consider for predictors
surface_bottom <- 10
meso_top <- 200
meso_bottom <- 1000
# bathy is anything below meso

## UVP
# Depth above which to keep UVP objects
max_depth_uvp <- 1000

# Minimum number of objects in a UVP profile to consider it
n_min_uvp <- 10

# Number of best predictors to plot partial dependence plots for
n_pdp <- 3



## Colour palettes ----
#--------------------------------------------------------------------------#

## Colour palettes for image.plot
col_temp   <- cmocean("thermal")(100)
col_sal    <- cmocean("haline")(100)
col_dens   <- cmocean("dense")(100)
col_oxy    <- brewer_colors(100, "Blues")
col_aou    <- brewer_colors(100, "Purples")
col_nit    <- cmocean("tempo")(100)
col_phos   <- brewer_colors(100, "BuPu")
col_sil    <- brewer_colors(100, "PuBu")
col_chl    <- cmocean("algae")(100)
col_irr    <- cmocean("solar")(100)
col_depth  <- cmocean("deep")(100)
col_alk    <- brewer_colors(100, "RdPu")
col_bbp    <- cmocean("turbid")(100)
col_npp    <- brewer_colors(100, "YlGn")
col_psd    <- brewer_colors(100, "Greens")
col_fmicro <- brewer_colors(100, "Greys")
col_dic    <- brewer_colors(100, "YlOrBr")
col_pic    <- brewer_colors(100, "OrRd")
col_poc    <- cmocean("matter")(100)
col_misc   <- viridis_colors(100)



## Colour palettes for ggplot
div_pal <- scale_colour_gradient2(low = "#4575b4", mid = "#ffffbf", high = "#d73027") # diverging palette centered at 0


## Get world map data ----
#--------------------------------------------------------------------------#
# For plots
world <- fortify(map_data('world', wrap = c(-180, 180))) %>% rename(lon = long)
world_sf <- ne_countries(scale = "medium", type = "countries", returnclass = "sf")

# To compute inland points
coast <- read_csv(here("data/raw/gshhg_world_c.csv"), col_types = cols())


## Function to increment file numbers in batch ----
#--------------------------------------------------------------------------#
increment_files <- function(to_increment){
  ## Scripts
  # List files with this number
  to_rename <- list.files(pattern = paste0(to_increment, ".*"))

  # Compute new number and generate new file names
  new_nb <- as.character(as.numeric(to_increment) + 1) %>% str_pad(width = 2, pad = "0")
  renamed <- str_replace_all(to_rename, pattern = to_increment, replacement = new_nb)

  # Rename files
  file.rename(from = to_rename, to = renamed)

  ## Data
  # List files with this number
  to_rename <- list.files(path = "data", pattern = paste0(to_increment, ".*"))

  # Compute new number and generate new file names
  new_nb <- as.character(as.numeric(to_increment) + 1) %>% str_pad(width = 2, pad = "0")
  renamed <- str_replace_all(to_rename, pattern = to_increment, replacement = new_nb)

  # Rename files
  file.rename(from = to_rename, to = renamed)
}


## Compute percentage of NA in a vector ----
#--------------------------------------------------------------------------#
#' Compute percentage of NA in a vector
#' @param x input vector
#' @return number in [0,1] = the percentage of NA; returns 1 is x is empty.
percent_na <- function(x) {
  if (length(x) == 0) {
    res <- 1
  } else {
    res <- sum(is.na(x)) / length(x)
  }
  return(res)
}


## Read mat files from Cael ----
#--------------------------------------------------------------------------#
#' Read mat file and extract variable of interest from Cael climatology
#'
#' @param file path to .mat file
#' @param var name of variable to extract
#' @param yearly whether to compute yearly average, if `FALSE` monthly data is returned
#'
#' @return an array
#'
#' @examples read_clim(file = "data/raw/clim4cael.mat", var = "cafes")
read_clim <- function(file, var, yearly = TRUE){
  # Read file
  clim <- readMat(file)

  # Extract variable of interest
  data <- clim[[var]]

  # Swap lon and lat dimensions
  data <- aperm(data, c(2,1,3))
  # Reorder lat
  data <- data[, ncol(data):1,]

  # Compute yearly from monthly
  if (yearly){
    data <- apply(data, c(1,2), mean, na.rm = TRUE)
  }
  return(data)
}

## Convert monthly data to seasonal data ----
#--------------------------------------------------------------------------#
#' Convert monthly data to seasonal data
#'
#' Given a dataset that is a 3D array (with 3rd dimensions as month) of observations in lon and lat,
#' returns a list of 4 (i.e. seasonal) of observations in lon and lat. 
#' For each season in "[1:4]", months can be computed as "(seas * 3 - 2):(seas * 3)"
#'  eg: season 1 is months 1 to 3, season 2 is months 4 to 6…
#'
#' @param data a 3D array, with 3rd dimension as month
#' @param to_array whether to return an array
#'
#' @return a list of 4 OR a 3D array
#'
#' @examples npp_seas <- month_to_seas(npp_month)
month_to_seas <- function(data, to_array = FALSE){
  # List seasons
  seas <- c(1, 2, 3, 4)
  
  # For each season, compute mean across appropriate months
  res <- mclapply(seas, function(s) {
    apply(data[,,(s * 3 - 2):(s * 3)], c(1,2), mean, na.rm = TRUE)
  }, mc.cores = 4)
  
  if (to_array) {
    res <- do.call(abind, list(res, along = 3))
  }
  
  return(res)
}


## Plot a nice ggplot map ----
#--------------------------------------------------------------------------#
#' Plot a map, whether raster or points
#'
#' Generate a ggplot raster map from a dataframe and a variable name.
#' Palette can be inferred from the name of the variable or provided in the arguments.
#' If no palette is found, it defaults to viridis.
#' Land is plotted by default but this can be changed in arguments.
#'
#' @param df a dataframe, must contain at least 3 columns: `lon`, `lat` and var to plot
#' @param var a character, name of variable to plot
#' @param type a character, defining type of map (raster or points)
#' @param land a boolean, whether to plot land or not
#' @param palette a filling palette, will be generated if none is
#'
#' @return a ggplot object
#'
#' @examples ggmap(df, var = "temperature", type = "raster")
ggmap <- function(df, var, type = c("raster", "point"), land = TRUE, palette = NULL) {
  ## Check args
  # df is a dataframe containing "lon", "lat" and var
  checkmate::assert_data_frame(df)
  checkmate::assert_names(names(df), must.include = c("lon", "lat", var))
  # var is a character
  checkmate::assert_character(var)
  # type is either "raster" or "point"
  checkmate::assert_choice(type, c("raster", "point"))
  # land is a boolean
  checkmate::assert_logical(land)
  # palette is a palette or NULL
  checkmate::assert_multi_class(palette, c("ScaleContinuous", "Scale", "NULL"))

  ## Palettes
  # To look for palette, remove "_mean" in case we are working with annual climatology values
  var_pal <- str_remove(var, "_mean")
  # Remove also anything after "."
  var_pal <- str_split(var_pal, "\\.")[[1]][1]

  # If no palette is supplied, generate one
  if(is.null(palette)){

    # List of palettes for common variables
    pals <- tribble(
      ~vars, ~raster, ~point,
      "temperature", scale_fill_cmocean(name = "thermal", na.value = NA),                    scale_colour_cmocean(name = "thermal", na.value = NA),
      "salinity",    scale_fill_cmocean(name = "haline", na.value = NA),                     scale_colour_cmocean(name = "haline", na.value = NA),
      "density",     scale_fill_cmocean(name = "dense", na.value = NA),                      scale_colour_cmocean(name = "dense", na.value = NA),
      "oxygen",      scale_fill_distiller(palette = "Blues", na.value = NA, direction = 1),  scale_colour_distiller(palette = "Blues", na.value = NA, direction = 1),
      "nitrate",     scale_fill_cmocean(name = "tempo", na.value = NA),                      scale_colour_cmocean(name = "tempo", na.value = NA),
      "phosphate",   scale_fill_distiller(palette = "BuPu", na.value = NA, direction = 1),   scale_colour_distiller(palette = "BuPu", na.value = NA, direction = 1),
      "silicate",    scale_fill_distiller(palette = "PuBu", na.value = NA, direction = 1),   scale_colour_distiller(palette = "PuBu", na.value = NA, direction = 1),
      "chl",         scale_fill_cmocean(name = "algae", na.value = NA),                      scale_colour_cmocean(name = "algae", na.value = NA),
      "par",         scale_fill_cmocean(name = "solar", na.value = NA),                      scale_colour_cmocean(name = "solar", na.value = NA),
      "mld",         scale_fill_cmocean(name = "deep", na.value = NA),                       scale_colour_cmocean(name = "deep", na.value = NA),
      "mld_argo",    scale_fill_cmocean(name = "deep", na.value = NA),                       scale_colour_cmocean(name = "deep", na.value = NA),
      "thermo",      scale_fill_cmocean(name = "deep", na.value = NA),                       scale_colour_cmocean(name = "deep", na.value = NA),
      "pycno",       scale_fill_cmocean(name = "deep", na.value = NA),                       scale_colour_cmocean(name = "deep", na.value = NA),
      "z_eu",        scale_fill_cmocean(name = "deep", na.value = NA),                       scale_colour_cmocean(name = "deep", na.value = NA),
      "s_cline",     scale_fill_cmocean(name = "deep", na.value = NA),                       scale_colour_cmocean(name = "deep", na.value = NA),
      "p_cline",     scale_fill_cmocean(name = "deep", na.value = NA),                       scale_colour_cmocean(name = "deep", na.value = NA),
      "n_cline",     scale_fill_cmocean(name = "deep", na.value = NA),                       scale_colour_cmocean(name = "deep", na.value = NA),
      "alkalinity",  scale_fill_distiller(palette = "PuRd", na.value = NA, direction = 1),   scale_colour_distiller(palette = "PuRd", na.value = NA, direction = 1),
      "bbp",         scale_fill_cmocean(name = "turbid", na.value = NA),                     scale_colour_cmocean(name = "turbid", na.value = NA),
      "npp",         scale_fill_distiller(palette = "YlGn", na.value = NA, direction = 1),   scale_colour_distiller(palette = "YlGn", na.value = NA, direction = 1),
      "psd",         scale_fill_distiller(palette = "Greens", na.value = NA, direction = 1), scale_colour_distiller(palette = "Greens", na.value = NA, direction = 1),
      "fmicro",      scale_fill_distiller(palette = "Greys", na.value = NA, direction = 1),  scale_colour_distiller(palette = "Greys", na.value = NA, direction = 1),
      "dic",         scale_fill_distiller(palette = "YlOrBr", na.value = NA, direction = 1), scale_colour_distiller(palette = "YlOrBr", na.value = NA, direction = 1),
      "pic",         scale_fill_distiller(palette = "OrRd", na.value = NA, direction = 1),   scale_colour_distiller(palette = "OrRd", na.value = NA, direction = 1),
      "poc",         scale_fill_cmocean(name = "matter", na.value = NA),                     scale_colour_cmocean(name = "matter", na.value = NA)
    )

    # Find the matching palette for variable to plot
    pal <- pals %>% filter(vars == var_pal) %>% pull(all_of(type))

    # If no palette is found, use viridis
    if (length(pal) == 0 & type == "raster"){
      pal <- scale_fill_viridis_c(na.value = NA)
    } else if (length(pal) == 0 & type == "point"){
      pal <- scale_colour_viridis_c(na.value = NA)
    }

  } else { # If palette is supplied, use it!
    pal <- palette
  }

  ## Plot
  # Base plot
  p <- ggplot(df) +
    coord_quickmap(expand = FALSE) +
    theme_minimal()

  # add raster or point layer
  if (type == "raster"){
    p <- p + geom_raster(aes(x = lon, y = lat, fill = .data[[var]]), na.rm = TRUE)
  } else {
    p <- p + geom_point(aes(x = lon, y = lat, colour = .data[[var]]), na.rm = TRUE, size = 0.5)
  }

  # add land if required
  if (land){p <- p + geom_polygon(data = world, aes(x = lon, y = lat, group = group), fill = "gray")}

  # add palette
  p <- p + pal

  ## Return plot
  return(p)
}


## Preprocess images for morphr ----
#--------------------------------------------------------------------------#
# Function to pre-process images
preprocess <- function(x) {
  x |>
    # remove 31 pixels from the bottom (=the scale bar)
    img_chop(bottom=31) |>
    # change the gamma to see the light objects better
    img_adjust_gamma(gamma=0.7)
}

## Draw a circle ----
#--------------------------------------------------------------------------#
# To draw a circle
circleFun <- function(center = c(0,0),diameter = 1, npoints = 100){
  r = diameter / 2
  tt <- seq(0,2*pi,length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}


## Plot variable importance ----
#--------------------------------------------------------------------------#
ggplot_imp <- function(...) {
  obj <- list(...)
  metric_name <- attr(obj[[1]], "loss_name")
  metric_lab <- paste(metric_name,
                      "after permutations\n(higher indicates more important)")

  full_vip <- bind_rows(obj) %>%
    filter(variable != "_baseline_")

  perm_vals <- full_vip %>%
    filter(variable == "_full_model_") %>%
    group_by(label) %>%
    summarise(dropout_loss = mean(dropout_loss))

  p <- full_vip %>%
    filter(variable != "_full_model_") %>%
    mutate(variable = forcats::fct_reorder(variable, dropout_loss)) %>%
    ggplot(aes(dropout_loss, variable))
  if(length(obj) > 1) {
    p <- p +
      facet_wrap(vars(label)) +
      geom_vline(data = perm_vals, aes(xintercept = dropout_loss, color = label),
                 linewidth = 1.4, lty = 2, alpha = 0.7) +
      geom_boxplot(aes(color = label, fill = label), alpha = 0.2)
  } else {
    p <- p +
      geom_vline(data = perm_vals, aes(xintercept = dropout_loss),
                 linewidth = 1.4, lty = 2, alpha = 0.7) +
      geom_boxplot(fill = "#91CBD765", alpha = 0.4)

  }
  p +
    theme(legend.position = "none") +
    labs(x = metric_lab,
         y = NULL,  fill = NULL,  color = NULL)
}


## Plot partial dependence plot ----
#--------------------------------------------------------------------------#
# Generate plot from explainer and variable
plot_pdp <- function(explainer, variable, unlog = FALSE){
  xgb_pdp_1 <- model_profile(explainer = xgb_explain, variables = variable)
  ggplot_pdp(xgb_pdp_1, variable, unlog = unlog) + labs(x = variable, y = "Logged predicted POC")
}

# Plot pdp from model profile
ggplot_pdp <- function(obj, x, unlog = FALSE) {

  # Get cp profiles
  df <- as_tibble(obj$cp_profiles)

  # Unlog predictions
  if (unlog){
    df <- df %>% mutate(`_yhat_` = exp(`_yhat_`))
  }

  # Compute mean across cp profiles for each x value
  # Compute sd across centered cp profiles
  df <- df %>%
    group_by(`_ids_`) %>%
    mutate(yhat_cent = `_yhat_` - mean(`_yhat_`)) %>% # center cp profiles
    ungroup() %>%
    group_by(across(all_of(x))) %>%
    summarise(
      yhat_loc = mean(`_yhat_`), # compute mean of profiles
      yhat_spr = sd(yhat_cent) # compute sd of cp profiles
    ) %>%
    ungroup() %>%
    arrange(all_of(x)) %>%
    setNames(c("x", "yhat_loc", "yhat_spr"))


  # Plot
  p <- ggplot(df) +
    geom_ribbon(aes(x = x, ymin = yhat_loc - yhat_spr, ymax = yhat_loc + yhat_spr), fill = "gray80", alpha = 0.5) +
    geom_path(aes(x = x, y = yhat_loc), color = "midnightblue", linewidth = 1.2)

  p
}


## Average pdp across CV folds ----
#--------------------------------------------------------------------------#
#' Average pdp across CV folds.
#' 
#' Steps as follows for each cv_type and each variable
#' 1- compute the mean and spread of cp profiles within each fold
#' 2- interpolate yhat value and spread within each fold using a common set of x values
#' 3- perform a weighted average of yhat value and spread, using 1/var as weights
#'
#' @param cp dataframe of cp profiles
#' @param vars dataframe with variables to use
#'
#' @return
#'
prep_pdp <- function(cp, vars){
  # Get names of folds, for later use
  folds <- sort(unique(cp$fold))
  
  # Apply on each cv_type and variable
  mean_pdp <- lapply(1:nrow(vars), function(r){
    
    # Get variable and cvtype
    var_name <- vars[r,]$variable
    cv_type_name <- vars[r,]$cv_type
    
    ## Get corresponding CP profiles, compute mean and spread for each fold (step 1)
    d_pdp <- cp %>% 
      filter(cv_type == cv_type_name & `_vname_` == var_name) %>% 
      select(cv_type, fold, `_yhat_`, `_vname_`, `_ids_`, all_of(var_name)) %>% 
      arrange(`_ids_`, across(all_of(var_name))) %>% 
      # center each cp profiles across fold, variable and ids
      group_by(cv_type, fold, `_vname_`, `_ids_`) %>%
      mutate(yhat_cent = `_yhat_` - mean(`_yhat_`)) %>% # center cp profiles
      ungroup() %>%
      # compute mean and sd of centered cp profiles for each fold and value of the variable of interest
      group_by(cv_type, fold, across(all_of(var_name))) %>%
      summarise(
        yhat_loc = mean(`_yhat_`), # compute mean of profiles
        yhat_spr = sd(yhat_cent), # compute sd of cp profiles
        .groups = "keep"
      ) %>%
      ungroup() %>% 
      setNames(c("cv_type", "fold", "x", "yhat_loc", "yhat_spr"))
    
    ## Interpolate yhat values and spread on a common x distribution (step 2)
    # Regularise across folds: need a common x distribution, and interpolate y on this new x
    new_x <- quantile(d_pdp$x, probs = seq(0, 1, 0.01), names = FALSE)
    # x is different within each fold, so interpolation is performed on each fold
    
    int_pdp <- lapply(1:length(folds), function(i){
      # Get data corresponding to this fold
      fold_name <- folds[i]
      this_fold <- d_pdp %>% filter(fold == fold_name)
      
      # Extract original x values
      x <- this_fold$x
      # Extract values to interpolate (yhat_loc and yhat_spr)
      yhat_loc <- this_fold$yhat_loc
      yhat_spr <- this_fold$yhat_spr
      # Interpolate yhat_loc and yhat_spr on new x values
      int <- tibble(
        x = new_x,
        yhat_loc = castr::interpolate(x = x, y = yhat_loc, xout = new_x),
        yhat_spr = castr::interpolate(x = x, y = yhat_spr, xout = new_x),
      ) %>% 
        mutate(
          cv_type = cv_type_name,
          fold = fold_name,
          var_name = var_name,
          .before = x
        )
      # Return the result
      return(int)
      
    }) %>% 
      bind_rows()
    
    ## Across fold, compute the weighted mean, using 1/var as weights (step 3)
    mean_pdp <- int_pdp %>% 
      group_by(cv_type, var_name, x) %>% 
      summarise(
        yhat_loc = wtd.mean(yhat_loc, weights = 1/(yhat_spr)^2),
        yhat_spr = wtd.mean(yhat_spr, weights = 1/(yhat_spr)^2),
        .groups = "drop"
      ) %>% 
      arrange(x)
    
    # Return the result
    return(mean_pdp)
  }) %>% 
    bind_rows()
  
  return(mean_pdp)
}

