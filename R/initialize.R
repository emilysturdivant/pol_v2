# v2: using GBIF points that I downloaded
# Input: points downloaded from GBIF in folder input_data/GBIF/family_order_query
# Output: files in folder for current iteration of RF: data_out/sdm/rfX_params/query_term

# Load libraries ----
library('patchwork')
library('terra')
library('dismo')
library('raster') # use raster for compatibility with dismo
library('stars')
library('sf')
library('tidyverse')
library('purrr')

# Parameters ----
unq_cells = TRUE
excludep = TRUE
rf_vers <- 4
contin_vars_only <- TRUE
contin_vars_and_lc <- TRUE
vars <- c('lc', 'alt', 'wc')
filt_dates = FALSE
date_range <- c(2000, 2020)

# Function to create pred_dir path
get_pred_dir <- function(unq_cells = TRUE,
                         excludep = TRUE,
                         rf_vers = 1,
                         vars = c('biom', 'lc', 'alt', 'wc'),
                         filt_dates = FALSE) {
  
  exclude_code <- ifelse(excludep, 'excl', 'ovrlp')
  unq_code <- ifelse(unq_cells, 'unqc', 'unqpt')
  dfilt_code <- ifelse(filt_dates, '2000to2020', 'alldates')
  var_code <- str_c(c( str_to_title(vars)), collapse = '')
  
  rf_name <- str_c('rf', str_c(rf_vers, exclude_code, unq_code, dfilt_code, var_code, sep='_'))
  fp_tail <- file.path('sdm', rf_name)
  file.path('data', 'data_out', fp_tail)
}

# directory paths ----
var_code <- str_c(c( str_to_title(vars)), collapse = '')
(pred_dir <- get_pred_dir(unq_cells = unq_cells,
                         excludep = excludep,
                         rf_vers = rf_vers,
                         vars = vars,
                         filt_dates = filt_dates))
rf_fig_dir <- pred_dir
dir.create(pred_dir, recursive=T, showWarnings = F)

pts_dir <- 'data/tidy/pollinator_points'

# Environment variables
crop_dir <- file.path('data', 'input_data', 'environment_variables', 'cropped')

# Mexico for masking and mapping
mex <- raster::getData('GADM', country='MEX', level=1, 
                       path='data/input_data/context_Mexico') %>% 
  st_as_sf()
mex0 <- st_union(mex)

# Set extent for testing
ext <- raster::extent(mex)

# Functions ----

prep_predictor_stack <- function(pred_grd, crop_dir, vars, mex0, overwrite = FALSE) {
  
  # Laad raster and return if it already exists
  if(file.exists(pred_grd) & !overwrite) {
    return(raster::stack(pred_grd))
  }
  
  # List variables for removal
  drop_lst <- c('biomes_CVEECON2', 'biomes_CVEECON1', 'biomes_CVEECON4',
                'ESACCI.LC.L4.LC10.Map.10m.MEX_2016_2018', 
                'usv250s6gw_USV_SVI')
  
  # Include categorical layers for removal
  if(!'biom' %in% vars) drop_lst <- c(drop_lst, 'biomes_CVEECON3')
  if(!'lc' %in% vars)   drop_lst <- c(drop_lst, 'ESACCI.LC.L4.LCCS.Map.300m.P1Y.2015.v2.0.7cds')
  
  # List files to include
  fps <- list.files(crop_dir, 'tif$', full.names=T)
  drop_fps <- drop_lst %>% purrr::map_chr(~ str_subset(fps, .x))
  fps <- fps[!(fps %in% drop_fps)]
  
  # Use raster
  mex0r <- raster::rasterize(as_Spatial(mex0), raster::raster(fps[[1]]))
  pred <- raster::stack(fps) %>% raster::mask(mex0r)
  names(pred) <- names(pred) %>% 
    str_remove_all('wc2\\.1_30s_') %>% 
    str_replace('ESA.*', 'landcover') %>% 
    str_replace('MEX_msk_alt', 'elevation')
  
  # Save for next time
  pred %>% writeRaster(pred_grd, overwrite = overwrite)
  
  # Return
  return(pred)
}

#' @export
model_species_rf <- function(sp_df,
                             pred, 
                             pred_dir,
                             sp_name,
                             mutually_exclusive_pa=TRUE,
                             unq_cells=TRUE, 
                             return_obj = 'objects') {
  
  # File paths
  sp_nospc <- str_replace(sp_name, ' ', '_')
  model_fp <- file.path(pred_dir, 'models', str_c(sp_nospc, '.rds'))
  eval_fp <- file.path(pred_dir, 'model_evals', str_c(sp_nospc, '.csv'))
  erf_fp <- file.path(pred_dir, 'model_evals', str_c(sp_nospc, '.rds'))
  
  if (file.exists(model_fp)) {
    print(str_c(sp_name, ': Model already created.'))
    rf <- readRDS(model_fp)
    erf <- readRDS(erf_fp)
    if(return_obj == 'paths') return(list(rf=model_fp, erf=erf_fp))
    return(list(rf=rf, erf=erf))
  }
  
  # Presence points
  # Convert to coords DF
  sp1 <- tryCatch({
    mutate(sp_df,
           lon = unlist(map(sp_df$geometry, 1)), 
           lat = unlist(map(sp_df$geometry, 2)))
    },
    warning = function(w) {
      mutate(sp_df,
             lon = unlist(map(sp_df$geom, 1)), 
             lat = unlist(map(sp_df$geom, 2)))
    }
    ) %>%
    st_drop_geometry() %>% 
    dplyr::select(lon, lat)
  
  # Background points
  set.seed(10)
  backg <- dismo::randomPoints(mask = pred, n=1000, 
                               p = as_Spatial(sp_df), 
                               excludep = excludep) %>% 
    as_tibble() %>% 
    rename(lon = 'x', lat = 'y')
  
  # Split presence into training and testing 
  set.seed(0)
  group <- kfold(sp1, 5)
  train_1 <- sp1[group != 1, ]
  test_1 <- sp1[group == 1, ]
  
  # Split background into training and testing
  set.seed(0)
  group <- kfold(backg, 5)
  train_0 <- backg[group != 1, ]
  test_0 <- backg[group == 1, ]
  
  # Extract environmental data 
  # Training dataset 
  train <- bind_rows(train_1, train_0)
  envtrain1 <- raster::extract(pred, train, cellnumbers=T, df=T)
  
  pb_train <- c(rep(1, nrow(train_1)), rep(0, nrow(train_0)))
  envtrain1 <- data.frame( cbind(pa = pb_train, envtrain1) ) %>% 
    mutate(across(starts_with(c('biomes', 'ESA', 'usv')), as.factor)) %>% 
    dplyr::select(-ID)
  
  # Remove duplicated cells
  if(unq_cells) {
    envtrain1 <- envtrain1 %>% distinct()
  }
  
  envtrain <- envtrain1 %>% dplyr::select(-cells)
  
  # Exit if there are fewer than 25 true presence points in the training set.
  if(nrow(filter(envtrain, pa == 1)) < 25) {
    print(str_c(sp_name, ': Fewer than 25 true non-duplicated presence points in the training set'))
    return(list(rf = NA, erf = NA))
  }
  
  # Testing datasets - get predictors for test presence and background points
  testpres <- data.frame( raster::extract(pred, test_1) ) %>%
    mutate(across(starts_with(c('biomes', 'ESA', 'usv')), as.factor))
  
  testbackg <- data.frame( raster::extract(pred, test_0) ) %>%
    mutate(across(starts_with(c('biomes', 'ESA', 'usv')), as.factor))
  
  # Standardize factor levels
  vars <- envtrain %>% dplyr::select(where(is.factor)) %>% tbl_vars
  for(var in vars){
    levs <- tribble(~name, ~levels,
                    'train', levels(envtrain[[var]]),
                    'testb', levels(testbackg[[var]]), 
                    'testp', levels(testpres[[var]])) %>% 
      mutate(n =  purrr::map_dbl(levels, length))
    
    longest_levels <- levs %>% slice_max(n) %>% 
      dplyr::select(levels) %>% flatten %>% deframe
    
    levels(testpres[[var]]) <- longest_levels
    levels(testbackg[[var]]) <- longest_levels
    levels(envtrain[[var]]) <- longest_levels
  }
  
  # Random forest model
  suppressWarnings(
  rf <- randomForest::randomForest(
    pa ~ . -pa,
    data = envtrain,
    na.action = na.exclude,
    importance=T, 
    ntree=1000)
  )
  
  # Save
  dir.create(dirname(model_fp), recursive = T, showWarnings = F)
  saveRDS(rf, model_fp)
  # rf <- readRDS(model_fp)
  
  # Filenames
  dir.create(dirname(eval_fp), recursive = T, showWarnings = F)
  
  # Evaluate model with test data
  erf <- dismo::evaluate(testpres, testbackg, rf)
  
  # Save model evaluation
  saveRDS(erf, erf_fp)
  
  # Save simple statistics in CSV
  spc_eval <- tibble(
    species = sp_name,
    N_unq_pts = nrow(sp_df),
    N_unq_cells = nrow(filter(envtrain, pa == 1)),
    np = erf@np, 
    na = erf@na, 
    auc = erf@auc,
    cor = erf@cor, 
    pcor = erf@pcor, 
    thresh_spec_sens = dismo::threshold(erf, "spec_sens"), 
    thresh_kappa = dismo::threshold(erf, "kappa"), 
    thresh_no_omission = dismo::threshold(erf, "no_omission"), 
    thresh_prevalence = dismo::threshold(erf, "prevalence"))
  
  spc_eval %>% write_csv(eval_fp)
  
  if(!is.null(pred_dir)){
    # Get plot directory
    plot_fp <- file.path(pred_dir, 'var_importance', str_c(sp_nospc, '.png'))
    dir.create(dirname(plot_fp), recursive = T, showWarnings = F)
    
    # Save plot
    png(plot_fp)
    randomForest::varImpPlot(rf, type=1, sort=F, 
                             main=sp_name,
                             pt.cex=1,
                             bg='black')
    dev.off()
  }
  
  # Return
  if(return_obj == 'paths') {
    return(list(rf=model_fp, erf=erf_fp))
  }
  
  return(list(rf=rf, erf=erf))
}

#' @export
predict_distribution_rf <- function(rf_list, 
                                    sp_name, 
                                    pred_dir, 
                                    ext, 
                                    predictors,
                                    write_binned = FALSE, 
                                    thresh = 'spec_sens') {
  
  rf <- rf_list$rf
  erf <- rf_list$erf
  
  # Filepaths  
  sp_nospc <- str_replace(sp_name, ' ', '_')
  likelihood_fp <- file.path(pred_dir, 'likelihood', str_glue(sp_nospc, '.tif'))
  dir.create(dirname(likelihood_fp), recursive = TRUE, showWarnings = FALSE)
  binned_fp <- file.path(pred_dir, str_c('binned_', thresh), 
                         str_glue(sp_nospc, '.tif'))
  
  if(file.exists(likelihood_fp) & (file.exists(binned_fp) | !write_binned)) {
    if(!write_binned) {
      print(str_c(likelihood_fp, ' already exists.'))
    } else {
      print(str_c(likelihood_fp, ' and ', binned_fp, ' already exist.'))
    }
    return(likelihood_fp)
  } 
  
  if( !file.exists(likelihood_fp) ) {
    
    # Load model
    if(is.character(rf)){
      rf <- readRDS(rf)
    } else if (is.na(rf)) {
      print(str_c(sp_name, ': NA provided for model.'))
      return(NA)
    }
    
    # Create map and interpolate to fill holes
    pr_rf1 <- dismo::predict(predictors, rf, ext=ext)
    pr_rf1 <- raster::focal(pr_rf1, 
                            w=matrix(1,nrow=3, ncol=3), 
                            fun=mean, 
                            NAonly=TRUE, 
                            na.rm=TRUE) 
    pr_rf1 <- round(pr_rf1, 3)
    
    # Save likelihood raster
    writeRaster(pr_rf1, likelihood_fp, overwrite=T, 
                options=c("dstnodata=-99999"), 
                wopt=list(gdal='COMPRESS=LZW',
                          datatype = ''))
    
  } else if(write_binned & !file.exists(binned_fp)) {
    # If likelihood tiff already exists, but now we want to get the binary...
    pr_rf1 <- raster(likelihood_fp)
    
  } else {
    # If likelihood tiff already exists and we don't need to create the binary...
    print(str_c(likelihood_fp, ' already exists and binary not requested.'))
    return(likelihood_fp)
  }
  
  # Presence-absence map
  dir.create(dirname(binned_fp), recursive = TRUE, showWarnings = FALSE)
  
  # Apply threshold from max TPR+TNR and save
  if(is.character(erf)) {
    erf <- readRDS(erf)
  } 
  
  tr <- threshold(erf, thresh)
  pa_rf1 <- pr_rf1 > tr
  
  # Save
  dir.create(dirname(binned_fp), recursive = T, showWarnings = F)
  writeRaster(pa_rf1, binned_fp, overwrite=T,
              options=c("dstnodata=-99999"), wopt=list(gdal='COMPRESS=LZW'))
  
  if(file.exists(likelihood_fp)) return(likelihood_fp)
  
  return(NA)
}

get_accuracy_metrics <- function(thresh, erf) {
  
  # Get threshold value
  threshold <- threshold(erf, thresh, sensitivity = 0.99)
  
  # Create table with predicted and observed presences and absences
  pa_tbl <- bind_rows(tibble(value = erf@presence, true_pres = 1),
                      tibble(value = erf@absence, true_pres = 0)) %>% 
    mutate(pred_pres = cut(value, 
                           breaks = c(-Inf, threshold, Inf), 
                           labels = c(0, 1)),
           result = case_when(true_pres == 1 & pred_pres == 1 ~ 'tp', 
                              true_pres == 1 & pred_pres == 0 ~ 'fn',
                              true_pres == 0 & pred_pres == 0 ~ 'tn',
                              true_pres == 0 & pred_pres == 1 ~ 'fp') %>% 
             factor(levels = c('tp', 'fp', 'fn', 'tn'))
    )
  
  # Calculate accuracy metrics
  cols <- c(tp = 0, fp = 0, fn = 0, tn = 0)
  results <- pa_tbl %>% 
    count(result) %>% 
    pivot_wider(names_from = result, values_from = n) %>% 
    add_column(!!!cols[!names(cols) %in% names(.)]) %>% 
    mutate(overall_acc = (tp + tn) / (tp + tn + fp + fn),
           sensitivity = tp / (tp + fn), # a.k.a. recall
           # ** more important for us because we're using pseudo-absences; 
           # What proportion of true positives was correctly identified? 
           # Inverse of omission error
           specificity = tn / (tn + fp),
           PPV = tp / (tp + fp), # aka. precision 
           # What proportion of positive predictions was correct?
           NPV = tn / (tn + fn), 
           PLR = sensitivity / (1 - specificity),
           NLR = (1 - sensitivity) / specificity, 
           TSS = sensitivity + specificity - 1, 
           odds_ratio = tp*tn / (fn*fp),
           ae = ((tp + fn) * (tp + fp) + (tn + fn) * (tn + fp)) / (tp + tn + fp + fn)^2,
           Kappa = (overall_acc - ae) / (1 - ae)) %>% 
    bind_cols(tibble(threshold = thresh), .)
}

#' @export
stack_sdms <- function(sp_fps, rich_tif_fp, rich_plot_fp, mex){
  
  pol_stack <- stack(sp_fps)
  
  # Sum layers and save 
  pol_rich <- sum(pol_stack, na.rm=T)
  pol_rich_msk <- mask(pol_rich, as_Spatial(mex))
  writeRaster(pol_rich_msk, rich_tif_fp, overwrite=T,
              options=c("dstnodata=-99999"), wopt=list(gdal='COMPRESS=LZW'))
  
  # Plot richness
  pol_rich_stars <- st_as_stars(pol_rich_msk)
  rich_plot  <- ggplot() +
    geom_stars(data=pol_rich_stars) +
    geom_sf(data = mex, fill = "transparent", size = 0.2, color = alpha("lightgray", 0.2)) +
    colormap::scale_fill_colormap(str_glue("Richness\n(N = {nlayers(pol_stack)})"), 
                                  na.value = "transparent", 
                                  colormap = colormap::colormaps$viridis) +
    ggthemes::theme_hc() +
    theme(legend.position=c(.95, 1), legend.title.align=0, legend.justification = c(1,1)) +
    labs(x = NULL, y = NULL)
  
  # Save
  ggsave(rich_plot_fp, rich_plot, width=9, height=5.7, dpi=120)
  
  # Facet individual PA maps
  # pa_facets  <- rasterVis::gplot(pol_stack) +
  #   geom_tile(aes(fill = value)) +
  #   colormap::scale_fill_colormap(str_glue("Occupancy likelihood"), na.value = "transparent",
  #                                 colormap = colormap::colormaps$viridis) +
  #   theme_minimal() +
  #   theme(legend.position='bottom') +
  #   labs(x = NULL, y = NULL) +
  #   coord_equal() +
  #   facet_wrap(~ variable, nrow=3)
  # 
  # # Save
  # ggsave(file.path(rf_fig_dir, str_glue('Likhd_{nlayers(pol_stack)}species.png')), pa_facets, width=9, height=5)
  
}

#' @export
stack_sdms_terra <- function(sp_fps, rich_tif_fp, rich_plot_fp, mex){
  
  interval <- 30
  idx <- seq(to = length(sp_fps), by = interval)
  
  temp_dir <- file.path(dirname(rich_tif_fp), 'temp')
  unlink(temp_dir, recursive = TRUE)
  dir.create(temp_dir, recursive = TRUE)
  
  for(i in idx) {
    fps_sub <- sp_fps[i:(i+interval-1)]
    fps_sub <- fps_sub[!is.na(fps_sub)]
    pol_stack <- terra::rast(fps_sub)
    
    temp_fp <- file.path(dirname(rich_tif_fp), 'temp', 
                         str_c(i, basename(rich_tif_fp)))
    # Sum layers and save 
    pol_rich <- sum(pol_stack, na.rm=T)
    pol_rich_msk <- terra::mask(pol_rich, terra::vect(mex), filename = temp_fp,
                                overwrite = TRUE)
  }

  # Sum temp layers
  temp_fps <- list.files(temp_dir, basename(rich_tif_fp))
  pol_stack <- terra::rast(fps_sub)
  pol_rich <- sum(pol_stack, na.rm=T)
  
  # Convert back to decimal and save
  pol_rich_msk <- terra::mask(pol_rich, terra::vect(mex), filename = rich_tif_fp,
                              overwrite = TRUE)
  
  # Plot richness
  pol_rich_stars <- st_as_stars(pol_rich_msk)
  rich_plot  <- ggplot() +
    geom_stars(data=pol_rich_stars) +
    geom_sf(data = mex, fill = "transparent", size = 0.2, color = alpha("lightgray", 0.2)) +
    colormap::scale_fill_colormap(str_glue("Richness\n(N = {length(sp_fps)})"), 
                                  na.value = "transparent", 
                                  colormap = colormap::colormaps$viridis) +    
    theme_minimal() +
    theme(legend.position = c(.95, 1), 
          legend.title.align = 0,
          legend.justification = c(1,1),
          plot.background = element_rect(fill = 'white'),
          panel.grid = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank())
  
  # Save
  ggsave(rich_plot_fp, rich_plot, width=9, height=5.7, dpi=120)
  
  # Facet individual PA maps
  # pa_facets  <- rasterVis::gplot(pol_stack) +
  #   geom_tile(aes(fill = value)) +
  #   colormap::scale_fill_colormap(str_glue("Occupancy likelihood"), na.value = "transparent",
  #                                 colormap = colormap::colormaps$viridis) +
  #   theme_minimal() +
  #   theme(legend.position='bottom') +
  #   labs(x = NULL, y = NULL) +
  #   coord_equal() +
  #   facet_wrap(~ variable, nrow=3)
  # 
  # # Save
  # ggsave(file.path(rf_fig_dir, str_glue('Likhd_{nlayers(pol_stack)}species.png')), pa_facets, width=9, height=5)
  
}

plot_lklhd_map <- function(lklhd_tif, filename = NULL){
  # Load as stars
  lklhd_stars <- read_stars(lklhd_tif)
  
  # Plot
  p <- ggplot() +
    geom_stars(data = lklhd_stars) +
    geom_sf(data = mex, fill = "transparent", size = 0.2, color = "gray70") +
    colormap::scale_fill_colormap("Occupancy\nlikelihood", na.value = "transparent", 
                                  colormap = colormap::colormaps$viridis) +
    theme_minimal() +
    theme(legend.position = c(.95, 1), 
          legend.title.align = 0,
          legend.justification = c(1,1),
          panel.grid.minor = element_blank(),
          panel.grid = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank())
  
  # Save
  if(is.character(filename)) ggsave(filename, p, width=9, height=5.7, dpi=120)
  
  # Return
  return(p)
}

plot_binned_map <- function(bin_tif, filename = NULL) {
  # Load as stars
  bnnd_stars <- read_stars(bin_tif) %>% 
    mutate(across(everything(), ~ as.factor(.x)), 
           across(everything(), ~ na_if(.x, 0)))
  
  # Plot
  p <- ggplot() +
    geom_stars(data = bnnd_stars) +
    geom_sf(data = mex, fill = "transparent", size = 0.2, color = "gray70") +
    scale_fill_manual(values = c('mediumblue'), na.value = "transparent", 
                      labels = c('likely present')) +
    theme_minimal() +
    theme(legend.position = c(.95, 1), 
          legend.title = element_blank(),
          legend.justification = c(1,1),
          panel.grid = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank())
  
  # Save
  if(is.character(filename)) ggsave(filename, p, width=9, height=5.7, dpi=120)
  
  # Return
  return(p)
}

plot_qc_maps <- function(sp_row, fig_dir) {
  
  # Get data for species
  sp_df <- sp_row$data[[1]]
  sp_ch <- sp_row$convhull[[1]]
  sp_name <- sp_row$species
  sp_nospc <- str_replace(sp_name, ' ', '_')
  group <- str_c(sp_row$common_group, sp_row$subgroup_a, sep = ', ')
  title <- bquote(.(group)~': '~italic(.(sp_name)))
  lklhd_fp <- file.path(pred_dir, 'likelihood', str_c(sp_nospc, '.tif'))
  bnnd_fp <- file.path(pred_dir, 'binned_spec_sens', str_c(sp_nospc, '.tif'))
  
  # Filepaths
  plot_fp <- file.path(fig_dir, str_c(sp_nospc, '_maps.png'))
  
  # Get likelihood raster
  # Create TIFF if it doesn't already exist
  if (!file.exists(lklhd_fp)) {
    # Filepaths 
    rf_fp <- file.path(pred_dir, 'models', str_c(sp_nospc, '.rds'))
    erf_fp <- file.path(pred_dir, 'model_evals', str_c(sp_nospc, '.rds'))
    
    # Make TIFs of likelihood and presence/absence
    predict_distribution_rf(rf_fp, erf_fp, sp_name, pred_dir, ext, pred, write_binned = TRUE)
    
  }
  
  # Plot
  likelihood_plot <- plot_lklhd_map(lklhd_fp)
  
  # Add presence points to likelihood plot
  like_plot <- likelihood_plot +
    scale_fill_gradient(low = 'gray50', high = 'white', na.value = "transparent") +
    geom_sf(data = sp_df, color = "firebrick3", size = .5, alpha = 0.5,
            show.legend = TRUE) +
    # scale_color_manual(name = "presence") +
    guides(fill = 'none')
  
  # Create P/A raster if doesn't exist
  if (!file.exists(bnnd_fp)) {
    # Filepaths 
    rf_fp <- file.path(pred_dir, 'models', str_c(sp_nospc, '.rds'))
    erf_fp <- file.path(pred_dir, 'model_evals', str_c(sp_nospc, '.rds'))
    
    # Make TIFs of likelihood and presence/absence
    predict_distribution_rf(rf_fp, erf_fp, sp_name, pred_dir, ext, pred, write_binned = TRUE)
    
  }
  
  # Plot
  binned_map <- plot_binned_map(bnnd_fp)
  
  # Get model evaluation values
  erf_csv <- file.path(pred_dir, 'model_evals', str_c(sp_nospc, '.csv'))
  erf_rds <- file.path(pred_dir, 'model_evals', str_c(sp_nospc, '.rds'))
  
  erf_tbl <- read_csv(erf_csv) %>% 
    select(where(is.numeric)) %>% 
    mutate(across(where(is.numeric), ~ signif(.x, digits = 4)))
  
  # Table with model statistics
  stat_grob <- erf_tbl %>% 
    select(-starts_with('thresh')) %>% 
    pivot_longer(everything(), 
                 names_to = 'Statistic', 
                 values_to = 'Value') %>% 
    column_to_rownames("Statistic") %>% 
    gridExtra::tableGrob()
  
  # Table with threshold values
  thresh_grob <- erf_tbl %>% 
    select(starts_with('thresh')) %>% 
    pivot_longer(everything(), 
                 names_to = 'Statistic', 
                 names_prefix = 'thresh_',
                 values_to = 'Threshold') %>% 
    column_to_rownames("Statistic") %>% 
    gridExtra::tableGrob()
  
  # Table with accuracy values after binning at given thresholds
  erf <- readRDS(erf_rds)
  acc_tbl <- get_accuracy_metrics('spec_sens', erf) %>% 
    mutate(across(where(is.numeric), ~ round(.x, digits = 2))) 
  acc_grob <- acc_tbl %>% 
    select(sensitivity, TSS, Kappa) %>% 
    pivot_longer(everything(), 
                 names_to = 'Statistic', 
                 values_to = 'Accuracy') %>% 
    column_to_rownames("Statistic") %>% 
    gridExtra::tableGrob() 
  
  # Layout maps
  layout <- '
  AAABBB
  AAABBB
  CDEGGG
  CDEGGG
  '
  layout <- '
  AAAABBBB
  AAAABBBB
  CCDDGGGG
  CCDDGGGG
  '
  maps <- wrap_plots(A = like_plot, 
                     B = likelihood_plot, 
                     C = stat_grob, D = thresh_grob, 
                     # E = acc_grob, 
                     G = binned_map + 
                       inset_element(acc_grob, left = 0.05, bottom = 0, 
                                     right = 0.3, top = 0.35
                                     ), 
                     design = layout) +
    plot_annotation(
      title = title,
      # tag_levels = 'a',
      # tag_suffix = ')', 
      caption = 'Top row shows relative likelihood of species presence. Red points are observations. 
      Bottom right are areas with likelihood greater than the spec_sens threshold. 
      Right table shows accuracy metrics for the binary map.'
    )
  
  # Save
  ggsave(plot_fp, maps, width=12, height=8, dpi=120)
}
