# v2: using GBIF points that I downloaded
# Input: points downloaded from GBIF in folder input_data/GBIF/family_order_query
# Output: files in folder for current iteration of RF: data_out/sdm/rfX_params/query_term

# Load libraries ----
# library(tools)
# library(maptools)
# library(mapview)
# library(units)
# library(rgbif)
# library(scrubr)
# library('terra')
# library(tmap)
# tmap::tmap_mode('view')
library('patchwork')
library('dismo')
library('raster') # use raster for compatibility with dismo
library('stars')
library('sf')
library('tidyverse')

# Initialize -----
unq_cells = TRUE
mutually_exclusive_pa = TRUE
nspecies <- 15
rf_vers <- 1
contin_vars_only <- TRUE

# Date range
filt_dates = FALSE
date_range <- c(2000, 2020)

pts_dir <- 'data/tidy/pollinator_points'
# noApis = TRUE
# apis_code <- ifelse(order == 'Hymenoptera', ifelse(noApis, '_noApis', ''), '')

# directory paths ----
unq_code <- ifelse(unq_cells, 'unq_cells', 'unq_pts')
unq_code <- ifelse(mutually_exclusive_pa, 'excl', unq_code)
dfilt_code <- ifelse(filt_dates, '2000to2020', 'alldates')
cont_var_code <- ifelse(contin_vars_only, '_contpred', '')

rf_name <- str_c('rf', str_c(rf_vers, unq_code, dfilt_code, sep='_'), cont_var_code)
fp_tail <- file.path('sdm', rf_name)
pred_dir <- file.path('data', 'data_out', fp_tail)
rf_fig_dir <- pred_dir
dir.create(pred_dir, recursive=T, showWarnings = F)

# Functions ----
#' @export
model_species_rf <- function(sp_df,
                             pred, 
                             pred_dir,
                             sp_name,
                             mutually_exclusive_pa=TRUE,
                             unq_cells=TRUE,
                             rf_fig_dir=NULL) {
  
  # File paths
  sp_nospc <- str_replace(sp_name, ' ', '_')
  model_fp <- file.path(pred_dir, 'models', str_c(sp_nospc, '.rds'))
  eval_fp <- file.path(pred_dir, 'model_evals', str_c(sp_nospc, '.csv'))
  erf_fp <- file.path(pred_dir, 'model_evals', str_c(sp_nospc, '.rds'))
  
  if (file.exists(model_fp)) {
    print(str_c(sp_name, ': Model already created.'))
    rf <- readRDS(model_fp)
    erf <- readRDS(erf_fp)
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
  pres_pts <- as_Spatial(sp_df)
  
  set.seed(10)
  if(mutually_exclusive_pa) {
    
    backg <- dismo::randomPoints(mask = pred, n=1000, p = pres_pts
                          # prob = T, # use mask as sampling bias grid
    )
    
    # bg <- spatSample(pred, 1000, "random", na.rm = TRUE, xy = TRUE)
    
  } else {
    
    backg <- dismo::randomPoints(pred, n=1000) 
  }
  
  backg <- backg %>% 
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
  
  # Testing datasets - get predictors for test presence and background points
  testpres <- data.frame( raster::extract(pred, test_1) ) %>%
    mutate(across(starts_with(c('biomes', 'ESA', 'usv')), as.factor))
  
  testbackg <- data.frame( raster::extract(pred, test_0) ) %>%
    mutate(across(starts_with(c('biomes', 'ESA', 'usv')), as.factor))
  
  # Set factor levels for test DFs to match training data
  vars <- envtrain %>% dplyr::select(where(is.factor)) %>% tbl_vars
  for(var in vars){
    levels(testpres[[var]]) <- levels(envtrain[[var]])
    levels(testbackg[[var]]) <- levels(envtrain[[var]])
  }
  
  # Random forest model
  rf <- randomForest::randomForest(
    pa ~ . -pa,
    data = envtrain,
    na.action = na.exclude,
    importance=T, 
    ntree=1000)
  
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
    species=sp_name,
    N_unq_pts=nrow(sp_df),
    N_unq_cells=nrow(filter(envtrain, pa == 1)),
    np=erf@np, 
    na=erf@na, 
    auc=erf@auc,
    cor=erf@cor, 
    pcor=erf@pcor, 
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
  return(list(rf=rf, erf=erf))
}

#' @export
predict_distribution_rf <- function(rf, 
                                    erf, 
                                    sp_name, 
                                    pred_dir, 
                                    ext, 
                                    predictors,
                                    write_binned = FALSE, 
                                    thresh = 'spec_sens') {
  
  # Filepaths  
  sp_nospc <- str_replace(sp_name, ' ', '_')
  likelihood_fp <- file.path(pred_dir, 'likelihood', str_glue(sp_nospc, '.tif'))
  dir.create(dirname(likelihood_fp), recursive = TRUE, showWarnings = FALSE)
  binned_fp <- file.path(pred_dir, str_c('binned_', thresh), 
                         str_glue(sp_nospc, '.tif'))
  
  if(file.exists(likelihood_fp) & (file.exists(binned_fp) | !write_binned)) {
    return()
  } 
  
  if( !file.exists(likelihood_fp) ) {
    
    # Load model
    if(is.character(rf)){
      rf <- readRDS(rf)
    } 
    
    # Create map and interpolate to fill holes
    pr_rf1 <- dismo::predict(predictors, rf, ext=ext)
    pr_rf1 <- raster::focal(pr_rf1, 
                            w=matrix(1,nrow=3, ncol=3), 
                            fun=mean, 
                            NAonly=TRUE, 
                            na.rm=TRUE) 
    
    # Save likelihood raster
    dir.create(dirname(likelihood_fp), recursive = T, showWarnings = F)
    writeRaster(pr_rf1, likelihood_fp, overwrite=T, 
                options=c("dstnodata=-99999"), wopt=list(gdal='COMPRESS=LZW'))
    
  } else if(write_binned & !file.exists(binned_fp)) {
    pr_rf1 <- raster(likelihood_fp)
    
  } else {
    return()
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

# Mexico
mex <- raster::getData('GADM', country='MEX', level=1, 
                       path='data/input_data/context_Mexico') %>% 
  st_as_sf()

# Load environment variables ----
crop_dir <- file.path('data', 'input_data', 'environment_variables', 'cropped')
predictors <- raster::stack(list.files(crop_dir, 'tif$', full.names=T))

# Remove layers from predictors
drop_lst <- c('biomes_CVEECON2', 'biomes_CVEECON1', 'biomes_CVEECON4',
              'ESACCI.LC.L4.LC10.Map.10m.MEX_2016_2018', 
              'usv250s6gw_USV_SVI')

if(contin_vars_only) {
  drop_lst <- c(drop_lst, 
                'biomes_CVEECON3', 'ESACCI.LC.L4.LCCS.Map.300m.P1Y.2015.v2.0.7cds')
}

pred <- predictors[[- which(names(predictors) %in% drop_lst) ]]

# Set extent for testing
ext <- raster::extent(mex)

# Load filtered points ----
filt_pts_rds <- file.path(pts_dir, str_c('points_nested_species_filt.rds'))
pol_df2 <- readRDS(filt_pts_rds)
pol_df2 <- pol_df2 %>% 
  filter(genus != "", species != "") %>% 
  filter(!is.na(genus), !is.na(species))

# # Get species list
# sp_counts <- pol_df2 %>%
#   select(species, genus, nobs) %>% 
#   arrange(desc(nobs))
# 
# if(noApis) {
#   sp_counts <- sp_counts %>%
#     filter(!str_detect(species, 'Apis mellifera'))
# }
# 
# sp_list <- sp_counts %>%
#   slice(1:nspecies) %>%
#   dplyr::select(species) %>%
#   deframe
# sp_nospc_list <- str_replace(sp_list, ' ', '_')

# ~ RF model for each species ----
df <- pol_df2 %>% slice(760:762)
stop <- nrow(df)
for (i in seq(1, stop)) {
  
  # Get data
  sp_dat <- df %>% slice(i)
  sp_name <- sp_dat$species
  sp_nospc <- sp_name %>% str_replace(' ', '_')
  print(str_glue("\n\n{i} of {stop}: {sp_name}"))
    
  # Get observation points
  sp_df <- sp_dat$data[[1]]
  
  # Save points
  pts_fp <- file.path(pred_dir, 'species_points', str_c(sp_nospc, '.gpkg'))
  if(!file.exists(pts_fp)) sp_df %>% st_write(pts_fp)
  
  # Run random forest and model evaluation
  mod_rf <- model_species_rf(sp_df,
                          pred, 
                          pred_dir,
                          sp_name,
                          mutually_exclusive_pa,
                          unq_cells)
  
  # Use model to create prediction rasters
  predict_distribution_rf(mod_rf$rf, mod_rf$erf, sp_name, pred_dir, ext, pred, write_binned = TRUE)
  
}

# # Add threshold values to model eval CSV ----
# expand_erf <- function(erf_fp) {
#   # File paths
#   sp_nospc <- tools::file_path_sans_ext(basename(erf_fp))
#   eval_fp <- file.path(pred_dir, 'model_evals', str_c(sp_nospc, '.csv'))
# 
#   # Load model evaluation
#   erf <- readRDS(erf_fp)
#   if(file.exists(eval_fp)){
#     spc_eval <- read_csv(eval_fp)
#   }
#   
#   # Save simple statistics in CSV
#   spc_eval <- bind_cols(
#     spc_eval, 
#     tibble(
#       thresh_kappa = dismo::threshold(erf, "kappa"), 
#       thresh_no_omission = dismo::threshold(erf, "no_omission"), 
#       thresh_prevalence = dismo::threshold(erf, "prevalence"))
#   )
#   
#   # Save
#   spc_eval %>% write_csv(eval_fp)
#   
# }
# 
# # Add threshold values to model eval CSV 
# fps <- list.files(file.path(pred_dir, 'model_evals'), '*.rds$', full.names = T)
# fps %>% walk(expand_erf)

# Look at model statistics together
fps <- list.files(file.path(pred_dir, 'model_evals'), '*.csv$', full.names = T)
if(length(fps) > 0) {
  eval_tbl <- fps %>% purrr::map_dfr(read.csv)
  eval_tbl_fp <- file.path(pred_dir, str_c('model_evals_', length(fps), 'species.csv'))
  eval_tbl %>% write_csv(eval_tbl_fp)
}

# Save TIFs of likelihood and presence/absence ----
# filter filepaths to species list
fps <- list.files(file.path(pred_dir, 'models'), '*.rds', full.names = T)
sp_fps <- sp_nospc_list %>% 
  map(~str_subset(.y, str_glue("{.x}.rds")), fps) %>% 
  flatten_chr()

for(rf_fp in sp_fps){
  
  # Filepaths  
  sp_nospc <- tools::file_path_sans_ext(basename(rf_fp))
  erf_fp <- file.path(pred_dir, 'model_evals', str_c(sp_nospc, '.rds'))
  sp_name <- sp_nospc %>% str_replace('_', '')
  
  # Make TIFs of likelihood and presence/absence
  predict_distribution_rf(rf_fp, erf_fp, sp_name, pred_dir, ext, pred)
  
}

# Look at model objects
rf_fp <- fps[[1]]
# sp_row <- pol_df2 %>% filter(species == 'Caria stillaticia')
erf_fp <- file.path(pred_dir, 'model_evals', str_c(sp_nospc, '.rds'))
rf_fp <- file.path(pred_dir, 'models', str_c(sp_nospc, '.rds'))
erf <- readRDS(erf_fp)
rf <- readRDS(rf_fp)

pa_tbl <- bind_rows(tibble(value = erf@presence, true_pres = 1),
                    tibble(value = erf@absence, true_pres = 0)) %>% 
  mutate(pred_pres = cut(value, 
                         breaks = c(-Inf, 0.0787, Inf), 
                         labels = c(0, 1)),
         result = case_when(true_pres == 1 & pred_pres == 1 ~ 'tp', 
                         true_pres == 1 & pred_pres == 0 ~ 'fn',
                         true_pres == 0 & pred_pres == 0 ~ 'tn',
                         true_pres == 0 & pred_pres == 1 ~ 'fp') %>% 
           factor(levels = c('tp', 'fp', 'fn', 'tn'))
         )
results <- pa_tbl %>% 
  count(result) %>% 
  pivot_wider(names_from = result, values_from = n) %>% 
  mutate(precision = tp / (tp + fp), # What proportion of positive predictions was correct?
         recall = tp / (tp + fn),# What proportion of true positives was correctly identified? ** more important for us because we're using pseudo-absences
         oa = (tp + tn) / (tp + tn + fp + fn),
         sens = tp / (tp + fn),
         spec = tn / (tn + fp),
         ppv = tp / (tp + fp),
         npv = tn / (tn + fn), 
         plr = sens / (1 - spec),
         nlr = (1 - sens) / spec, 
         tss = sens + spec - 1, 
         or = tp*tn / (fn*fp),
         ae = ((tp + fn) * (tp + fp) + (tn + fn) * (tn + fp)) / (tp + tn + fp + fn)^2,
         kappa = (oa - ae) / (1 - ae))

erf@TPR
erf@confusion
plot(erf, 'kappa')
boxplot(erf)
density(erf)

rf$predicted

# Create PNGs of all maps ----
(sp_row <- pol_df2 %>% sample_n(1))

# testing
for (sp_row in pol_df2) {
  
  # Get data for species
  sp_df <- sp_row$data[[1]]
  sp_ch <- sp_row$convhull[[1]]
  sp_name <- sp_row$species
  sp_nospc <- str_replace(sp_name, ' ', '_')
  group <- str_c(sp_row$common_group, sp_row$subgroup_a, sep = ', ')
  title <- str_c(group, ': ', sp_name)
  lklhd_fp <- file.path(pred_dir, 'likelihood', str_c(sp_nospc, '.tif'))
  bnnd_fp <- file.path(pred_dir, 'binned_spec_sens', str_c(sp_nospc, '.tif'))
  
  # Filepaths
  plot_fp <- file.path(rf_fig_dir, 'map_predictions', str_c(sp_nospc, '_maps.png'))
  dir.create(dirname(plot_fp), recursive = T, showWarnings = F)
  
  # Get likelihood raster
  # Create TIFF if it doesn't already exist
  if (!file.exists(lklhd_fp)) {
    # Filepaths 
    rf_fp <- file.path(pred_dir, 'models', str_c(sp_nospc, '.rds'))
    erf_fp <- file.path(pred_dir, 'model_evals', str_c(sp_nospc, '.rds'))
    
    # Make TIFs of likelihood and presence/absence
    predict_distribution_rf(rf_fp, erf_fp, sp_name, pred_dir, ext, pred, write_binned = TRUE)
    
  }
  
  # Load as stars
  lklhd_stars <- read_stars(lklhd_fp)
  
  # Plot
  likelihood_plot <- ggplot() +
    geom_stars(data = lklhd_stars) +
    geom_sf(data = mex, fill = "transparent", size = 0.2, color = "gray70") +
    colormap::scale_fill_colormap("Occupancy\nlikelihood", na.value = "transparent", 
                                  colormap = colormap::colormaps$viridis) +
    ggthemes::theme_hc() +
    theme(legend.position=c(.95, 1), legend.title.align=0, legend.justification = c(1,1)) +
    labs(x = NULL, y = NULL)
  
  like_plot <- likelihood_plot +
    # geom_sf(data = sp_ch, fill = "transparent", size = 0.2, color = "white") +
    geom_sf(data = sp_df, color = "firebrick3", size = 1)
  
  # Get P/A raster as stars object
  if (!file.exists(bnnd_fp)) {
    # Filepaths 
    rf_fp <- file.path(pred_dir, 'models', str_c(sp_nospc, '.rds'))
    erf_fp <- file.path(pred_dir, 'model_evals', str_c(sp_nospc, '.rds'))
    
    # Make TIFs of likelihood and presence/absence
    predict_distribution_rf(rf_fp, erf_fp, sp_name, pred_dir, ext, pred, write_binned = TRUE)
    
  }
  
  # Load as stars
  bnnd_stars <- read_stars(bnnd_fp) %>% 
    mutate(across(everything(), ~ as.factor(.x)), 
           across(everything(), ~ na_if(.x, 0)))

  # Plot
  binned_map <- ggplot() +
    geom_stars(data = bnnd_stars) +
    geom_sf(data = mex, fill = "transparent", size = 0.2, color = "gray70") +
    # scale_fill_manual("likely present", values = c('white', 'blue4'), na.value = "transparent") +
    scale_fill_manual(values = c('mediumblue'), na.value = "transparent", 
                      labels = c('likely present')) +
    ggthemes::theme_hc() +
    theme(legend.position=c(.95, 1), 
          # legend.title.align=0, 
          legend.title = element_blank(),
          legend.justification = c(1,1))+
    labs(x = NULL, y = NULL) 
  
  # Get model evaluation values
  erf_fp <- file.path(pred_dir, 'model_evals', str_c(sp_nospc, '.csv'))
  erf <- read_csv(erf_fp)
  grob1 <- erf %>% 
    select(where(is.numeric)) %>% 
    mutate(across(where(is.numeric), ~ signif(.x, digits = 4))) %>% 
    t() %>% 
    gridExtra::tableGrob()
  
  # Layout maps
  maps <- likelihood_plot + like_plot + grob1 + binned_map +
    plot_layout(ncol = 2, nrow = 2) +
    plot_annotation(
      title = title
    )
  
  # Save
  ggsave(plot_fp, maps, width=12, height=8, dpi=120)
  
}

# Create PNGs of likelihood maps ----
fps <- list.files(file.path(pred_dir, 'likelihood'), '*.tif', full.names=T)

# filter filepaths to species list
sp_fps <- sp_nospc_list %>% 
  map(~str_subset(.y, str_glue("{.x}.tif")), fps) %>% 
  flatten_chr()

# testing
lklhd_fp <- fps[[1]]
for (lklhd_fp in fps) {
  
  sp_nospc <- tools::file_path_sans_ext(basename(lklhd_fp))
  
  # Filepaths
  plot_fp <- file.path(rf_fig_dir, 'map_predictions', str_c(sp_nospc, '_likelihood.png'))
  dir.create(dirname(plot_fp), recursive = T, showWarnings = F)
  
  # Get observation points
  sp_dat <- pol_df2 %>% filter(str_detect(species, str_replace(sp_nospc, '_', ' ')))
  sp_df <- sp_dat$data[[1]]
  sp_ch <- sp_dat$convhull[[1]]
  
  # Load TIFF as stars
  pr_rf1_stars <- read_stars(lklhd_fp)

  # Plot
  likelihood_plot <- ggplot() +
    geom_stars(data=pr_rf1_stars) +
    geom_sf(data = mex, fill = "transparent", size = 0.2, color = "black") +
    colormap::scale_fill_colormap("Occupancy\nlikelihood", na.value = "transparent", 
                                  colormap = colormap::colormaps$viridis) +
    ggthemes::theme_hc() +
    theme(legend.position=c(.95, 1), legend.title.align=0, legend.justification = c(1,1)) +
    labs(x = NULL, y = NULL)
  
  like_plot <- likelihood_plot +
    geom_sf(data = sp_ch, fill = "transparent", size = 0.2, color = "white") +
    geom_sf(data = sp_df)
  
  # Save
  ggsave(plot_fp, like_plot, width=9, height=5.7, dpi=120)
  
}

# Create PNGs of presence maps ----
fps <- list.files(file.path(pred_dir, 'binned_spec_sens'), '*.tif', full.names=T)
for (fp in fps) {
  
  # Filepaths
  fn <- basename(fp)
  sp_nospc <- tools::file_path_sans_ext(fn)
  plot_fp <- file.path(rf_fig_dir, 'map_predictions', str_c(sp_nospc, '_bin_specsens.png'))
  dir.create(dirname(plot_fp), recursive = T, showWarnings = F)
  
  # Load TIFF as stars
  pa_rf1_stars <- read_stars(fp)
  
  # Plot
  binned_map <- ggplot() +
    geom_stars(data=pa_rf1_stars) +
    geom_sf(data = mex, fill = "transparent", size = 0.2, color = "black") +
    colormap::scale_fill_colormap("likely present", na.value = "transparent", 
                                  colormap = colormap::colormaps$viridis) +
    ggthemes::theme_hc() +
    theme(legend.position=c(.95, 1), legend.title.align=0, legend.justification = c(1,1))+
    labs(x = NULL, y = NULL) 
  
  # Save
  ggsave(plot_fp, binned_map, width=9, height=5.7, dpi=120)
  
}

# ~ COMBINE Sum likelihood maps ----
# list species (such as nocturnal butterflies)
if(pol_group == 'Mariposas'){
  sp_groups <- pol_df1 %>% st_drop_geometry %>% distinct(species, nocturna)
  list1 <- sp_groups %>%
    filter(nocturna=='nocturna') %>%
    transmute(species = str_replace_all(species, ' ', '_')) %>% 
    deframe
  
  l1_pattern <- str_c(list1, collapse='|')
  
  fps <- list.files(file.path(pred_dir, 'likelihood'), 'tif$', full.name=T) %>% 
    str_subset(str_glue('{l1_pattern}'))
}

# List TIFs and filter to species list
fps <- list.files(file.path(pred_dir, 'likelihood'), 'tif$', full.name=T)
sp_fps <- sp_nospc_list %>% 
  map(~str_subset(.y, str_glue("{.x}.tif")), fps) %>% 
  flatten_chr()

if(length(sp_fps) > 0) {
  
  rich_fn <- str_glue('richness_{name}_{rf_name}_{length(sp_fps)}species{apis_code}_lkhd')
  rich_plot_fp <- file.path(rf_fig_dir, str_glue('{rich_fn}.png'))
  rich_tif_fp <- file.path(pred_dir, str_glue('{rich_fn}.tif'))
  
  stack_sdms(sp_fps, rich_tif_fp, rich_plot_fp, mex)
  
}

# Sum presence/absence maps ----
# Use raster package
fps <- list.files(file.path(pred_dir, 'binned_spec_sens'), 'tif$', full.name=T)
sp_fps <- sp_nospc_list %>% 
  map(~str_subset(.y, str_glue("{.x}.tif")), fps) %>% 
  flatten_chr()

if(length(fps) > 0) {
  
  fn <- str_glue('richness_{name}_{rf_name}_{length(sp_fps)}species_binned')
  binned_rich_plot_fp <- file.path(rf_fig_dir, str_glue('{fn}.png'))
  binned_rich_tif_fp <- file.path(pred_dir, str_glue('{fn}.tif'))

  stack_sdms(sp_fps, binned_rich_tif_fp, binned_rich_plot_fp, mex)

}



# CULTIVO - sum to richness ----
# Likelihood

# List TIFs 
fps <- list.files(file.path(pred_dir, 'likelihood'), 'tif$', full.name=T)

# Filter TIFs to species list
sp_fps <- sp_nospc_list %>% 
  map(~str_subset(.y, str_glue("{.x}.tif")), fps) %>% 
  flatten_chr()

cult_dir <- 'data/data_out/pollinator_points/sdms_by_crop/Cucurbita_argyrosperma'
# cult_dir <- pred_dir

if(length(sp_fps) > 0) {
  
  rich_fp_prefix <- file.path(cult_dir, 
                              str_glue('Likhd_rich_{dfilt_code}_{length(sp_fps)}species'))
  rich_plot_fp <-str_glue('{rich_fp_prefix}.png')
  rich_tif_fp <- str_glue('{rich_fp_prefix}.tif')
  
  stack_sdms(sp_fps, rich_tif_fp, rich_plot_fp, mex)
  
}

# Presence/Absence
# List TIFs and filter to species list
fps <- list.files(file.path(pred_dir, 'binned_spec_sens'), 'tif$', full.name=T)
sp_fps <- sp_nospc_list %>% 
  map(~str_subset(.y, str_glue("{.x}.tif")), fps) %>% 
  flatten_chr()

# Run for crop
if(length(fps) > 0) {
  
  rich_fp_prefix <- file.path(cult_dir, 
                              str_glue('PA_rich_{dfilt_code}_{length(sp_fps)}species'))
  rich_tif_fp <- str_glue('{rich_fp_prefix}.tif')
  rich_plot_fp <-str_glue('{rich_fp_prefix}.png')
  
  stack_sdms(sp_fps, rich_tif_fp, rich_plot_fp, mex)
  
}



# COMPARE model versions ----
taxa <- 'Apis_mellifera'
fps <- list.files('data/data_out/sdm', str_glue('{taxa}\\.csv'),
           full.names = TRUE, recursive = TRUE) %>% 
  str_subset('model_evals')
fps

evals <- fps %>% map_dfr(read_csv)
# evals$fp <- fps %>% str_extract("(?<=\\/)rf[^\\/]*")
evals$fp <- fps %>% str_extract("(?<=\\/)rf.*(?=\\/)")

taxa <- 'Bombus_ephippiatus'
fps <- list.files('data/data_out/sdm', str_glue('{taxa}\\.csv'),
                  full.names = TRUE, recursive = TRUE) %>% 
  str_subset('model_evals')
evals <- fps %>% map_dfr(read_csv)
# evals$fp <- fps %>% str_extract("(?<=\\/)rf[^\\/]*")
evals$fp <- fps %>% str_extract("(?<=\\/)rf.*(?=\\/)")

taxa <- 'Partamona_bilineata'
fps <- list.files('data/data_out/sdm', str_glue('{taxa}\\.csv'),
                  full.names = TRUE, recursive = TRUE) %>% 
  str_subset('model_evals')
evals <- fps %>% map_dfr(read_csv)
# evals$fp <- fps %>% str_extract("(?<=\\/)rf[^\\/]*")
evals$fp <- fps %>% str_extract("(?<=\\/)rf.*(?=\\/)")
