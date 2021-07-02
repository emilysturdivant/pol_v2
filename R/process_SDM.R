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
library('dismo')
library('raster') # use raster for compatibility with dismo
library('stars')
library('sf')
library('tidyverse')

# Initialize -----
unq_cells = TRUE
mutually_exclusive_pa = TRUE
filt_dates = FALSE
nspecies <- 15
rf_vers <- 2
name <- order <- 'Hymenoptera'
contin_vars_only <- TRUE
# contin_vars_only <- FALSE
noApis = TRUE

# Date range
date_range <- c(2000, 2020)
# pol_groups <- c('Abejas', 'Avispas', 'Colibries', 'Mariposas', 'Moscas', 'Murcielagos')

# Functions ----
#' @export
model_species_rf <- function(sp_df,
                             pred, 
                             model_fp, 
                             sp_name,
                             eval_fp,
                             erf_fp,
                             mutually_exclusive_pa=TRUE,
                             unq_cells=TRUE,
                             rf_fig_dir=NULL) {
  
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
  rf1 <- randomForest::randomForest(
    pa ~ . -pa,
    data = envtrain,
    na.action = na.exclude,
    importance=T, 
    ntree=1000)
  
  # Save
  dir.create(dirname(model_fp), recursive = T, showWarnings = F)
  saveRDS(rf1, model_fp)
  # rf1 <- readRDS(model_fp)
  
  # Filenames
  dir.create(dirname(eval_fp), recursive = T, showWarnings = F)
  
  # Evaluate model with test data
  erf <- dismo::evaluate(testpres, testbackg, rf1)
  
  # Save model evaluation
  saveRDS(erf, erf_fp)
  
  # Save simple statistics in CSV
  spc_eval <- tibble(
    species=sp_name,
    N_unq_pts=nrow(sp_df),
    N_unq_cells=nrow(filter(envtrain, pa == 1)),
    np=erf@np, na=erf@na, auc=erf@auc,
    cor=erf@cor, pcor=erf@pcor, 
    spec_sens=dismo::threshold(erf, "spec_sens"))
  
  spc_eval %>% write_csv(eval_fp)
  
  if(!is.null(rf_fig_dir)){
    # Get plot directory
    plot_fp <- file.path(rf_fig_dir, 'var_importance', str_c(sp_nospc, '.png'))
    dir.create(dirname(plot_fp), recursive = T, showWarnings = F)
    
    # Save plot
    png(plot_fp)
    randomForest::varImpPlot(rf1, type=1, sort=F, 
                             main=sp_name,
                             pt.cex=1,
                             bg='black')
    dev.off()
  }
  
  # Return
  return(erf)
}

#' @export
predict_distribution_rf <- function(rf_fp, erf_fp, likelihood_fp, binned_fp, ext){
  
  # Load model
  if(is.character(rf_fp)){
    rf1 <- readRDS(rf_fp)
  } else {
    rf1 <- rf_fp
  }
  
  # Create map and interpolate to fill holes
  pr_rf1 <- dismo::predict(predictors, rf1, ext=ext)
  pr_rf1 <- raster::focal(pr_rf1, 
                          w=matrix(1,nrow=3, ncol=3), 
                          fun=mean, 
                          NAonly=TRUE, 
                          na.rm=TRUE) 
  
  # Save likelihood raster
  dir.create(dirname(likelihood_fp), recursive = T, showWarnings = F)
  writeRaster(pr_rf1, likelihood_fp, overwrite=T, 
              options=c("dstnodata=-99999"), wopt=list(gdal='COMPRESS=LZW'))
  
  # Apply threshold from max TPR+TNR and save
  erf <- readRDS(erf_fp)
  tr <- threshold(erf, 'spec_sens')
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

# directory paths ----
unq_code <- ifelse(unq_cells, 'unq_cells', 'unq_pts')
unq_code <- ifelse(mutually_exclusive_pa, 'exclusive', unq_code)
dfilt_code <- ifelse(filt_dates, '2000to2020', 'alldates')
cont_var_code <- ifelse(contin_vars_only, '_continuouspredictors', '')
apis_code <- ifelse(order == 'Hymenoptera', ifelse(noApis, '_noApis', ''), '')

rf_name <- str_c('rf', str_c(rf_vers, unq_code, dfilt_code, sep='_'), cont_var_code)
fp_tail <- file.path('sdm', rf_name, name)
pred_dir <- file.path('data', 'data_out', fp_tail)
rf_fig_dir <- file.path('figures', fp_tail)
rf_fig_dir <- pred_dir
dir.create(pred_dir, recursive=T, showWarnings = F)
dir.create(rf_fig_dir, recursive=T, showWarnings = F)

# Input GBIF points
pts_fp <- str_glue("data/input_data/GBIF/family_order_query/gbif_{order}.rds")

# Output filtered points
filt_pts_rds <- file.path(pred_dir, str_c('filtered_pts_', name, '.rds'))

# Load and filter GBIF points for given name ----
if( !file.exists(filt_pts_rds) ) {
  
  # Load raw GBIF points
  dat <- readRDS(pts_fp)
  
  # Data tidying steps: drop imprecise coordinates, drop duplicates
  vars <- c('species', 'genus', 'tribe', 'family', 'superfamily', 'order', 'class', 'nocturna', 
            'decimalLongitude', 'decimalLatitude', 
            'eventDate', 'coordinateUncertaintyInMeters', 'habitat', 
            'basisOfRecord', 'country', 'stateProvince', 'institutionCode')
  
  df <- dat %>% 
    filter(coordinateUncertaintyInMeters < 1000 | is.na(coordinateUncertaintyInMeters)) %>% 
    filter(genus != "") %>% 
    dplyr::select(matches(vars))
  
  # drop duplicates
  pol_df1 <- df %>% 
    distinct %>% 
    st_as_sf(x = .,                         
             coords = c("decimalLongitude", "decimalLatitude"),
             crs = 4326)
  
  # Optionally filter to date range
  if(filt_dates) {
    # Dates
    date_min <- as.POSIXct(str_c(date_range[[1]], "-01-01"))
    date_max <- as.POSIXct(str_c(date_range[[2]], "-12-31"))
    
    pol_df1 <- pol_df1 %>% 
      filter(eventDate >= date_min & eventDate <= date_max)
  }
  
  # Remove species with less than 25 distinct observations (based on Koch et al. 2017)
  pol_df2 <- pol_df1 %>% 
    st_transform(st_crs(predictors)) %>% 
    mutate(species = ifelse(species == "", genus, species)) %>% 
    group_by(species, genus) %>% 
    filter(n() > 24) %>% 
    ungroup()
  
  # Save
  pol_df2 %>% saveRDS(filt_pts_rds)
  
}

# # Merge filtered points into one file
# pol_dir_out <- file.path(pol_dir, str_c(dfilt_code, '_gt24perSpecies'))
# (fps <- list.files(pol_dir_out, 'gpkg$', full.names=T))
# combo <- fps %>% map_dfr(st_read)
# combo %>% st_write(file.path(pol_dir, 'combined', str_c(dfilt_code, '_gt24perSpecies.gpkg')))

# Load filtered points ----
pol_df2 <- readRDS(filt_pts_rds)

# Get species list 
sp_counts <- pol_df2 %>% 
  st_drop_geometry %>% 
  count(species, genus) %>% 
  arrange(desc(n))

if(noApis) {
  sp_counts <- sp_counts %>% 
    filter(!str_detect(species, 'Apis mellifera'))
}

sp_list <- sp_counts %>% 
  slice(1:nspecies) %>% 
  dplyr::select(species) %>% 
  deframe
sp_nospc_list <- str_replace(sp_list, ' ', '_')

# ~ Make RF model for each species ----
# sp_name <- sp_list[[2]]
for(sp_name in sp_list) {
  
  print(sp_name)
  
  # Filepath
  sp_nospc <- str_replace(sp_name, ' ', '_')
  model_fp <- file.path(pred_dir, 'models', str_c(sp_nospc, '.rds'))
  if (file.exists(model_fp)) {
    print(str_c(sp_name, ': Model already created.'))
    next
  }
  
  # Filter to species
  sp_df <- pol_df2 %>% filter(species == sp_name)
  
  if(nrow(sp_df) < 1) {
    print(str_c(sp_name, ':No rows for the given species.'))
    next
  }
  
  # File paths
  eval_fp <- file.path(pred_dir, 'model_evals', str_c(sp_nospc, '.csv'))
  erf_fp <- file.path(pred_dir, 'model_evals', str_c(sp_nospc, '.rds'))
  
  erf <- model_species_rf(sp_df,
                          pred, 
                          model_fp, 
                          sp_name,
                          eval_fp,
                          erf_fp,
                          mutually_exclusive_pa,
                          unq_cells,
                          rf_fig_dir)
}

# Look at model statistics together
fps <- list.files(file.path(pred_dir, 'model_evals'), '*.csv$', full.names = T)
if(length(fps) > 0) {
  eval_tbl <- fps %>% purrr::map_dfr(read.csv)
  eval_tbl_fp <- file.path(pred_dir, str_c('model_evals_', length(fps), 'species.csv'))
  eval_tbl %>% write_csv(eval_tbl_fp)
}

# Save TIFs of likelihood and presence/absence ----
fps <- list.files(file.path(pred_dir, 'models'), '*.rds', full.names = T)

# filter filepaths to species list
sp_fps <- sp_nospc_list %>% 
  map(~str_subset(.y, str_glue("{.x}.rds")), fps) %>% 
  flatten_chr()

for(rf_fp in sp_fps){
  
  # Filepaths  
  sp_nospc <- tools::file_path_sans_ext(basename(rf_fp))
  print(sp_nospc)
  
  likelihood_fp <- file.path(pred_dir, 'likelihood', str_glue(sp_nospc, '.tif'))
  binned_fp <- file.path(pred_dir, 'binned_spec_sens', str_glue(sp_nospc, '.tif'))
  
  if(file.exists(likelihood_fp) & file.exists(binned_fp)){
    next
  } 
  
  erf_fp <- file.path(pred_dir, 'model_evals', str_c(sp_nospc, '.rds'))
  
  # Make TIFs of likelihood and presence/absence
  predict_distribution_rf(rf_fp, erf_fp, likelihood_fp, binned_fp, ext)
  
}

# Create PNGs of likelihood maps ----
fps <- list.files(file.path(pred_dir, 'likelihood'), '*.tif', full.names=T)

# filter filepaths to species list
sp_fps <- sp_nospc_list %>% 
  map(~str_subset(.y, str_glue("{.x}.tif")), fps) %>% 
  flatten_chr()

for (fp in fps) {
  
  # Filepaths
  sp_nospc <- tools::file_path_sans_ext(basename(fp))
  plot_fp <- file.path(rf_fig_dir, 'map_predictions', str_c(sp_nospc, '_likelihood.png'))
  dir.create(dirname(plot_fp), recursive = T, showWarnings = F)
  
  # Load TIFF as stars
  pr_rf1_stars <- read_stars(fp)

  # Plot
  likelihood_plot <- ggplot() +
    geom_stars(data=pr_rf1_stars) +
    geom_sf(data = mex, fill = "transparent", size = 0.2, color = "black") +
    colormap::scale_fill_colormap("Occupancy\nlikelihood", na.value = "transparent", 
                                  colormap = colormap::colormaps$viridis) +
    ggthemes::theme_hc() +
    theme(legend.position=c(.95, 1), legend.title.align=0, legend.justification = c(1,1)) +
    labs(x = NULL, y = NULL)
  
  # Save
  ggsave(plot_fp, likelihood_plot, width=9, height=5.7, dpi=120)
  
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
