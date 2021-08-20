
# setwd('/share/Part1/esturdivant/pol_v2')

# Initialize
source('R/initialize.R')

# Load filtered points ----
filt_pts_rds <- file.path(pts_dir, str_c('points_nested_species_filt.rds'))
pol_df2 <- readRDS(filt_pts_rds) %>% 
  filter(genus != "", species != "") %>% 
  filter(!is.na(genus), !is.na(species))

# Load environment variables ----
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

# pred_tif <- file.path('data/tidy/environment_variables', str_c('pred_', var_code, '.tif'))
# pred %>% writeRaster(pred_tif)
# pred2 <- raster::stack(pred_tif)
# 
# pred_ncdf <- file.path('data/tidy/environment_variables', str_c('pred_', var_code, '.nc'))
# pred %>% writeRaster(pred_ncdf, overwrite = TRUE, zname = 'Name')
# pred2 <- raster::stack(pred_ncdf)
# 
# pred_ncdf <- file.path('data/tidy/environment_variables', str_c('pred2_', var_code, '.nc'))
# predb <- raster::brick(pred)
# names(predb) <- names(pred)
# predb %>% writeRaster(pred_ncdf, varname = names(predb))
# pred2 <- raster::stack(pred_ncdf)

# RF model for each species ----
pol_df3 <- pol_df2 %>% 
  slice(110:120) %>% 
  mutate(
    mod_fp = purrr::map2(
      data, species, 
      ~ model_species_rf(.x,
                         pred,
                         pred_dir,
                         .y,
                         excludep,
                         unq_cells, 
                         return_obj = 'paths')),
    lklhd_tif = purrr::map2_chr(
      mod_fp, species, 
      ~ predict_distribution_rf(.x, 
                                .y, 
                                pred_dir, 
                                ext,
                                pred, 
                                write_binned = FALSE))
  )

# Optionally save points ----
save_pts <- FALSE
if(save_points) {
  
  # Save species points in separate gpk and add filepath to df
  save_pts <- function(.x, .y) {
    fp <- file.path(pred_dir, 'species_points', 
                    str_c(str_replace(.y, ' ', '_'), '.gpkg'))
    if(!file.exists(fp)) st_write(.x, fp)
    return(fp)
  }
  
  # Run
  dir.create(file.path(pred_dir, 'species_points'), 
             recursive = TRUE, showWarnings = FALSE)
  df %>% mutate(pts_fp = purrr::map2_chr(data, species, save_pts))
  
}
