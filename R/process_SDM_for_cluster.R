# print(getwd())
setwd('/share/Part1/esturdivant/pol_v2')

# Initialize
source('R/initialize.R')

# Load filtered points ----
filt_pts_rds <- file.path(pts_dir, str_c('points_nested_species_filt.rds'))
pol_df2 <- readRDS(filt_pts_rds) %>% 
  filter(genus != "", species != "") %>% 
  filter(!is.na(genus), !is.na(species))

# Load environment variables ----
pred_grd <- file.path('data/tidy/environment_variables', str_c('pred_', var_code, '.grd'))
pred <- prep_predictor_stack(pred_grd, crop_dir, vars, mex0, overwrite = TRUE)

# RF model for each species ----
pol_df3 <- pol_df2 %>% 
  slice(120:122) %>% 
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
