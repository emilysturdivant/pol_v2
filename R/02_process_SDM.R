# print(getwd())

# Set working directory if on LANASE cluster
sysinfo <- Sys.info()
if(sysinfo[['nodename']] == 'lanase.cluster') {
  setwd('/share/Part1/esturdivant/pol_v2')
}

# Initialize
source('R/00_initialize.R')
save_pts <- FALSE

# Load filtered points ----
pol_df2 <- readRDS(filt_pts_rds) 

# Load environment variables ----
pred_grd <- file.path('data/tidy/environment_variables', str_c('pred_', var_code, '.grd'))
mex <- st_read(mex_fp)
mex0 <- st_union(mex)
pred <- prep_predictor_stack(pred_grd, crop_dir, vars, mex0, overwrite = FALSE)

# RF model for each species ----
# Add filenames for likelihood tifs 
fps <- list.files(file.path(pred_dir, 'models'), '*.rds$', full.names = T) %>% 
  as_tibble_col('mod_fp') %>% 
  mutate(species = basename(tools::file_path_sans_ext(mod_fp)) %>% 
           str_replace('_', ' '))
pol_df3 <- pol_df2 %>% 
  left_join(fps, by = 'species') %>% 
  filter(is.na(mod_fp)) %>% 
  slice(1:3) %>%
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

# Look at model statistics together
fps <- list.files(file.path(pred_dir, 'model_evals'), '*.csv$', full.names = T)
eval_tbl_fp <- file.path(pred_dir, str_c('model_evals_', length(fps), 'species.csv'))
if(!file.exists(eval_tbl_fp)){
  if(length(fps) > 0) {
    eval_tbl <- fps %>% purrr::map_dfr(read.csv)
    eval_tbl %>% write_csv(eval_tbl_fp)
  }
}
eval_tbl <- read_csv(eval_tbl_fp)
pol_df4b <- pol_df3 %>% left_join(eval_tbl, by = 'species')

pol_df4b %>% 
  select(species, nobs, mod_fp, lklhd_tif, N_unq_pts, N_unq_cells, np, na, auc, cor, pcor, thresh_spec_sens) %>% 
  filter()

# Optionally save points ----
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
