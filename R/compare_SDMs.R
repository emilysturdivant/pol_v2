# Input: points downloaded from GBIF in folder input_data/GBIF/family_order_query
# Output: files in folder for current iteration of RF: data_out/sdm/rfX_params/query_term

# Load libraries ----
source('R/00_initialize.R')

# Load filtered points ----
# filt_pts_rds <- file.path(pts_dir, str_c('points_nested_species_filt.rds'))
# pol_df2 <- readRDS(filt_pts_rds) %>% 
#   filter(genus != "", species != "") %>% 
#   filter(!is.na(genus), !is.na(species))
# 
# # Look at model statistics together
# fps <- list.files(file.path(pred_dir, 'model_evals'), '*.csv$', full.names = T)
# eval_tbl_fp <- file.path(pred_dir, str_c('model_evals_', length(fps), 'species.csv'))
# if (file.exists(eval_tbl_fp)) {
#   eval_tbl <- read_csv(eval_tbl_fp)
#   pol_df3 <- pol_df2 %>% left_join(eval_tbl, by = 'species')
# }
# 
# # Add filenames for likelihood tifs 
# fps <- list.files(file.path(pred_dir, 'likelihood'), '*.tif$', full.names = T) %>% 
#   as_tibble_col('lklhd_tif') %>% 
#   mutate(species = basename(tools::file_path_sans_ext(lklhd_tif)) %>% 
#            str_replace('_', ' '))
# pol_df3 <- pol_df3 %>% left_join(fps, by = 'species')
# pol_df3 %>% saveRDS(file.path(pred_dir, 'pol_df3.rds'))

# Load
pred_dir1 <- get_pred_dir(contin_vars_only = TRUE)
rf1 <- readRDS(file.path(pred_dir1, 'pol_df3.rds')) %>% 
  mutate(success = if_else(!is.na(lklhd_tif), 'success', 'fail'))

pred_dir3 <- get_pred_dir(contin_vars_only = FALSE, rf_vers = 3)
rf3 <- readRDS(file.path(pred_dir3, 'pol_df3.rds')) %>% 
  mutate(success = if_else(!is.na(lklhd_tif), 'success', 'fail'))

# Look at overview of success by group
rf1_success_counts <- rf1 %>% 
  count(common_group, subgroup_a, success) %>% 
  pivot_wider(names_from = success, values_from = n)

rf3_success_counts <- rf3 %>% 
  count(common_group, subgroup_a, success) %>% 
  pivot_wider(names_from = success, values_from = n)

full_join(rf1_success_counts, rf3_success_counts, 
          by = c('common_group', 'subgroup_a'),
          suffix = c('_rf1', '_rf3'))

# List species that failed in rf3
rf3 %>% tbl_vars
sp_list <- rf3 %>% 
  filter(success != 'success',
         common_group == 'Abejas', 
         subgroup_a == 'Eusocial') %>% 
  select(species) %>% 
  deframe()

rf1 %>% 
  filter(species %in% sp_list) %>% 
  select(species, nobs, N_unq_pts, N_unq_cells)

sp_dlist <- rf1 %>% 
  filter(species %in% sp_list) %>% 
  select(species, data) 

sp_dlist %>% purrr::walk(~ print(.x$species))
spt <- sp_dlist[[1]]
  
function(sp_df) {
  sp_name <- sp_dat$species
  sp_nospc <- sp_name %>% str_replace(' ', '_')
  
  pts_fp <- file.path(pred_dir, 'species_points', str_c(sp_nospc, '.gpkg'))
  dir.create(dirname(pts_fp), recursive = TRUE, showWarnings = FALSE)
  
  # Get observation points
  sp_df <- sp_dat$data[[1]]
  
  if(!file.exists(pts_fp)) sp_df %>% st_write(pts_fp)
}
