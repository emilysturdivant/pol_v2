# v2: using GBIF points that I downloaded
# Input: points downloaded from GBIF in folder input_data/GBIF/family_order_query
# Output: files in folder for current iteration of RF: data_out/sdm/rfX_params/query_term

# Load libraries ----
source('R/initialize.R')

# Load filtered points ----
pol_df2 <- readRDS(filt_pts_rds)

dat <- readRDS(pts_nested_rds)

# Look at model statistics together ----
eval_tbl_fp <- file.path(pred_dir, str_c('model_evals_', length(fps), 'species.csv'))
if(!file.exists(eval_tbl_fp)){
  if(length(fps) > 0) {
    eval_tbl <- fps %>% purrr::map_dfr(read.csv)
    eval_tbl %>% write_csv(eval_tbl_fp)
  }
}

# Get filepath for model_evals table with most species
eval_tbl_fp <- list.files(pred_dir, 'model_evals.*\\.csv$', full.names = T) %>% 
    as_tibble_col() %>% 
    mutate(ct = as.integer(str_extract_all(value, '\\d+(?=species)'))) %>% 
    slice_max(ct) %>% 
    nth(1)
eval_tbl <- read_csv(eval_tbl_fp)
pol_df3 <- pol_df2 %>% left_join(eval_tbl, by = 'species')

# ~ Sum likelihood maps ----
# Add filenames for likelihood tifs 
fps <- list.files(file.path(pred_dir, 'likelihood'), '*.tif$', full.names = T) %>% 
  as_tibble_col('lklhd_tif') %>% 
  mutate(species = basename(tools::file_path_sans_ext(lklhd_tif)) %>% 
           str_replace('_', ' '))
pol_df3 <- pol_df3 %>% left_join(fps, by = 'species')

# Sum richness for each pollinator group ----
iterate_stack_sdms <- function(fps, name, mex) {
  
  # Load reference polygons
  if(is.character(mex)) mex <- st_read(mex)
  # mex <- st_read(mex_fp)
  
  fps <- fps[[1]]
  
  # Output filenames
  rich_fn <- str_glue('richness_{name}_rf{rf_vers}_{length(fps)}species_lkhd')
  rich_plot_fp <- file.path(pred_dir, 'richness', str_glue('{rich_fn}.png'))
  rich_tif_fp <- file.path(pred_dir, 'richness',  str_glue('{rich_fn}.tif'))
  
  # Sum likelihood rasters
  stack_sdms(fps, rich_tif_fp, rich_plot_fp, mex)
}

# Create nested DF
pol_nested <- pol_df3 %>% 
  select(group, lklhd_tif) %>% 
  group_by(group) %>% 
  nest(fps = c(lklhd_tif)) %>% 
  ungroup()

purrr::walk2(pol_nested$data, pol_nested$group, iterate_stack_sdms, mex)


# pol_df3 %>% distinct(group, bee_sociality, bee_nesting)
pol_df3 %>% distinct(group)

sp_fps <- pol_df3 %>% 
  filter(group == 'bee', !is.na(lklhd_tif)) %>% 
  dplyr::select(lklhd_tif) %>% 
  deframe()

# Create nested DF
pol_nested <- pol_df3 %>% 
  filter(!is.na(bee_sociality)) %>% 
  select(bee_sociality, lklhd_tif) %>% 
  group_by(bee_sociality) %>% 
  nest(fps = c(lklhd_tif)) %>% 
  ungroup()

purrr::walk2(pol_nested$fps, pol_nested$bee_sociality, iterate_stack_sdms, mex)

# Create nested DF
pol_nested <- pol_df3 %>% 
  filter(!is.na(bee_nesting)) %>% 
  select(bee_nesting, lklhd_tif) %>% 
  group_by(bee_nesting) %>% 
  nest(fps = c(lklhd_tif)) %>% 
  ungroup()

purrr::walk2(pol_nested$fps, pol_nested$bee_nesting, iterate_stack_sdms, mex)




# pngs for QC ----
(sp_row <- pol_df3 %>% sample_n(1))
(sp_row <- pol_df3 %>% filter(species == 'Alypiodes bimaculata'))
(sp_row <- pol_df3 %>% arrange(auc) %>% slice(1))
(sp_row <- pol_df3 %>% slice_max(N_unq_cells, n = 6) %>% sample_n(1))
(sp_row <- pol_df3 %>% filter(N_unq_cells < 25) %>% sample_n(1))

dir.create(file.path(pred_dir, 'map_predictions'), recursive = T, showWarnings = F)
plot_qc_maps(sp_row, file.path(pred_dir, 'map_predictions'))

# # Compare accuracy metrics by threshold ----
# fps <- list.files(file.path(pred_dir, 'models'), '*.rds', full.names = T)
# rf_fp <- fps[[10]]
# sp_nospc <- tools::file_path_sans_ext(basename(rf_fp))
# erf_fp <- file.path(pred_dir, 'model_evals', str_c(sp_nospc, '.rds'))
# erf <- readRDS(erf_fp)
# 
# c('spec_sens', 'kappa', 'no_omission', 'sensitivity') %>% 
#   map_dfr(get_accuracy_metrics, erf)

# 
# # pngs of likelihood maps ----
# fps <- list.files(file.path(pred_dir, 'likelihood'), '*.tif', full.names=T)
# 
# # filter filepaths to species list
# sp_fps <- sp_nospc_list %>% 
#   map(~str_subset(.y, str_glue("{.x}.tif")), fps) %>% 
#   flatten_chr()
# 
# # testing
# lklhd_fp <- fps[[1]]
# for (lklhd_fp in fps) {
#   
#   sp_nospc <- tools::file_path_sans_ext(basename(lklhd_fp))
#   
#   # Filepaths
#   plot_fp <- file.path(pred_dir, 'map_predictions', str_c(sp_nospc, '_likelihood.png'))
#   dir.create(dirname(plot_fp), recursive = T, showWarnings = F)
#   
#   # Get observation points
#   sp_dat <- pol_df2 %>% filter(str_detect(species, str_replace(sp_nospc, '_', ' ')))
#   sp_df <- sp_dat$data[[1]]
#   sp_ch <- sp_dat$convhull[[1]]
#   
#   # Plot
#   likelihood_plot <- plot_lklhd_map(lklhd_fp)
#   
#   # Add convex hull and presence points
#   like_plot <- likelihood_plot +
#     geom_sf(data = sp_ch, fill = "transparent", size = 0.2, color = "white") +
#     geom_sf(data = sp_df)
#   
#   # Save
#   ggsave(plot_fp, like_plot, width=9, height=5.7, dpi=120)
#   
# }
# 
# # pngs of presence maps ----
# fps <- list.files(file.path(pred_dir, 'binned_spec_sens'), '*.tif', full.names=T)
# for (fp in fps) {
#   
#   # Filepaths
#   fn <- basename(fp)
#   sp_nospc <- tools::file_path_sans_ext(fn)
#   plot_fp <- file.path(pred_dir, 'map_predictions', str_c(sp_nospc, '_bin_specsens.png'))
#   dir.create(dirname(plot_fp), recursive = T, showWarnings = F)
#   
#   binned_map <- plot_binned_map(fp, binned_map)
#   
# }

