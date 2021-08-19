# v2: using GBIF points that I downloaded
# Input: points downloaded from GBIF in folder input_data/GBIF/family_order_query
# Output: files in folder for current iteration of RF: data_out/sdm/rfX_params/query_term

# Load libraries ----
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

# ~ RF model for each species ----
# Model each species
save_pts <- FALSE
df <- pol_df2 #%>% filter(common_group == 'Moscas')
# stop <- nrow(df)
# stop <- 2
# for (i in seq(1, stop)) {
#   
#   # Get data
#   sp_dat <- df %>% slice(i)
#   sp_name <- sp_dat$species
#   sp_nospc <- sp_name %>% str_replace(' ', '_')
#   print(str_glue("\n\n{i} of {stop}: {sp_name}"))
#     
#   # Get observation points
#   sp_df <- sp_dat$data[[1]]
#   
#   # Save points
#   if(save_pts) {
#     pts_fp <- file.path(pred_dir, 'species_points', str_c(sp_nospc, '.gpkg'))
#     dir.create(dirname(pts_fp), recursive = TRUE, showWarnings = FALSE)
#     if(!file.exists(pts_fp)) sp_df %>% st_write(pts_fp)
#   }
# 
#   # Run random forest and model evaluation
#   mod_rf <- model_species_rf(sp_df,
#                           pred, 
#                           pred_dir,
#                           sp_name,
#                           excludep,
#                           unq_cells)
#   
#   if(is_empty(mod_rf)) next
#   
#   # Use model to create prediction rasters
#   predict_distribution_rf(mod_rf$rf, mod_rf$erf, sp_name, pred_dir, ext, pred, write_binned = TRUE)
#   
# }

df1 <- df %>% 
  slice(7:20) %>% 
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
pol_df3 <- pol_df2 %>% left_join(eval_tbl, by = 'species')

# pngs for QC ----
(sp_row <- pol_df2 %>% sample_n(1))
(sp_row <- pol_df2 %>% filter(species == 'Alypiodes bimaculata'))
(sp_row <- pol_df3 %>% arrange(auc) %>% slice(1))
(sp_row <- pol_df3 %>% arrange(N_unq_cells) %>% slice(6))
(sp_row <- pol_df3 %>% filter(N_unq_cells < 20) %>% sample_n(1))

dir.create(file.path(rf_fig_dir, 'map_predictions'), recursive = T, showWarnings = F)
plot_qc_maps(sp_row, file.path(rf_fig_dir, 'map_predictions'))

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
#   plot_fp <- file.path(rf_fig_dir, 'map_predictions', str_c(sp_nospc, '_likelihood.png'))
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
#   plot_fp <- file.path(rf_fig_dir, 'map_predictions', str_c(sp_nospc, '_bin_specsens.png'))
#   dir.create(dirname(plot_fp), recursive = T, showWarnings = F)
#   
#   binned_map <- plot_binned_map(fp, binned_map)
#   
# }

# ~ COMBINE Sum likelihood maps ----
# Add filenames for likelihood tifs 
fps <- list.files(file.path(pred_dir, 'likelihood'), '*.tif$', full.names = T) %>% 
  as_tibble_col('lklhd_tif') %>% 
  mutate(species = basename(tools::file_path_sans_ext(lklhd_tif)) %>% 
           str_replace('_', ' '))
pol_df3 <- pol_df3 %>% left_join(fps, by = 'species')

# pol_df3 %>% distinct(common_group, subgroup_a, subgroup_b)
pol_df3 %>% distinct(common_group, subgroup_a)

# Get file list for given group
pol_group <- 'Abejas'
subgrp_a <- 'solitaria'
apis_code <- ''

name <- str_c(str_to_title(str_sub(c(pol_group, subgrp_a), end=4)), collapse = '')

sub_df <- pol_df3 %>% 
  filter(common_group == pol_group, !is.na(lklhd_tif))

if(subgrp_a != '') {
  sub_df <- sub_df %>% 
    filter(str_detect(subgroup_a, regex(subgrp_a, ignore_case = TRUE))) 
}

sp_fps <- sub_df %>% 
  dplyr::select(lklhd_tif) %>% 
  deframe()

# Sum likelihood rasters
if(length(sp_fps) > 0) {
  
  rich_fn <- str_glue('richness_{name}_rf{rf_vers}_{length(sp_fps)}species{apis_code}_lkhd')
  rich_plot_fp <- file.path(rf_fig_dir, 'richness', str_glue('{rich_fn}.png'))
  rich_tif_fp <- file.path(pred_dir, 'richness',  str_glue('{rich_fn}.tif'))
  
  stack_sdms_terra(sp_fps, rich_tif_fp, rich_plot_fp, mex)
  
}



# IN PROGRESS - Perform for all (group_by) ----
pol_df4 <- pol_df3 %>%
  mutate(solitaria = str_detect(subgroup_a, regex('solitaria', ignore_case = TRUE)),
         social = str_detect(subgroup_a, regex('social', ignore_case = TRUE)),
         eusocial = str_detect(subgroup_a, regex('eusocial', ignore_case = TRUE)),
         diurna =  str_detect(subgroup_a, regex('diurna', ignore_case = TRUE)),
         nocturna = str_detect(subgroup_a, regex('nocturna', ignore_case = TRUE)))

pol_df4 %>% 
  group_by(common_group, solitaria) %>% 
  nest() %>% 
  ungroup()

row <- pol_df4 %>% slice(1)

pol_df4 %>% #tbl_vars()
  filter(common_group == pol_group, !is.na(lklhd_tif)) %>%
  pivot_longer(cols = contains(subgrp_a)) %>%
  filter(value) %>%
  dplyr::select(lklhd_tif)

