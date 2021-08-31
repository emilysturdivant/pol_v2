
# Load libraries ----
source('R/00_initialize.R')

# Load filtered points ----
pol_df2 <- readRDS(filt_pts_rds)
pol_df2 %>% distinct(bee_sociality)
pol_df2 <- pol_df2 %>%
  mutate(bee_sociality = str_replace_all(bee_sociality, 'solitary/social', 'solitary-social'))

dat <- readRDS(pts_nested_rds)

# Look at model statistics together ----
fps <- list.files(file.path(pred_dir, 'model_evals'), '*.csv$', full.names = T)
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
stack_sdms_by_group <- function(df, out_dir, grp_var, mex, rf_vers) {

  rich_df <- df %>% 
    select(.data[[grp_var]], lklhd_tif) %>% 
    group_by(.data[[grp_var]]) %>% 
    nest(fps = c(lklhd_tif)) %>% 
    ungroup()
  
  rich_df <- rich_df %>% 
    mutate(rich_tif = purrr::map2_chr(fps, .data[[grp_var]], 
                                      stack_sdms,
                                      out_dir, 
                                      mex, 
                                      rf_vers),
           species_ct = purrr::map_dbl(fps, nrow))
  
  acc_df <- df %>% 
    group_by(.data[[grp_var]]) %>% 
    summarize(auc_mean = mean(auc, na.rm = TRUE))
  
  full_join(rich_df, acc_df)
}

# Richness for all pollinator groups ----
out_dir <- file.path(pred_dir, 'richness')
pol_richness <- stack_sdms_by_group(pol_df3, out_dir, 'group', mex, rf_vers)

purrr::walk2(pol_richness$rich_tif, 
             pol_richness$species_ct,
             plot_richness,
             mex)

# Richness by bee sociality ----
out_dir <- file.path(pred_dir, 'richness', 'bee_subgroups')
bee_richness_soc <- pol_df3 %>% 
  filter(!is.na(bee_sociality)) %>% 
  stack_sdms_by_group(out_dir, 'bee_sociality', mex, rf_vers)

purrr::walk2(bee_richness_soc$rich_tif, 
            bee_richness_soc$species_ct,
              plot_richness,
              mex)

# Richness by bee nesting ----
out_dir <- file.path(pred_dir, 'richness', 'bee_subgroups')
bee_richness_nest <- pol_df3 %>% 
  filter(!is.na(bee_nesting)) %>% 
  stack_sdms_by_group(out_dir, 'bee_nesting', mex, rf_vers)

purrr::walk2(bee_richness_nest$rich_tif, 
             bee_richness_nest$species_ct,
             plot_richness,
             mex)

# Specific richness map ----
# Plot richness
anps <- st_read(anp_terr_fp)
species_ct <- 243
rich_tif <- list.files(file.path(pred_dir, 'richness', 'bee_subgroups'),
                       'stems.*\\.tif', full.names = TRUE)
rich_plot_fp <- file.path(pred_dir, 'richness', 'requested_maps', 'rf4_stemnesting_bees24lkhd_anps.png')
plot_richness(rich_tif, species_ct, anps, rich_plot_fp)


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

