# Newer version of richness_by_zones notebooks
# Compare richness differences by zones

# Load libraries ----
source('R/initialize.R')

buffer_distance = 10
rich_code = 'Likhd'
samp_n = 350


# Zones raster and lookup table
zones_fn <- 'stack_anps3_biomes7_usv3'
zone_polys_dir <- file.path('data', 'data_out', 'zone_polys')
zones_ras_fp <- file.path(zone_polys_dir, str_c(zones_fn, '.tif'))
zones_lu_fp <- file.path(zone_polys_dir, str_c(zones_fn, '_lu.rds'))

if(!file.exists(zones_ras_fp)){
  
  # File paths
  anps_fp <- file.path(zone_polys_dir, 'ANPs_buff10km.gpkg')
  biom_diss_fp <- file.path(zone_polys_dir, 'ecoregions_diss7.gpkg')
  usv_fn <- 'usv250s6gw_diss_3class.gpkg'
  usv_fp <- file.path(zone_polys_dir, usv_fn)
  
  # Input pollinator richness 
  pol_group <- 'Abejas'
  fp_tail <- file.path('sdm', str_c(unq_code, '_', dfilt_code), rf_code, pol_group)
  pred_dir <- here::here('data', 'data_out', fp_tail)
  rich_fn <- str_glue('{params$rich_code}_rich_{pol_group}_{dfilt_code}_{params$nspecies}species')
  rich_tif_fp <- here::here(pred_dir, str_glue('{rich_fn}.tif'))
  rich_ras <- terra::rast(rich_tif_fp)
  
  # ANP polys
  anps_ras_fp <- str_c(tools::file_path_sans_ext(anps_fp), '.tif')
  if(!file.exists(anps_ras_fp)){
    # ANP polys
    anps_polys <- terra::vect(anps_fp) 
    anps_diss <- terra::aggregate(anps_polys, by='zone')
    
    # Get lookup table
    anp_lu <- terra::values(anps_diss) %>% rownames_to_column('id') %>% 
      select(id, name=zone)
    anp_lu <- setNames(anp_lu$name, anp_lu$id)
    
    anp_lu
    
    # Rasterize and save
    anps_ras <- terra::rasterize(anps_diss, rich_ras)
    anps_ras %>% terra::writeRaster(anps_ras_fp, overwrite=T, 
                                    options=c("dstnodata=-99999"),
                                    wopt=list(gdal='COMPRESS=LZW'))
  } else {
    anps_ras <- terra::rast(anps_ras_fp)
  }
  
  # Biome polys
  biom_ras_fp <- str_c(tools::file_path_sans_ext(biom_diss_fp), '.tif')
  if(!file.exists(biom_ras_fp)){
    # Load polygons and dissolve
    biom_polys <- terra::vect(biom_diss_fp) 
    biom_diss <- terra::aggregate(biom_polys, by='DESECON1')
    
    # Get lookup table
    biom_lu <- terra::values(biom_diss) %>% rownames_to_column('id') %>% 
      select(id, name=DESECON1)
    biom_lu <- setNames(biom_lu$name, biom_lu$id)
    
    # Rasterize and save
    biom_ras <- terra::rasterize(biom_diss, rich_ras)
    biom_ras %>% terra::writeRaster(biom_ras_fp, overwrite=T, 
                                    options=c("dstnodata=-99999"),
                                    wopt=list(gdal='COMPRESS=LZW'))
  } else {
    biom_ras <- terra::rast(biom_ras_fp)
  }
  
  # Land use polys
  usv_ras_fp <- str_c(tools::file_path_sans_ext(usv_fp), '.tif')
  if(!file.exists(usv_ras_fp)){
    # Load polygons
    usv_polys <- terra::vect(usv_fp)
    
    # Get lookup table
    usv_lu <- terra::values(usv_polys) %>% 
      rownames_to_column('id') %>% 
      select(id, name=tipo)
    usv_lu <- setNames(usv_lu$name, usv_lu$id)
    
    # Rasterize and save
    usv_ras <- terra::rasterize(usv_polys, rich_ras)
    usv_ras %>% terra::writeRaster(usv_ras_fp, 
                                   overwrite=T, options=c("dstnodata=-99999"), 
                                   wopt=list(gdal='COMPRESS=LZW'))
  } else {
    usv_ras <- terra::rast(usv_ras_fp)
  }
  
  # Combine zone rasters
  combo_ras <- c(anps_ras, biom_ras, usv_ras)
  names(combo_ras) <- c('anp_zone', 'biome', 'usv')
  
  combo_ras %>% terra::writeRaster(zones_ras_fp, 
                                   overwrite=T, options=c("dstnodata=-99999"), 
                                   wopt=list(gdal='COMPRESS=LZW'))
  
  saveRDS(list(anp=anp_lu, biom=biom_lu, usv=usv_lu), file=zones_lu_fp)
} 

# Directory for output richness dataframes
# div_dir <- file.path('data', 'data_out', 'diversity_by_zones', 'from_richness', zones_fn)
div_dir <- file.path('data', 'data_out', 'diversity_by_zones', rf_name, zones_fn)
dir.create(div_dir, recursive=TRUE, showWarnings = F)

# Richness data frames ----
# Convert rasters to dataframes of richness by ANP zone, biome, land cover, and pollinator group. 
# Load stacked zones 
combo_ras <- terra::rast(zones_ras_fp)
lu <- readRDS(zones_lu_fp)

# Sum likelihood maps ----
# Get pollinators dataframe - Quesada updated, top 10 crops
lu_df <- 'data/tidy/Quesada/crop_pllntrs_updated_top10.csv'
lu_df <- 'data/tidy/Quesada/crop_pllntrs_from_ap2_updated_tidied.csv'
pol_df <- read_csv(lu_df)

# See values in subgroup_a (equivalent to sociabilidad for Hymenoptera)
pol_df %>% filter(order == order) %>% distinct(subgroup_a)

# Get Millard pollinators list -----
# HERE NOW -----
# 2019 output
pllntrs_millard_csv <- 'data/input_data/Millard/05. genus_aggregations.csv'
pllntrs_millard <- read_csv(pllntrs_millard_csv)
pllntrs_millard %>% 
  filter(str_detect(unique_name, regex('Mexico', ignore_case = TRUE)))

# 2021 supp data 1 
# KEPT = genus initially identified as a pollinator and confirmed by panel of experts
pllntrs_millard_csv <- 'data/input_data/Millard/41467_2021_23228_MOESM4_ESM.csv'
pllntrs_millard <- read_csv(pllntrs_millard_csv)
pllntrs_millard <- pllntrs_millard %>% 
  filter(str_detect(status, regex('KEPT', ignore_case = TRUE)))

# subset of pollinating species in the PREDICTS database
# pllntrs_millard_rds <- 'data/input_data/Millard/PREDICTS_pollinators_8_exp.rds'
# pllntrs_predicts <- readRDS(pllntrs_millard_rds)

# Compare Millard pollinators to our list
anti_join(pllntrs_millard, pol_df, by = c(Genus = 'genus'))
semi_join(pllntrs_millard, pol_df, by = c(Genus = 'genus'))

genus_list_q <- pol_df %>% distinct(genus)
genus_list_m <- pllntrs_millard %>% distinct(Genus)

# Genera in our list that aren't found in Millard's - 122
anti_join(genus_list_q, genus_list_m, by = c(genus = 'Genus'))
# Genera common to both our list and Millard's - 71
semi_join(genus_list_q, genus_list_m, by = c(genus = 'Genus'))

pols_distinct <- pol_df %>% select(species:subgenus) %>% distinct
genus_list_m <- pllntrs_millard %>% distinct(Genus)

# Genera common to both our list and Millard's - 71
pols_filt <- semi_join(pols_distinct, genus_list_m, by = c(genus = 'Genus'))
pols_filtx <- anti_join(pols_distinct, genus_list_m, by = c(genus = 'Genus'))
t1 <- pols_distinct %>% count(order) %>% rename(total = 'n')
t2 <- pols_filt %>% count(order) %>% rename(in_Millard = 'n')
t3 <- pols_filtx %>% count(order) %>% rename(not_Millard = 'n')

sum_tbl <- plyr::join_all(list(t1, t2, t3))


# list species (such as nocturnal butterflies)
if(order == 'Lepidoptera'){
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
  
  rich_fn <- str_glue('richness_{name}_{rf_name}_{length(sp_fps)}species_lkhd')
  rich_plot_fp <- file.path(rf_fig_dir, str_glue('{rich_fn}.png'))
  rich_tif_fp <- file.path(pred_dir, str_glue('{rich_fn}.tif'))
  
  stack_sdms(sp_fps, rich_tif_fp, rich_plot_fp, mex)
  
}

# ALT: List richness tifs ----
tif_query <- str_glue(
  '{rich_code}_rich_.*_{dfilt_code}_{nspecies}species.tif')

apis_code <- ifelse(order == 'Hymenoptera', ifelse(noApis, '_noApis', ''), '')
tif_query <- str_glue(
  'richness_.*_{dfilt_code}_.*_{nspecies}species{apis_code}_lkhd.tif')
(fps <- list.files(pred_dir, pattern=tif_query, recursive=T, full.names=T))

for(rich_tif_fp in fps){
  
  # Resulting path for CSV
  rich_fn <- basename(tools::file_path_sans_ext(rich_tif_fp))
  div_df_fp <- file.path(div_dir, str_c(rich_fn, '.csv'))
  
  if(!file.exists(div_df_fp)){
    
    # Load raster and add to SpatRaster
    rich_ras <- terra::rast(rich_tif_fp)
    names(rich_ras) <- 'richness'
    combo_all <- c(combo_ras, rich_ras)
    
    # Convert to DF and recode
    div_df <- terra::as.data.frame(combo_all, cells=TRUE) %>% 
      mutate(pol_group = name)
    
    # Save
    div_df %>% write_csv(div_df_fp)
    
  }
}

# Combine all richness DFs into one ----
# Load DFs and row bind
fps <- list.files(div_dir, pattern='csv$', full.name=T) %>% 
  str_subset(str_glue('/richness_.*{dfilt_code}'))
div_df1 <- fps %>% purrr::map_dfr(read_csv, col_types='ifffdf')

# Recode biome names
biom_codes_all <- tibble(id=names(lu$biom), 
                         name=lu$biom, 
                         code=c('CalMedit', 'Desierto', 'ElevSemiar', 'GranPlanic', 
                                'SelvCalHum', 'SelvCalSec', 'SierTemp'))
biom_codes <- setNames(biom_codes_all$code, biom_codes_all$id)

# Convert numeric code to factors
div_df <- div_df1 %>% 
  mutate(anp_zone = recode(anp_zone, !!!lu$anp) %>% 
           factor(levels = c('Outside buffer', 'Buffer 10 km', 'Inside NPA')),
         biome = recode(biome, !!!biom_codes) %>% 
           factor(),
         usv = recode(usv, !!!lu$usv) %>% 
           factor(levels = c('veg', 'ag', 'otro')))

# get_grp_diffs evaluates the statistical significance of groups and returns 
# medians, CIs, and significant groups for plotting. 
# It will be applied to each group of pollinator, biome, and either ANP zone or land cover type.
get_grp_diffs <- function(df, grp_var){
  
  # Get medians and confidence intervals
  f <- str_c('richness ~ ', grp_var)
  Sum <- do.call("groupwiseMedian", 
                 list(as.formula(f), 
                      data = as.name('df'),
                      conf       = 0.95,
                      R          = 5000,
                      percentile = TRUE,
                      bca        = FALSE,
                      digits     = 3))
  
  if(length(unique(df[[grp_var]])) == 1) {
    
    print('Only one group.')
    cld <- tibble(zone = unique(df[[grp_var]]), sig_group = 'a')
    
  } else {
    
    # Check whether any groups are significantly different
    kw <- kruskal.test(df[['richness']] ~ df[[grp_var]]) 
    if(kw$p.value > 0.05) print('No significantly different groups')
    
    # Post-hoc test to compare groups
    PT <- pairwise.wilcox.test(df[['richness']], df[[grp_var]], 
                               p.adjust.method = 'fdr')
    PT1 = rcompanion::fullPTable(PT$p.value)
    cld_list <- multcompView::multcompLetters(PT1,
                                              compare="<",
                                              threshold=0.05,  # p-value to use as significance threshold
                                              Letters=letters,
                                              reversed = FALSE)
    lv <- cld_list$Letters
    cld <- tibble(sig_group = factor(lv, levels = sort(unique(lv))), 
                  zone = factor(names(lv)))
    
  }
  
  # Join letter indicating significantly different groups to medians and CIs
  sum1 <- select(tibble(cld), zone, sig_group) %>% 
    left_join(tibble(Sum), 
              by=c(zone=grp_var))
  
  return(sum1)
}

plot_pollin_facet <- function(df, pol_grp, var2='usv', grp_var='anp_zone', 
                              patch.position='inner', legend_on=FALSE){
  
  xmax <- max(df$Percentile.upper)
  
  df_rect <- distinct(df, biome, .data[[var2]], pol_group)
  if(grp_var == 'anp_zone'){
    
    var_levels <- c( 'Outside buffer', 'Buffer 10 km', 'Inside NPA')
    facet_y <- 'Land cover'
    
  } else if (grp_var == 'usv') {
    
    var_levels <- c( 'otro', 'ag', 'veg')
    facet_y <- 'ANP zone'
    
  }
  
  if(var2 != 'pol_group'){
    p <- ggplot(df) +
      
      # shade facet by ANP zone
      geom_rect(data=df_rect, aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
                                  fill=.data[[var2]]), alpha=0.25, color='white') +
      scale_fill_grey(name = facet_y, start=0.05, end=0.75) 
    
  } else {
    
    p <- ggplot(df)
    
  }
  
  p <- p +
    
    # plot CIs
    geom_errorbar(aes(y = zone, x = Median, color=sig_group,
                      xmin = Percentile.lower, xmax = Percentile.upper),
                  width = 0.3, size  = 0.5) +
    
    # plot medians
    geom_point(aes(y = zone, x = Median, color=sig_group), 
               shape = 15, size  = 2) +
    xlab("Median richness with 95% CI") +
    
    scale_color_manual(values=c('a'='#F8766D', 'ab'='#C77CFF', 'b'='#00BFC4', 
                                'bc'='#7CAE00', 'c'='#DE8C00'), 
                       # limits = c('a', 'ab', 'b', 'bc', 'c'),
                       name='Group',
                       drop=FALSE) +
    
    scale_y_discrete(limits = var_levels) +
    scale_x_continuous(n.breaks= 3) +
    
    # Facet by x=pollinator group, y = ANP zone
    facet_grid(vars(.data[[var2]]), vars(biome)) +
    ggtitle(pol_grp) +
    theme_minimal() +
    theme(
      plot.title = element_text(size=rel(1), face='bold', hjust=0),
      panel.spacing = unit(0.5, 'mm'),
      panel.border=element_blank(),
      # panel.grid = element_line(color='darkgray'), 
      # panel.grid.minor.x = element_line(color='darkgray'),
      panel.grid.major.y=element_blank(),
      panel.ontop = FALSE,
      panel.background = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(), 
      axis.ticks.y = element_blank(),
      # axis.text.x = element_text(margin = margin(t = .2, unit = "mm")),
      strip.text = element_blank(),
      strip.background = element_blank(),
      strip.switch.pad.grid = unit(0.5, 'mm'),
      legend.position = 'none'
    )
  
  return(p)
  
}

# Sample DF ----
# Stratified random sample
fn <- str_c('samp', samp_n, '_', rich_code, '_', 
            dfilt_code, '_', nspecies, 'species_',
            length(levels(div_df$pol_group)), 'pols', apis_code, '.rds')
samp_fp <- file.path(div_dir, fn)

if(!file.exists(samp_fp)){
  # Sample the same cells for all the pollinator groups
  
  # Get all usable cells
  grp1 <- div_df %>% distinct(pol_group) %>% slice(1) %>% deframe
  samp_idx <- div_df %>% 
    filter(pol_group == as.character(grp1)) %>% 
    group_by(anp_zone, biome, usv) %>%
    filter(n() > samp_n)
  
  # Random sample stratified by group
  set.seed(1)
  samp_idx <- samp_idx %>% 
    do(slice_sample(., n = samp_n)) %>%
    ungroup
  
  df_samp <- div_df %>% filter(cell %in% samp_idx$cell)
  
  # Save
  df_samp %>% saveRDS(samp_fp)
  
} else {
  
  df_samp <- readRDS(samp_fp)
  
}

# Compare groups ----
# no USV - by ANP zone w/o land cover ----
grp_var <- 'anp_zone'
var2 <- NA

# File name
fn <- str_c('grpdiffs_samp', samp_n, '_', 
            grp_var, '_', 'noUSV_',
            rich_code, '_', 
            dfilt_code, '_', nspecies, 'species_',
            length(levels(div_df$pol_group)), 'pols',
            apis_code, '.rds')
grp_diffs_fp <- file.path(div_dir, fn)

if(!file.exists(grp_diffs_fp)){
  
  # Get group differences for each combination of biome and ANP zone
  if(is.na(var2)){
    df_nest <- df_samp %>% group_by(pol_group, biome) %>% nest()
  } else {
    df_nest <- df_samp %>% group_by(pol_group, biome, .data[[var2]]) %>% nest()
  }
  medians_anp <- df_nest %>%
    mutate(sum = purrr::map(data, get_grp_diffs, grp_var)) %>%
    select(-data) %>% 
    unnest(cols=sum)
  
  sg_levels <- sort(as.character(unique(medians_anp$sig_group)))
  medians_anp <- medians_anp %>% mutate(sig_group = factor(sig_group, 
                                                           levels = sg_levels))
  
  # Save 
  medians_anp %>% saveRDS(grp_diffs_fp)
  
} else {
  
  medians_anp <- readRDS(grp_diffs_fp)
  
}

# Plot
medians_anp %>% 
  plot_pollin_facet('All', var2='pol_group', grp_var=grp_var) +
  scale_x_continuous(n.breaks= 4) +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_text(margin= margin(t=6, unit='pt')),
        # axis.text.x = element_text(hjust=c(0, 0.5, 1)), 
        plot.title = element_blank(), 
        legend.position = 'bottom', 
        legend.box.spacing = unit(1, 'pt'),
        strip.text = element_text(size=rel(0.8)),
        panel.background = element_rect(fill=alpha('#A6A6A6', 0.3), 
                                        color='white'),
        panel.grid = element_line(color='lightgray'),
        panel.grid.minor.x = element_line(color='lightgray')
  ) +
  labs(caption = "Different colors within one facet indicate statistically significant (p < 0.05) difference between groups.")

# Save figure
# fig_dir <- here::here('figures', 'diversity_by_zones')
# fig_fp <- file.path(fig_dir, dfilt_code,
#                     str_c('grpdiffs_samp', samp_n, '_', 
#                           grp_var, '_', 'noUSV_',
#                           rich_code, '_', 
#                           dfilt_code, '_', 
#                           nspecies, 'species_',
#                           length(levels(div_df$pol_group)), 'pols', 
#                           apis_code, '.png'))
fig_fp <- file.path(div_dir,
                    str_c('grpdiffs_samp', samp_n, '_', 
                          grp_var, '_', 'noUSV_',
                          rich_code, '_', 
                          dfilt_code, '_', 
                          nspecies, 'species_',
                          length(levels(div_df$pol_group)), 'pols',
                          apis_code, '.png'))
if(!file.exists(fig_fp)) ggsave(fig_fp, width=9, height=4)
ggsave(fig_fp, width=7, height=3)
