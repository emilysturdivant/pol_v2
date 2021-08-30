# 
# Create polygons to delimit zones (ANP: inside, buffer, outside; biomes; land use)
# Create diversity table from those zones.
# Script created after prep_create_ANP... and explore_subset_by_biome

# Load libraries ---------------------------------------------------------------
# library(tidyverse)
# library(sf)
# library(units)
# library(tmap)
# tmap_mode('view')

# Initialize -------------------------------------------------------------------
source('R/initialize.R')

# Functions --------------------------------------------------------------------
get_diversity_for_zone_polys <- function(df, zone_polys, zone_str=NA, polys_id_fld, 
                                         polys_area_fld='area_ha'){
  
  df <- df %>% 
    # Group pollinator points by intersecting ANP 
    st_transform(crs = st_crs(zone_polys)) %>% 
    st_join(zone_polys, left=F) %>% {
      if(nrow(.) > 0) {
        
        # Get diversity and abundance of species
        get_diversity_metrics(., polys_id_fld, polys_area_fld)
        
      } else {
        
        # Make empty tibble if there are no pollinators in the given zone
        tibble(!!polys_id_fld := character(), !!polys_area_fld := numeric())
        
      } 
    }
  
  # Get just the polygon ID and area
  all_ids <- zone_polys %>% 
    st_drop_geometry() %>% 
    select(.data[[polys_id_fld]], .data[[polys_area_fld]])
  
  # Join so that we have NA anywhere statistics are missing
  df <- full_join(all_ids, df)
  
  if(!is.na(zone_str)){
    
    df <- mutate(df, zone = zone_str)
    
  }
  
  # Return
  return(df)
  
}

get_pollinators_in_anp_zones <- function(name, date_range, anps_terr, 
                                         buffer_distance, 
                                         polys_id_fld='rowname', 
                                         polys_area_fld = 'area_ha',
                                         pol_dir = 'data/data_out/pollinator_points'){
  
  # Zone codes
  inside_str <- str_c('Inside NPA')
  buff_str <- str_c('Buffer ', buffer_distance, ' km')
  outside_str <- str_c('Outside buffer')
  
  # Convert date range
  date_min <- as.POSIXct(str_c(date_range[[1]], "-01-01"))
  date_max <- as.POSIXct(str_c(date_range[[2]], "-12-31"))
  
  # Load pollinator file
  pol_fp <- file.path(pol_dir, str_c(name, '.geojson'))
  df <- st_read(pol_fp) %>% 
    # Filter to date range
    filter(eventDate >= date_min & eventDate <= date_max) 
  
  # Get diversity and abundance of pollinators within NPA ----
  anps_inside <- get_diversity_for_zone_polys(df, anps_terr, inside_str, 
                                              polys_id_fld, polys_area_fld)
  
  # NPA buffer ----
  anps_buff <- get_buffer_conditional(buffer_distance, anps_terr, area_fld='area_ha')
  
  # Get diversity and abundance of pollinators in NPA buffer
  anps_buffer <- get_diversity_for_zone_polys(df, anps_buff, buff_str, 
                                              polys_id_fld, polys_area_fld)
  
  # Voronoi polygons ----
  vpols <- get_voronois_conditional(anps_terr, buffer_distance)
  
  # Get diversity and abundance of pollinators outside NPA buffer
  anps_outside <- get_diversity_for_zone_polys(df, vpols, outside_str, 
                                               polys_id_fld, polys_area_fld)
  
  # Concatenate the three DFs
  anps_df <- bind_rows(anps_inside, anps_buffer, anps_outside)
  
  # Set pollinator group name
  anps_df$pol_group <- name
  
  # Save
  fp_out <- file.path('data/data_out/ANPs', 
                      str_c('anps_', name, '_', strftime(date_min, format="%Y"), 
                            '_to_', strftime(date_max, format="%Y"), '_buffer', 
                            buffer_distance, 'km.csv'))
  anps_df %>% write_csv(fp_out)
  
  # Return
  return(anps_df)
}

# Newer functions ----
#' @export
clean_anps_conditional <- function(anp_fp, anp_terr_fp, crs, mex, overwrite = FALSE) {
  
  if(file.exists(anp_terr_fp) & !overwrite) return(st_read(anp_terr_fp))

  # Load polygons, clean, and convert to single-part
  anps_proj <- st_read(anp_fp) %>% 
    select(ID_ANP) %>% 
    st_transform(crs=crs) %>% 
    st_set_precision(1e5)  %>% 
    st_make_valid %>% 
    st_collection_extract('POLYGON') %>% 
    st_simplify(dTolerance=40, preserveTopology=T) %>% 
    st_cast("MULTIPOLYGON") %>% 
    st_cast("POLYGON") %>%
    rownames_to_column %>% 
    st_make_valid()
  
  # Clip to terrestrial portions using Mexico boundary
  anps_terr <- anps_proj %>% 
    st_intersection(mex)
  
  # Get area
  anps_terr$area_ha <- anps_terr %>% 
    st_area %>% units::set_units('ha') %>% units::set_units(NULL)
  
  # Save
  anps_terr %>% st_write(anp_terr_fp, delete_dsn=T)
  
  return(anps_terr)
}

#' @export
get_buffer_conditional <- function(buffer_distance, anps_proj, area_fld='buff_area_ha'){
  # NPA buffer 
  buffer_fp <- str_c('data/tidy/analysis_zones/ANPs/anps_buffer_', buffer_distance,'km.gpkg')
  
  if(file.exists(buffer_fp)) {
    
    # Load buffer
    anps_buff <- st_read(buffer_fp)
    return(anps_buff)
    
  } 
    
  print(str_glue('Creating ANP buffer of {buffer_distance} km...'))
  
  # Create buffer around each ANP part
  anps_buff <- anps_proj %>% 
    st_buffer(buffer_distance, validate=T) %>% 
    st_difference(st_union(anps_proj))
  
  # Buffer 
  anps_buff2 <- anps_buff %>% 
    group_by(rowname, ID_ANP) %>% 
    summarise() %>% 
    ungroup() %>%
    st_cast('MULTIPOLYGON')
  
  # Get area
  anps_buff2[[area_fld]] <- anps_buff2 %>% 
    st_area %>% 
    units::set_units('ha') %>% 
    units::set_units(NULL)
  
  # Save
  anps_buff2 %>% st_write(buffer_fp, delete_dsn=T)
  
  return(anps_buff2)
}

#' @export
make_voronois <- function(zone_polys, region_poly){
  # Make Voronoi polygons from ANP polygons 
  
  zone_polys <- zone_polys %>% 
    st_transform(crs=st_crs(region_poly)) %>% 
    ungroup()
  
  # Dissolve polys into one
  zone_union <- zone_polys %>% 
    group_by(ID_ANP) %>% 
    summarize() %>% 
    st_filter(region_poly) %>% 
    st_simplify(dTolerance=40) %>% 
    st_union() 
  
  # Load polygons, clean, and convert to single-part
  zone_union <- zone_union %>% 
    st_buffer(units::set_units(0.1, 'km')) %>% 
    st_buffer(-units::set_units(0.1, 'km'))
  
  # Voronoi polygons
  vpols1 <- zone_union %>% 
    st_voronoi %>% 
    st_cast %>% 
    st_sf
  
  # join feature IDs (rowname) from ANPs
  vpols2 <- vpols1 %>% 
    st_join(zone_polys, join=st_nearest_feature) %>% 
    st_make_valid()
  
  # Dissolve by IDs
  vpols3 <- vpols2 %>% 
    group_by(rowname) %>% 
    summarize %>% 
    ungroup
  
  # Clip out ANPs - processing is faster with difference as the last step
  vpols4 <- vpols3 %>% 
    st_difference(zone_union)
  
  # Clip to land
  vpols5 <- vpols4 %>% 
    st_intersection(region_poly) %>% 
    st_collection_extract('POLYGON') %>% 
    group_by(rowname) %>% 
    summarize %>% 
    ungroup
  
  # Get area
  vpols <- vpols5 %>% 
    mutate(area_ha = st_area(geometry) %>% 
             units::set_units('ha') %>% units::set_units(NULL))
    
  return(vpols)
}

#' @export
get_voronois_conditional <- function(anps_proj, buffer_distance, vrnoi_fp=NA, region_poly){
  # Save
  if(is.na(vrnoi_fp)){
    
    vrnoi_fp <- str_c('data/tidy/analysis_zones/ANPs/voronoi_outside_anps_and_buffer_', 
                      buffer_distance,'km.gpkg')
    
  }

  if(file.exists(vrnoi_fp)){
    
    # Load file if it already exists 
    vpols <- st_read(vrnoi_fp)
    return(vpols)
    
  }
    
  # Create buffer polygons
  anps_buff <- anps_proj %>% 
    st_buffer(buffer_distance, validate=T) 
  
  # Make polygons
  vpols <- make_voronois(zone_polys=anps_buff, region_poly)
  
  # Save
  vpols %>% st_write(vrnoi_fp, delete_dsn=T)
  
  return(vpols)
}

#' @export
get_biomes_conditional <- function(biom_diss_fp, biom_dir, mex) {
  
  if(file.exists(biom_diss_fp)) {
    
    biom_diss <- st_read(biom_diss_fp)
    return(biom_diss)
    
  }
  
  # Filepath
  shp_fp <- list.files(biom_dir, '.shp$', full.names = T, recursive = T)
  
  # Load biomes shapefile
  biom <- st_read(shp_fp) %>% 
    st_transform(st_crs(mex)) %>% 
    st_simplify(dTolerance=40)
  
  # Dissolve to ecorregiones (7 in Mexico)
  biom_diss <- biom %>% 
    group_by(DESECON1) %>% 
    summarise() %>% 
    st_cast("POLYGON") %>% 
    st_simplify(dTolerance=40)
  
  # Save
  biom_diss %>% st_write(biom_diss_fp)
  
  return(biom_diss)
}

#' @export
merge_anp_zones_conditional <- function(bind_fp, anps_terr, anps_buff, vpols, buffer_distance) {
  
  if(file.exists(bind_fp)){
    
    anp_zones <- st_read(bind_fp)
    return(anp_zones)
    
  }
  
  # Zone codes
  inside_str <- str_c('Inside NPA')
  buff_str <- str_c('Buffer ', buffer_distance, ' km')
  outside_str <- str_c('Outside buffer')
  
  # Set zone names before binding
  anps_terr$zone <- inside_str
  anps_buff$zone <- buff_str
  vpols$zone <- outside_str
  
  # Create lookup list for recode
  lu <- setNames(c('anp', 'buff', 'out'), 
                 c(inside_str, buff_str, outside_str)) 
  
  # Bind SFCs
  anp_zones <- bind_rows(st_cast(anps_terr, 'MULTIPOLYGON'),
                         st_cast(anps_buff, 'MULTIPOLYGON'),
                         st_cast(vpols, 'MULTIPOLYGON')
                         ) %>% 
    mutate(zone_temp = recode(zone, !!!lu),
           name = str_c(zone_temp, rowname, sep='_')) %>% 
    select(name, ID_ANP, area_ha, zone)
  
  # Save
  anp_zones %>% st_write(bind_fp, delete_dsn=T)
  
  return(anp_zones)
}

get_zone_diversity <- function(pol_group, date_range, zone_polys, zone_id_fld) {
  # Dates
  date_min <- as.POSIXct(str_c(date_range[[1]], "-01-01"))
  date_max <- as.POSIXct(str_c(date_range[[2]], "-12-31"))
  
  # Load pollinator file
  pol_dir <- 'data/data_out/pollinator_points/with_duplicates'
  pol_fp <- file.path(pol_dir, str_c(pol_group, '.gpkg'))
  pol_df1 <- st_read(pol_fp) %>% 
    # Filter to date range
    filter(eventDate >= date_min & eventDate <= date_max) 
  
  zone_div <- pol_df1 %>%
    get_diversity_for_zone_polys(zone_polys, polys_id_fld=zone_id_fld) %>%
    left_join(st_drop_geometry(zone_polys))
  
  zone_div$pol_group <- pol_group
  
  return(zone_div)
}

# Load data --------------------------------------------------------------------
# # Biomes
# shp_fp <- 'data/input_data/environment_variables/CONABIO/ecort08gw.shp'
# biomes <- st_read(shp_fp)

# Mexico
mex <- st_read(mex_fp) %>% 
  st_transform(crs=crs) %>% 
  st_simplify(dTolerance=40, preserveTopology=T) 

mex0 <- st_union(mex)

# ANP polygons ----
anps_terr <- clean_anps_conditional(anps_in_fp, anp_terr_fp, crs, mex) %>% 
  st_transform(crs=crs)
anps_terr <- clean_anps_conditional(anps_in_fp, anp_terr_fp, crs, mex0) %>% 
  st_transform(crs=crs)

# ANP buffer 
anps_buff <- get_buffer_conditional(buffer_distance, 
                                    anps_proj = anps_terr, 
                                    area_fld = 'area_ha') %>% 
  st_transform(crs=crs)

# Voronoi polygons
vpols <- get_voronois_conditional(anps_terr, buffer_distance, region_poly=mex0) %>% 
  st_transform(crs=crs)

# Get merged ANP zones ----
anp_zones <- merge_anp_zones_conditional(bind_fp, anps_terr, anps_buff, 
                                           vpols, buffer_distance)

# Ecoregions --------------------------------------------------------
# Biomes CONABIO 
biom_diss <- get_biomes_conditional(biom_diss_fp, 
                                    biom_dir='data/input_data/environment_variables/CONABIO', 
                                    mex)

# # Land cover ----
# # Land cover from INEGI via CONABIO
# fp_usv <- "data/input_data/ag_INEGI_2017/other_sources/usv250s6gw.shp"
# 
# # Output files
# # grp_name
# usv_fname <- tools::file_path_sans_ext(basename(fp_usv))
# # crop_class_fp <- file.path("data/intermediate_data/ag_by_region", 
# #                            str_c(str_c(usv_fname, 'tipo', '7classes', 
# #                                        str_c(grp_name, collapse=''), 
# #                                        sep='_'), '.gpkg'))
# # 
# 
# # Load file
# usv1 <- st_read(fp_usv) %>% 
#   st_make_valid() %>% 
#   st_transform(crs)
# 
# usv2 <- usv1 %>% 
#   # st_crop(sub1) %>%
#   st_simplify(dTolerance=40, preserveTopology=T)
# 
# usv2 %>% object.size() %>% print(units='MB')
# 
# # classify into 7 classes
# usv3 <- usv2 %>% 
#   mutate(tipo = case_when(
#     str_detect(DESCRIPCIO, 'INDUCIDO$') ~ 'inducido',
#     str_detect(DESCRIPCIO, '^PASTIZAL') ~ 'pastizal',
#     str_detect(CVE_UNION, '^V|^B|^S|^P') ~ 'vegetacion',# bosque, selva, pastizal, sabana, palmar, etc.
#     str_detect(DESCRIPCIO, '^AGRICULTURA') ~ 'agricultura', # does not include shifting cultivation (nómada)
#     str_detect(CVE_UNION, '^ACUI|^H2O') ~ 'agua',
#     str_detect(CVE_UNION, '^ADV|^DV') ~ 'sin_veg',
#     str_detect(CVE_UNION, '^AH') ~ 'construido',
#     TRUE ~ 'otro'
#   ))
# 
# usv %>% object.size() %>% print(units='MB')
# 
# # # Save
# # usv %>% st_write(crop_class_fp, delete_dsn=T)
# # usv <- sf::st_read(crop_class_fp)
# 
# # Reclass and dissolve
# key <- c(inducido='veg', pastizal='veg', vegetacion='veg', agricultura='ag', 
#          agua='otro', sin_veg='otro', construido='otro')
# usv_diss1 <- usv3 %>% 
#   mutate(tipo =  recode(tipo, !!!key)) %>% 
#   group_by(tipo) %>% 
#   summarize()
# 
# usv_diss2 <- usv_diss1 %>% 
#   st_simplify(dTolerance=40, preserveTopology=T)
# 
# usv_diss3 <- usv_diss2 %>% 
#   nngeo::st_remove_holes(100000) 
# 
# usv_diss2 %>% object.size() %>% print(units='MB')
# 
# # Save dissolved USV ----
# usv_fname <- tools::file_path_sans_ext(basename(fp_usv))
# diss_fp <- file.path("data/intermediate_data",
#                      str_c(str_c(usv_fname, 'diss', '3class',
#                                  sep='_'), '.gpkg'))
# usv_diss2 %>% st_write(diss_fp, delete_dsn=T)
# 
# # Load dissolved land cover (3 classes) ----
# usv_diss <- st_read(diss_fp)
# 
# # Get shifting agriculture ----
# # nma_fp <- file.path("data/intermediate_data/ag_by_region", 'polys_nma_ChiapasTabasco.gpkg')
# # fn <- str_c(str_c(tools::file_path_sans_ext(basename(fp_usv)), 
# #                   'diss', '4class', 'nma', str_c(grp_name, collapse=''), sep='_'), '.gpkg')
# # diss_nma_fp <- file.path("data/intermediate_data/ag_by_region", fn)
# # 
# # # Get all agriculture subset to region
# # fn <- str_c(str_c('polys_ag_INEGI', str_c(grp_name, collapse=''), sep='_'), '.gpkg')
# # fp <- file.path("data/intermediate_data/ag_by_region", fn)
# # 
# # if(file.exists(fp)) {
# #   
# #   polys <- st_read(fp)
# #   
# # } else {
# #   
# #   polys_fp <- file.path("data/intermediate_data", 'polys_ag_INEGI.gpkg')
# #   polys <- st_read(polys_fp) %>% 
# #     st_filter(sub1)
# #   
# #   polys %>% st_write(fp)
# #   
# # }
# # 
# # # Dissolve shifting agriculture polygons
# # nma_diss <- polys %>% 
# #   filter(CLAVE == 'NMA') %>% 
# #   st_crop(sub1) %>% 
# #   group_by() %>% 
# #   summarize() %>% 
# #   transmute(tipo='ag_nma')
# # 
# # # Save
# # nma_diss %>% st_write(nma_fp, delete_dsn=T)
# # 
# # # Load
# # nma_diss <- st_read(nma_fp)
# # 
# # # Add shifting ag to land cover
# # diss_wo_nma <- st_difference(usv_diss, nma_diss)
# # usv_nma <- bind_rows(diss_wo_nma, nma_diss)
# # 
# # # Save
# # usv_nma %>% st_write(diss_nma_fp)
# 
# 
# # Biome subset ================================================
# # Select biome
# biom_list <- biom_diss %>% st_drop_geometry %>% distinct %>% deframe
# (biom_name <- biom_list[[7]])
# biom1 <- biom_diss %>% 
#   filter(DESECON1==biom_name)
# 
# # Intersect ANP zones with biome subset
# anp_zones_biom1 <- st_intersection(anp_zones, biom1) %>% 
#   st_make_valid() %>% 
#   st_cast('MULTIPOLYGON') %>% 
#   st_cast('POLYGON')
# anp_zones_biom1$area_ha <- st_area(anp_zones_biom1) %>% set_units('ha') %>% set_units(NULL)
# 
# 
# # Intersect ANP zones with biome subset
# diss_biom_fp <- file.path("data/intermediate_data", 
#                           str_replace_all(biom_name, ' |-', ''),
#                           str_c(str_c(usv_fname, 'diss', '3class', 
#                                       str_replace_all(biom_name, ' |-', ''),
#                                       sep='_'), '.gpkg'))
# 
# if(!file.exists(diss_biom_fp)) {
#   usv_biom1 <- st_intersection(usv_diss, biom1) %>% 
#     st_simplify(dTolerance=40, preserveTopology=T) %>% 
#     nngeo::st_remove_holes(1000000) 
#   usv_biom1$area_ha <- st_area(usv_biom1) %>% set_units('ha') %>% set_units(NULL)
#   
#   # Save
#   dir.create(dirname(diss_biom_fp))
#   usv_biom1 %>% st_write(diss_biom_fp, delete_dsn=T)
# }
# 
# usv_biom1 <- st_read(diss_biom_fp) %>% 
#   st_make_valid() %>% 
#   st_cast('MULTIPOLYGON') %>% 
#   st_cast('POLYGON')
# 
# usv_biom1 %>% object.size() %>% print(units='MB')
# 
# # Look
# tm_shape(usv_biom1) + tm_polygons(alpha=.5)
# tm_shape(anp_zones_biom1) + tm_polygons(alpha=.5)
# 
# # Intersect USV with ANP zones ----
# biom_zone_fp <- file.path("data/intermediate_data", 
#                           str_replace_all(biom_name, ' |-', ''),
#                           str_c(str_c('ixn_ANPs_usv3_diss', 
#                                       str_replace_all(biom_name, ' |-', ''),
#                                       str_glue('buff{buffer_distance}km'),
#                                       sep='_'), '.gpkg'))
# 
# if(!file.exists(biom_zone_fp)) {
#   biom1_zones1 <- st_intersection(anp_zones_biom1, select(usv_biom1, tipo))
#   
#   # Zone codes
#   inside_str <- str_c('Inside NPA')
#   buff_str <- str_c('Buffer ', buffer_distance, ' km')
#   outside_str <- str_c('Outside buffer')
#   
#   lu <- setNames(c('anp', 'buff', 'out'), 
#                  c(inside_str, buff_str, outside_str)) 
#   
#   biom1_zones2 <- biom1_zones1 %>% 
#     st_collection_extract('POLYGON') %>%
#     st_cast('MULTIPOLYGON') %>% 
#     st_cast('POLYGON') %>% 
#     group_by(zone, tipo) %>% 
#     summarize() %>% 
#     mutate(zone = recode(zone, !!!lu),
#            name = str_c(zone, tipo, sep='_'))
#   
#   biom1_zones2$area_ha <- st_area(biom1_zones2) %>% set_units('ha') %>% set_units(NULL)
#   
#   # Save
#   biom1_zones2 %>% st_write(biom_zone_fp, delete_dsn=T)
#   
# }
# 
# # Load
# biom1_zones <- st_read(biom_zone_fp)
# 
# tm_shape(biom1_zones) + tm_polygons(alpha=.5)

# ~~~Pollinator points~~~ ----
# Diversity dataframe for ANP zone, biome, land use sites ----
# Date range
date_range <- c(2000, 2020)
pol_groups <- c('Abejas', 'Avispas', 'Colibries', 'Mariposas', 'Moscas', 'Murcielagos')

# Biome subset
biom_list <- biom_diss %>% st_drop_geometry %>% distinct %>% deframe
fps <- list.files(path="data/intermediate_data", 
  str_glue("ixn_ANPs_usv3_diss_.+_buff{buffer_distance}km.gpkg"), recursive=T, 
  full.names=T)

for(biom_zone_fp in fps) {
  biom_code <- basename(dirname(biom_zone_fp))
  biom_name <- agrep(biom_code, biom_list, value=T)
  
  zone_polys <- st_read(biom_zone_fp)
  zone_id_fld <- 'name'
  polys_area_fld <- 'area_ha'
  
  # Pollinators subset
  pol_div <- pol_groups %>% 
    map_dfr(get_zone_diversity, date_range, zone_polys, zone_id_fld) %>% 
    mutate(rich_by_1kkm = richness_norm*1000, 
           biome = biom_name)
  
  pol_div_fp <- file.path("data/data_out/diversity_by_zones", 
                          str_c(str_c('div_ANPs_usv3', biom_code, 
                                      str_glue('buff{buffer_distance}km'),
                                      sep='_'), '.gpkg'))
  pol_div %>% write_csv(pol_div_fp)
  
  fig_dir <- 'figures/anps_and_pollinator_exploration/con_biomas_y_usv'
  
  pol_div %>% 
    # filter(pol_group %in% c('Abejas', 'Colibries', 'Mariposas', 'Murcielagos')) %>% 
    ggplot(., aes(fill=tipo, x=zone, y=rich_by_1kkm)) +
    geom_bar(position='dodge', stat='identity') +
    scale_fill_manual(values=c("#E69F00", "#CC79A7", "#009E73"), "Uso de suelo") +
    facet_wrap(facets=vars(pol_group)) +
    ggtitle(biom_name)
  
  ggsave(file.path(fig_dir, str_c('bar_usv_', biom_code, '_buff',
                                  buffer_distance, 
                                  'km.png')), 
         width = 5.15, height=5.03)
}
# 
# # Get Jaccard distance - diversity and abundance of pollinators ----
# pol_df2 <- pol_df1 %>% 
#   # Group pollinator points by intersecting ANP 
#   st_transform(crs = st_crs(zone_polys)) %>% 
#   st_join(zone_polys, left=F)
# 
# # Group by ANP and pivot wide (separate column with count for each species)
# spec_df1 <- pol_df2 %>% 
#   st_drop_geometry %>% 
#   group_by(.data[[zone_id_fld]], .data[[polys_area_fld]], species) %>% 
#   summarize(n = length(species)) %>% 
#   ungroup 
# 
# spec_df2 <- spec_df1 %>% 
#   pivot_wider(id_cols = all_of(c(zone_id_fld, polys_area_fld)), 
#               names_from = species, 
#               values_from = n) %>%
#   replace(is.na(.), 0) %>% 
#   select(-area_ha) %>% 
#   left_join(st_drop_geometry(zone_polys)) %>% 
#   # mutate(name = str_c(abbreviate(DESECON1, use.classes=T), '_', rowname)) %>% 
#   column_to_rownames('name') %>%
#   select(-any_of(c('rowname', 'ID_ANP', polys_area_fld, 'DESECON1', 'zone', 'tipo')))
# 
# # Get distance matrix (Jaccard)
# df.jaccard <- vegan::vegdist(spec_df2, method="jaccard")
# 
# # Tree plot
# plot(hclust(df.jaccard),
#      hang = .2,
#      main = str_glue("{biom_name} sites \nclustered by Jaccard similarity"),
#      axes = FALSE,
#      ylab = "")
# 
# ggdendro::ggdendrogram(hclust(df.jaccard), rotate=T, size=2)
#  
# # Euclidean distance
# df.euclidean <- dist(spec_df2)
# 
# # Display using non-metric multidimensional scaling
# mdsE <- vegan::metaMDS(spec_df2, distance='euc')
# plot(mdsE, display="sites", type="text")
# 
# # Bray-Curtis
# mdsB <- vegan::metaMDS(spec_df2, distance="bray", autotransform=FALSE, trace=0)
# plot(mdsB, display="sites", type="text")

# Combine diversity DFs into one ----
fps <- list.files(path="data/data_out/diversity_by_zones",
                  pattern=str_c('div_ANPs_usv3_[^singlepart].+_buff', 
                                buffer_distance, 'km.csv'), 
                  full.names = T)
div_df <- map_dfr(fps, read_csv)

# Plot
div_df %>% 
  filter(pol_group %in% c('Abejas', 'Colibries', 'Mariposas', 'Murcielagos')) %>%
  ggplot(., aes(fill=tipo, x=zone, y=rich_by_1kkm)) +
  geom_bar(position='dodge', stat='identity') +
  scale_fill_manual(values=c("#E69F00", "#CC79A7", "#009E73"), "Uso de suelo") +
  facet_grid(vars(pol_group), vars(biome), 
             scales='free_y')

ggsave(file.path(fig_dir, str_c('bar_usv_ALL4_buff',
                                buffer_distance, 
                                'km.png')), 
       width = 9.15, height=7.03)

# Plot
div_df %>% 
  filter(pol_group %in% c('Abejas', 'Colibries', 'Mariposas', 'Murcielagos')) %>%
  ggplot(., aes(fill=zone, x=tipo, y=rich_by_1kkm)) +
  geom_bar(position='dodge', stat='identity') +
  scale_fill_manual(values=c("#E69F00", "#CC79A7", "#009E73"), "Zona") +
  facet_grid(vars(pol_group), vars(biome), 
             scales='free_y')

ggsave(file.path(fig_dir, str_c('bar_zona_ALL4_buff',
                                buffer_distance, 
                                'km.png')), 
       width = 9.15, height=7.03)

# Get diversity for undissolved polygons ----
biom_list <- biom_diss %>% st_drop_geometry %>% distinct %>% deframe
zone_id_fld <- 'name'
polys_area_fld <- 'area_ha'

for(biom_name in biom_list){
  biom_zone_fp <- file.path("data/intermediate_data", 
                            str_replace_all(biom_name, ' |-', ''),
                            str_c(str_c('ixn_ANPs_usv3_diss', 
                                        str_replace_all(biom_name, ' |-', ''),
                                        str_glue('buff{buffer_distance}km'),
                                        sep='_'), '.gpkg'))
  biom1_zones <- st_read(biom_zone_fp)
  
  # Biome subset
  zone_polys <- biom1_zones %>% 
    st_cast('POLYGON')
  zone_polys$area_ha <- st_area(zone_polys) %>% set_units('ha') %>% set_units(NULL)
  
  # Pollinators subset
  pol_div <- pol_groups %>% 
    map_dfr(get_zone_diversity, date_range, zone_polys, zone_id_fld) %>% 
    mutate(rich_by_1kkm = richness_norm*1000, 
           biome = biom_name)
  
  pol_div_fp <- file.path("data/data_out/diversity_by_zones", 
                          str_c(str_c('div_ANPs_usv3_singlepart', 
                                      str_replace_all(biom_name, ' |-', ''),
                                      str_glue('buff{buffer_distance}km'),
                                      sep='_'), '.csv'))
  pol_div %>% write_csv(pol_div_fp)
  
}

# Combine DFs for all biomes (singlepart) ----
fps <- list.files("data/data_out/diversity_by_zones", 
           pattern=str_glue('div_ANPs_usv3_singlepart.+buff{buffer_distance}km.gpkg'), 
           full.names = T)
div_df <- map_dfr(fps, read_csv)

div_single_fp <- file.path("data/data_out/diversity_by_zones", 
                           str_glue('div_ANPs_usv3_singlepart_buff{buffer_distance}km.csv'))
div_df %>% write_csv(div_single_fp)

# Multivariate regression ----
div_df <- read_csv(div_single_fp)

div_df1 <- div_df %>% filter(!is.na(richness_norm)) 

model <- lm(rich_by_1kkm ~ biome + zone + tipo, data=div_df1)
summary(model)
car::Anova(model)

# Individual pollinator groups
div_df2 <- div_df1 %>% 
  filter(pol_group == 'Murcielagos') 
model <- lm(rich_by_1kkm ~ biome + zone + tipo, data=div_df2)
summary(model)
car::Anova(model)

# Individual pollinator groups
div_df2 <- div_df1 %>% 
  filter(biome == 'Selvas Calido-Humedas',
         pol_group == 'Abejas') 
model <- lm(rich_by_1kkm ~ zone + tipo, data=div_df2)
summary(model)
car::Anova(model)

# Individual pollinator groups
div_df2 <- div_df1 %>% 
  filter(biome == 'Sierras Templadas',
         pol_group == 'Abejas',
         zone == 'buff') 
model <- lm(rich_by_1kkm ~ tipo, data=div_df2)
summary(model)
car::Anova(model)


