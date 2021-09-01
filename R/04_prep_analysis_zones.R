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
source('R/00_initialize.R')

# Functions --------------------------------------------------------------------
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

# Load data --------------------------------------------------------------------
# Mexico
mex <- st_read(mex_fp) %>% 
  st_transform(crs=crs) %>% 
  st_simplify(dTolerance=40, preserveTopology=T) 

mex0 <- st_union(mex)

# ANP polygons -----------------------------------------------------------------
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
                                    biom_dir = biom_dir, 
                                    mex)

# Land cover ----
# Land cover from INEGI via CONABIO
# Load file
usv1 <- st_read(fp_usv) %>%
  st_make_valid() %>%
  st_transform(crs) %>%
  st_simplify(dTolerance=40, preserveTopology=T)

# classify into 7 classes
usv7 <- usv1 %>%
  mutate(tipo = case_when(
    str_detect(DESCRIPCIO, 'INDUCIDO$') ~ 'inducido',
    str_detect(DESCRIPCIO, '^PASTIZAL') ~ 'pastizal',
    str_detect(CVE_UNION, '^V|^B|^S|^P') ~ 'vegetacion',# bosque, selva, pastizal, sabana, palmar, etc.
    str_detect(DESCRIPCIO, '^AGRICULTURA') ~ 'agricultura', # does not include shifting cultivation (nómada)
    str_detect(CVE_UNION, '^ACUI|^H2O') ~ 'agua',
    str_detect(CVE_UNION, '^ADV|^DV') ~ 'sin_veg',
    str_detect(CVE_UNION, '^AH') ~ 'construido',
    TRUE ~ 'otro'
  ))

usv7 <- usv7 %>%
  group_by(tipo) %>%
  summarize() %>%
  st_simplify(dTolerance=40, preserveTopology=T) %>%
  nngeo::st_remove_holes(100000)

# Save dissolved USV 
usv7 %>% st_write(lc_7clas_fp, delete_dsn=T)

# Reclass and dissolve to Veg, Ag, and Other
key <- c(inducido='veg', pastizal='veg', vegetacion='veg', agricultura='ag',
         agua='otro', sin_veg='otro', construido='otro')
usv3 <- usv7 %>%
  mutate(tipo =  recode(tipo, !!!key)) %>%
  group_by(tipo) %>%
  summarize() %>%
  st_simplify(dTolerance=40, preserveTopology=T) %>%
  nngeo::st_remove_holes(100000)

# Save dissolved USV 
usv3 %>% st_write(lc_3clas_fp, delete_dsn=T)

# Convert to raster and stack ----
# Input pollinator richness 
fps <- list.files(file.path(pred_dir, 'richness'), '*\\.tif', full.names = TRUE)
rich_ras <- terra::rast(fps[[1]])

# Rasterize ANP zones
anps <- rasterize_zones(anp_zones_fp, rich_ras, 
                        field = 'zone', overwrite = FALSE)

# Rasterize ecoregions
biom <- rasterize_zones(biom_diss_fp, rich_ras, 
                        field = 'DESECON1', overwrite = TRUE)

# Rasterize landcover
usv <- rasterize_zones(usv_diss_fp, rich_ras, 
                        field = 'tipo', overwrite = TRUE)

# Combine zone rasters
combo_ras <- c(anps$r, biom$r, usv$r)
names(combo_ras) <- c('anp_zone', 'biome', 'landcover')

combo_ras %>% terra::writeRaster(zones_ras_fp, overwrite=T)

# Save lookup tables
saveRDS( list(anp = anps$lu, biom = biom$lu, landcover = usv$lu), file=zones_lu_fp)

# Rasterize ANP IDs
ras_fp <- str_c(tools::file_path_sans_ext(anp_terr_fp), '.tif')
anp_ids <- rasterize_zones(anps_terr, rich_ras, 
                           field = 'ID_ANP', 
                           ras_fp = ras_fp, 
                           overwrite = TRUE)

lu <- readRDS(zones_lu_fp)
lu$ID_ANP <- anp_ids$lu
saveRDS( lu, file=zones_lu_fp)
