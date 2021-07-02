# Pre-process GBIF data compiled by Mauricio Quesada's team
# parameters: https://www.gbif.org/developer/occurrence#parameters

# Load libraries ---------------------------------------------------------------
# library(tabulizer)
library(sf)
# library(rgdal)
library(tools)
# library(raster)
# library(scales)
library(taxize)
library(tidyverse)

# Initialize ----
data_dir <- 'data/input_data/Quesada_bioclim_pol_y_cultivos'
pts_dir <- 'data/tidy/pollinator_points'
dir.create(pts_dir, recursive = TRUE, showWarnings = FALSE)

# Load and tidy points for each pollinator group ----
pts_dir_in <- file.path(data_dir, 'Completos')
(fps <- list.files(pts_dir_in, '.csv', full.names = T))

for(fp in fps) { 
  # Extract pollinator group name from filename
  (name <- fp %>% 
     basename %>% file_path_sans_ext %>% 
     str_split('_') %>% last %>% last %>% 
     str_replace('Muercielagos', 'Murcielagos'))
  
  # file path for output tidied points
  pts_fp <- file.path(pts_dir, str_c(name, '.gpkg'))
  
  if(!file.exists(pts_fp)) {
    
    # Set temp folder
    temp_fp <- file.path(pts_dir_in, 'extras', str_c(name, '.csv'))
    dir.create(dirname(temp_fp), recursive = TRUE, showWarnings = FALSE)
    
    # Load data with additional taxonomic information
    if(!file.exists(temp_fp)){
      
      dat <- read_csv(fp)
      
      # List unique genuses
      gen_list <- dat %>% 
        dplyr::select(genus) %>% 
        distinct %>% 
        filter(!is.na(genus))
      
      # Use taxize to get tribe and superfamily
      new_taxons <- gen_list %>% 
        purrr::map_dfr(taxize::tax_name, get=c('tribe', 'superfamily'), db='both')
      
      # Tidy results
      new_taxons <- new_taxons %>% 
        group_by(query) %>% 
        fill(tribe, superfamily, .direction = 'down') %>% 
        fill(tribe, superfamily, .direction = 'up') %>% 
        slice(1) %>% 
        ungroup %>% 
        dplyr::select(genus = query, tribe, superfamily)
      
      # Join to df
      dat <- dat %>% 
        left_join(new_taxons) 
      
      if(name == 'Mariposas'){
        dat <- dat %>% 
          mutate(nocturna = case_when(
            superfamily %in% c('Papilionoidea', 'Hesperioidea', 'Hedyloidea') ~ 'diurna',
            !is.na(superfamily) ~ 'nocturna'
          ))
      }
      
      # Save
      dat %>% write_csv(temp_fp)
      
    } else {
      
      dat <- read_csv(temp_fp)
      
    }
    
    # Tidy points 
    # Extract certain variables
    vars <- c('species', 'genus', 'tribe', 'family', 'superfamily', 'order', 'class', 'nocturna', 
              'decimalLongitude', 'decimalLatitude', 
              'eventDate', 'coordinateUncertaintyInMeters', 'habitat', 
              'basisOfRecord', 'country', 'stateProvince', 'institutionCode')
    
    # Filter points: drop coords with NAs, drop imprecise
    df <- dat %>% 
      drop_na(decimalLongitude, decimalLatitude) %>% 
      filter(decimalLongitude != 0 & decimalLatitude != 0) %>% 
      filter(!str_detect(issues, 'COUNTRY_COORDINATE_MISMATCH')) %>% 
      filter(coordinateUncertaintyInMeters < 1000 | is.na(coordinateUncertaintyInMeters)) %>% 
      dplyr::select(matches(vars))
    
    # Drop duplicates and convert to points
    sf_df2 <- df %>% 
      distinct %>% 
      st_as_sf(x = .,                         
               coords = c("decimalLongitude", "decimalLatitude"),
               crs = 4326)
    
    # Save
    sf_df2 %>% st_write(pts_fp, append = FALSE)
  }
}

# Summarize points ----
# Function to load points and get summary values
load_pts <- function(fp) {
  # Get name of pollinator group
  pol_group <- basename(tools::file_path_sans_ext(fp))
  
  # Functions to apply to nested data
  convhull_fun <- function(df){
    st_union(df) %>% 
      st_convex_hull()
  }
  
  area_fun <- function(ch){
    st_area(ch) %>% 
      units::set_units('km^2') %>% 
      units::set_units(NULL) %>% 
      round(1)
  }
  
  # Read
  df <- st_read(fp)
  
  # Nest points by species
  pols_nest <- df %>% 
    select(-ends_with('Key'), -starts_with('verbatim')) %>% 
    group_by(species, genus, tribe, family, superfamily, order, class) %>% 
    nest()
  
  # Get summary values
  pols_hulls <- pols_nest %>% 
    mutate(common_group = pol_group,
           convhull = map(data, convhull_fun),
           area_km2 = map_dbl(convhull, area_fun), 
           nobs = map_int(data, nrow))
  
  # Return
  return(pols_hulls)
}

# file path for data file
df1_fp <- file.path(pts_dir, 'points_nested_species.rds')

if(!file.exists(df1_fp)) {
  
  # Load, nest, and bind pre-processed points
  fps <- list.files(pts_dir, 'gpkg$', full.names = TRUE)
  df1 <- fps %>% map_dfr(load_pts)
  
  # ~ Add supplemental information
  # Abejas: Sociabilidad and Anidación
  class_fp <- file.path(data_dir, 'datos_puros', 'Abejas_especies_unicasODC.xlsx')
  abeja_subgroups <- readxl::read_excel(class_fp, na = c("", "-"))
  abeja_subgroups <- abeja_subgroups %>% 
    select(species, subgroup_a = Sociabilidad, subgroup_b = Anidacion) %>% 
    mutate(species = str_remove_all(species, "\\(.*\\) "))
  
  # Join
  df2 <- df1 %>% 
    ungroup() %>% 
    left_join(abeja_subgroups, by = 'species')
  
  # Mariposas: diurna o nocturna?
  df2 <- df2 %>% 
    mutate(subgroup_a = case_when(
      common_group == 'Mariposas' & 
        superfamily %in% c('Papilionoidea', 'Hesperioidea', 'Hedyloidea') ~ 'diurna',
      common_group == 'Mariposas' ~ 'nocturna', 
      TRUE ~ subgroup_a
    )
    )
  
  # Save points
  df2 %>% saveRDS(df1_fp)
  
} else {
  
  df2 <- readRDS(df1_fp)
  
}

# Save summary table
summary_fp <- file.path('data/data_out', 'pollinator_points_summary.csv')
dir.create(dirname(summary_fp), recursive = TRUE, showWarnings = FALSE)
summary_tbl <- df2 %>% select(-data, -convhull) 
summary_tbl %>% write_csv(summary_fp)

# df2 %>% filter(str_detect(species, 'Leptonycteris'))
