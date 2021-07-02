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

# Load and tidy points for each pollinator group ----
pts_dir_in <- file.path(data_dir, 'Completos')
(fps <- list.files(pts_dir_in, '.csv', full.names = T))
for(fp in fps) { 
  # Extract pollinator group name from filename
  (name <- fp %>% 
     basename %>% file_path_sans_ext %>% 
     str_split('_') %>% last %>% last %>% 
     str_replace('Muercielagos', 'Murcielagos'))
  
  # Load data 
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
  fp_out <- file.path(pts_dir, str_c(name, '.gpkg'))
  dir.create(dirname(fp_out), recursive = TRUE, showWarnings = FALSE)
  sf_df2 %>% st_write(fp_out, append = FALSE)
  
  # spX <- sf_df2 %>% filter(str_detect(species, 'Lophornis brachylophus'))
  # fp_out <- file.path('data/data_out/pollinator_points/no_duplicates',
  #                     str_c(name, '_Lophornis_brachylophus.gpkg'))
  # spX %>% st_write(fp_out, append = FALSE)
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
  
  # # QC: Check for non-matching rows
  # abejas <- summary_tbl %>% filter(common_group == 'Abejas')
  # ab <- abeja_subgroups %>% 
  #   mutate(species = str_remove_all(species, "\\(.*\\) "))
  # t0 <- nrow(ab) - nrow(distinct(ab, species))
  # t1 <- anti_join(abejas, ab, by = 'species')
  # t2 <- anti_join(ab, abejas, by = 'species')
  # if(t0 + nrow(t1) + nrow(t2) > 0) print('WARNING: not an exact match between species in the two datasets.')
  
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
  dir.create(dirname(df1_fp), recursive = TRUE, showWarnings = FALSE)
  df2 %>% saveRDS(df1_fp)
  
} else {
  
  df2 <- readRDS(df1_fp)
  
}

# Save summary table
summary_fp <- file.path('data/data_out', 'pollinator_points_summary.csv')
dir.create(dirname(summary_fp), recursive = TRUE, showWarnings = FALSE)
summary_tbl <- df2 %>% select(-data, -convhull) 
summary_tbl %>% write_csv(summary_fp)

df2 %>% filter(str_detect(species, 'Leptonycteris'))

# Look at convex hulls ----
convhull <- df2 %>% 
  filter(nobs > 25) %>% 
  filter(str_detect(species, 'Lophornis')) %>% 
  select(species, convhull)
chs2 <- map(convhull$convhull, st_as_sf)
single_sf <- dplyr::bind_rows(chs2) %>%
  bind_cols(convhull$species) %>% 
  rename(species = `...2`) %>% 
  left_join(select(df2, -data, -convhull))

single_sf %>% st_write(file.path(pts_dir, 'convex_hulls_species.gpkg'), append = FALSE)
single_sf %>% filter(area_km2 < 500000) %>% 
  st_write(file.path(pts_dir, 'convex_hulls_species2.gpkg'), append = FALSE)

library(tmap)
tmap_mode('view')
tm_shape(single_sf) + tm_polygons(alpha = 0.2)



# List species with few points or small range ----
summary_tbl %>% 
  mutate(rare = case_when(nobs > 25 ~ 'n_comun', TRUE ~ 'n_rara')) %>% 
  count(common_group, subgroup_a, rare) %>% 
  pivot_wider(names_from = rare, values_from = n)

summary_tbl %>% 
  filter(nobs > 25) %>% 
  mutate(endemic = cut(area_km2, 
                       breaks = c(-Inf, 500, 1000, 2000, 4000, Inf), 
                       labels = c('a_lt500', 'b_500to1000', 'c_1000to2000', 'd_2000to4000', 'e_gt4000'))) %>% 
  count(common_group, subgroup_a, endemic) %>% 
  pivot_wider(names_from = endemic, names_sort = TRUE, values_from = n)


ggplot(summary_tbl, aes(x = area_km2)) +
  geom_histogram() +
  scale_x_log10()

rare_sps <- summary_tbl %>% filter(nobs < 26)
rare_sps %>% arrange(nobs)

# Mariposas. Most of the "rare" species are butterflies and of those most are nocturnal
rare_sps %>% 
  filter(common_group == 'Mariposas') %>% 
  ggplot(aes(x = nobs)) +
  geom_bar() +
  facet_wrap(vars(subgroup_a))

summary_tbl %>% 
  filter(subgroup_a == 'nocturna') %>% 
  ggplot(aes(x = nobs)) +
  geom_bar()

summary_tbl %>% 
  filter(subgroup_a == 'nocturna') %>% 
  ggplot(aes(x = nobs)) +
  geom_histogram() +
  facet_wrap(vars(superfamily))






rare_sps %>% ggplot(aes(x = nobs)) +
  geom_bar() +
  facet_wrap(vars(common_group))

endmc_sps <- summary_tbl %>% filter(area_km2 < 1000)
endmc_sps %>% ggplot(aes(x = area_km2)) +
  geom_histogram() +
  facet_wrap(vars(common_group))


















# ~ Look at Coleoptera -----
fp <- '/home/esturdivant/PROJECTS/polinizadores/data/input_data/GBIF/family_order_query/gbif_Coleoptera.rds'
df <- readRDS(fp)
df_unq <- df %>% distinct(species, genus, family, decimalLatitude, decimalLongitude)
ct_species_Coleoptera <- df_unq %>% count(species, genus, family)

df_unq %>% filter(species == '') %>% nrow
df_unq %>% filter(species != '') %>% count(species) %>% arrange(desc(n))
df_unq %>% filter(genus == '') %>% nrow
df_unq %>% filter(genus != '') %>% count(genus) %>% arrange(desc(n))
df_unq %>% filter(family != '') %>% count(family) %>% arrange(desc(n))
