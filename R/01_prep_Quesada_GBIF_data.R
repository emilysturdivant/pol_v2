# Pre-process GBIF data compiled by Dr. Quesada's team

# Initialize ----
source('R/initialize.R')
dir.create(pts_dir, recursive = TRUE, showWarnings = FALSE)
data_dir <- 'data/input_data/Quesada_bioclim_pol_y_cultivos'
pts_dir_in <- file.path(data_dir, 'Completos')

# Add tribe and superfamily to data ----
add_taxon_info <- function(fp) {
  # Extract pollinator group name from filename
  name <- str_extract(fp, '(?<=completo_)[:alpha:]+(?=\\.csv)') %>%
    str_replace('Muercielagos', 'Murcielagos')
  
  # file path for output tidied points
  temp_fp <- file.path(temp_dir, str_c(name, '.csv'))
  
  # Load data with additional taxonomic information
  dat <- read_csv(fp) %>% 
    mutate(genus = na_if(genus, ''), 
           species = na_if(species, ''), 
           genus = if_else(is.na(genus), 
                           str_extract(species, '[:alpha:]+'),
                           genus) )
  
  # Use taxize to get tribe and superfamily for each unique genus
  new_taxons <- dat %>%
    distinct(genus) %>%
    filter(!is.na(genus)) %>%
    purrr::map_dfr(taxize::tax_name, get=c('tribe', 'superfamily'), db='both')
  
  # Tidy results
  new_taxons <- new_taxons %>%
    group_by(query) %>%
    fill(tribe, superfamily, .direction = 'down') %>%
    fill(tribe, superfamily, .direction = 'up') %>%
    slice(1) %>%
    ungroup() %>%
    dplyr::select(genus = query, tribe, superfamily)
  
  # Join to df
  dat <- dat %>%
    left_join(new_taxons) %>%
    mutate(grupo = name)
  
  if(name == 'Mariposas'){
    dat <- dat %>%
      mutate(nocturna = case_when(
        superfamily %in% c('Papilionoidea', 'Hesperioidea', 'Hedyloidea') ~ 'diurna',
        !is.na(superfamily) ~ 'nocturna'
      ))
  }
  
  # Save
  dat %>% write_csv(temp_fp)
}

# Create dir
temp_dir <- file.path(pts_dir_in, 'extras')
dir.create(temp_dir, recursive = TRUE, showWarnings = FALSE)

# Run for every pollinator group
(fps <- list.files(pts_dir_in, '.csv', full.names = T))
fps %>% purrr::walk(add_taxon_info)

# Merge points ----
vars <- c('species', 'genus', 'tribe', 'family', 'superfamily', 'order', 'class', 
          'decimalLongitude', 'decimalLatitude', 'eventDate', 
          'coordinateUncertaintyInMeters', 'issues', 'habitat', 
          'basisOfRecord', 'country', 'stateProvince', 'institutionCode', 'nocturna')

# Load dataframe of all observation points, unfiltered
(fps <- list.files(file.path(pts_dir_in, 'extras'), '.csv', full.names = T))
df_all <- fps %>% 
  purrr::map_dfr( ~ {
    read_csv(.x) %>% 
      dplyr::select(any_of(vars)) %>% 
      mutate(grupo = str_extract(.x, '(?<=/)[:alpha:]+(?=\\.csv)'))
  }
  )

# Replace 0s with NA
df_all <- df_all %>% 
  mutate(across(starts_with('decimal'), ~ na_if(.x, 0))) %>% 
  mutate(genus = na_if(genus, ''), 
         species = na_if(species, ''), 
         genus = if_else(is.na(genus), 
                         str_extract(species, '[:alpha:]+'),
                         genus) )

# Save
all_fp <- file.path(temp_dir, 'all_unfiltered.rds')
df_all %>% saveRDS(all_fp)
# df_all <- readRDS(all_fp)

# Filter points ----
# drop coords with NAs, drop imprecise 
df <- df_all %>% 
  drop_na(decimalLongitude, decimalLatitude) %>% 
  filter(str_detect(issues, 'COUNTRY_COORDINATE_MISMATCH', negate = TRUE)) %>% 
  filter(coordinateUncertaintyInMeters < max_coord_uncertainty | 
           is.na(coordinateUncertaintyInMeters))

# Drop duplicates and convert to points
sf_df <- df %>% 
  distinct %>% 
  st_as_sf(x = .,                         
           coords = c("decimalLongitude", "decimalLatitude"),
           crs = 4326)

# Optionally filter to date range
if(filt_dates) {
  date_min <- as.POSIXct(str_c(date_range[[1]], "-01-01"))
  date_max <- as.POSIXct(str_c(date_range[[2]], "-12-31"))
  sf_df <- sf_df %>% filter(eventDate >= date_min & eventDate <= date_max)
}

# Save
pts_fp <- file.path(pts_dir, 'points_all.gpkg')
sf_df %>% st_write(pts_fp, append = FALSE)

# Nest points by species ----
# Functions to apply to nested data
convhull_fun <- function(df){
  suppressWarnings( st_convex_hull(st_union(df)) )
}

area_fun <- function(ch){
  st_area(ch) %>% 
    units::set_units('km^2') %>% 
    units::set_units(NULL) %>% 
    round(1)
}

# Read
sf_df <- st_read(pts_fp)

# Nest points by species
pol_df <- sf_df %>% 
  select(-country, -issues) %>% 
  group_by(grupo, species, genus, tribe, family, superfamily, order, class) %>% 
  nest() %>% 
  mutate(convhull = map(data, convhull_fun),
         area_km2 = map_dbl(convhull, area_fun), 
         nobs = map_int(data, nrow)) %>% 
  ungroup()

# Add species that do not have sufficient data
df_dropped <- df_all %>% 
  select(-country, -issues) %>% 
  distinct(grupo, species, genus, tribe, family, superfamily, order, class) %>% 
  anti_join(pol_df)

pol_df <- bind_rows(pol_df, df_dropped)

# Divide mariposas into butterflies and moths ----
pol_df <- pol_df %>% 
  mutate(
    nocturna = case_when(
      grupo == 'Mariposas' & 
        superfamily %in% c('Papilionoidea', 'Hesperioidea', 'Hedyloidea') ~ 'diurna', # Rhopalocera
      grupo == 'Mariposas' ~ 'nocturna', # Heterocera
      TRUE ~ as.character(NA)
    )
  )

# Create grouping variable
pol_df <- pol_df %>% 
  mutate( 
    group = case_when( 
      grupo == 'Mariposas' & nocturna == 'diurna' ~ 'butterfly', 
      grupo == 'Mariposas' & nocturna == 'nocturna' ~ 'moth',
      grupo == 'Abejas' ~ 'bee',
      grupo == 'Avispas' ~ 'wasp',
      grupo == 'Colibries' ~ 'hummingbird',
      grupo == 'Moscas' ~ 'fly',
      grupo == 'Murcielagos' ~ 'bat',
      TRUE ~ as.character(NA)
    )
  )

# Add bee sociality and nesting ----
# Load Abejas: Sociabilidad and Anidación
class_fp <- file.path(data_dir, 'datos_puros', 'Abejas_especies_unicasODC.xlsx')
abeja_groups_raw <- readxl::read_excel(class_fp, na = c("", "-"))
abeja_subgroups <- abeja_groups_raw %>% 
  select(species, genus, Sociabilidad, Anidacion) %>% 
  mutate(species = str_remove_all(species, "\\(.*\\) "),
         bee_sociality = case_when(
           Sociabilidad == 'Solitaria' ~ 'solitary',
           Sociabilidad %in% c('Eusocial', 'Social') ~ 'social',
           Sociabilidad %in% c('Solitarias/Sociales', 'Solitaria/semisocial') ~ 'solitary/social',
           TRUE ~ as.character(NA)
         ),
         bee_nesting = case_when(
           Anidacion == 'Oquedades' ~ 'cavities',
           Anidacion == 'Tierra' ~ 'ground',
           Anidacion %in% c('tallos', 'Tallos') ~ 'stems',
           Anidacion %in% c('Tierra/tallos', 'Tierra/tallos/aereos', 'Aereos') ~ 'above-ground',
           TRUE ~ as.character(NA)
         )
  )

soc_grps <- abeja_subgroups %>%
  distinct(genus, bee_sociality) %>%
  group_by(genus) %>%
  filter(n() < 2) %>% 
  ungroup()

nid_grps <- abeja_subgroups %>%
  distinct(genus, bee_nesting) %>%
  group_by(genus) %>%
  filter(n() < 2) %>% 
  ungroup()

# Add bee sociality and nesting
pol_df <- pol_df %>% 
  left_join(
    dplyr::select(abeja_subgroups, species, bee_sociality, bee_nesting), 
    by = 'species')

soc_patch <- pol_df %>% 
  distinct(species, genus, bee_sociality) %>% 
  filter(is.na(bee_sociality)) %>% 
  select(-bee_sociality) %>% 
  left_join(soc_grps, by = 'genus')

nest_patch <- pol_df %>% 
  distinct(species, genus, bee_nesting) %>% 
  filter(is.na(bee_nesting)) %>% 
  select(-bee_nesting) %>% 
  left_join(nid_grps, by = 'genus')

pol_df2 <- pol_df %>% 
  rows_patch(soc_patch, by = 'species') %>% 
  rows_patch(nest_patch, by = 'species')

# Save points ----
pol_df2 %>% saveRDS(pts_nested_rds)
pol_df2 <- readRDS(pts_nested_rds)

# Filter points for modeling ----
# Load and filter GBIF points for given name
if( !file.exists(filt_pts_rds) ) {
  dat <- pol_df2 %>% 
    filter(!is.na(nobs), nobs > 24,
           !is.na(genus), !is.na(species))
  dat %>% saveRDS(filt_pts_rds)
} 

# Save summary table ----
summary_fp <- file.path('data/data_out', 'pollinator_species_summary.csv')
dir.create(dirname(summary_fp), recursive = TRUE, showWarnings = FALSE)
summary_tbl <- pol_df %>% select(-data, -convhull) 
summary_tbl %>% write_csv(summary_fp)