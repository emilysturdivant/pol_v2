# Look at Coleoptera -----
library(tidyverse)

fp <- '/home/esturdivant/PROJECTS/polinizadores/data/input_data/GBIF/family_order_query/gbif_Coleoptera.rds'
df <- readRDS(fp)
df_unq <- df %>% distinct(species, genus, family, decimalLatitude, decimalLongitude)
ct_species_Coleoptera <- df_unq %>% count(species, genus, family)

df_unq %>% filter(species == '') %>% nrow
df_unq %>% filter(species != '') %>% count(species) %>% arrange(desc(n))
df_unq %>% filter(genus == '') %>% nrow
df_unq %>% filter(genus != '') %>% count(genus) %>% arrange(desc(n))
df_unq %>% filter(family != '') %>% count(family) %>% arrange(desc(n))

# Look at convex hulls ----
pts_dir <- 'data/tidy/pollinator_points'
df1_fp <- file.path(pts_dir, 'points_nested_species.rds')
df2 <- readRDS(df1_fp)

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

# Species with few points or small range ----
summary_fp <- file.path('data/data_out', 'pollinator_points_summary.csv')
summary_tbl <- read_csv(summary_fp)

raras <- summary_tbl %>% 
  mutate(rare = case_when(nobs > 25 ~ 'n_comun', TRUE ~ 'n_rara')) %>% 
  count(common_group, subgroup_a, rare) %>% 
  pivot_wider(names_from = rare, values_from = n)

exts <- summary_tbl %>% 
  filter(nobs > 25) %>% 
  mutate(endemic = cut(area_km2, 
                       breaks = c(-Inf, 50000, 100000, Inf), 
                       labels = c('ext_lt50k', 'ext_lt100k', 'ext_gt100k'))) %>% 
  count(common_group, subgroup_a, endemic) %>% 
  pivot_wider(names_from = endemic, names_sort = TRUE, values_from = n)

dists_tbl <- raras %>% left_join(exts)

ggplot(summary_tbl, aes(x = area_km2)) +
  geom_histogram() +
  scale_x_log10()

summary_tbl %>% 
  filter(nobs < 26) %>% 
  ggplot(aes(x = nobs)) +
  geom_bar() +
  facet_wrap(vars(common_group))

# Mariposas. Most of the "rare" species are butterflies and of those most are nocturnal
rare_sps %>% 
  filter(common_group == 'Mariposas') %>% 
  ggplot(aes(x = nobs)) +
  geom_bar() +
  facet_wrap(vars(subgroup_a))

summary_tbl %>% 
  filter(subgroup_a == 'nocturna') %>% 
  ggplot(aes(x = nobs)) +
  geom_histogram() +
  facet_wrap(vars(superfamily))

# Area
summary_tbl %>% 
  filter(area_km2 < 10000) %>% 
  ggplot(aes(x = area_km2)) +
  geom_histogram() +
  facet_wrap(vars(common_group))
