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

# RF with downsampling ----
mutually_exclusive_pa <- FALSE

rf_vers <- 2

unq_code <- ifelse(unq_cells, 'unq_cells', 'unq_pts')
unq_code <- ifelse(mutually_exclusive_pa, 'excl', unq_code)
rf_name <- str_c('rf', str_c(rf_vers, unq_code, dfilt_code, sep='_'), cont_var_code)
fp_tail <- file.path('sdm', rf_name)
pred_dir <- file.path('data', 'data_out', fp_tail)
pred_dir
# File paths
sp_nospc <- str_replace(sp_name, ' ', '_')
model_fp <- file.path(pred_dir, 'models', str_c(sp_nospc, '.rds'))
eval_fp <- file.path(pred_dir, 'model_evals', str_c(sp_nospc, '.csv'))
erf_fp <- file.path(pred_dir, 'model_evals', str_c(sp_nospc, '.rds'))

if (file.exists(model_fp)) {
  print(str_c(sp_name, ': Model already created.'))
  rf <- readRDS(model_fp)
  erf <- readRDS(erf_fp)
  return(list(rf=rf, erf=erf))
}

# Presence points

# Convert to coords DF
sp1 <- tryCatch({
  mutate(sp_df,
         lon = unlist(map(sp_df$geometry, 1)), 
         lat = unlist(map(sp_df$geometry, 2)))
},
warning = function(w) {
  mutate(sp_df,
         lon = unlist(map(sp_df$geom, 1)), 
         lat = unlist(map(sp_df$geom, 2)))
}
) %>%
  st_drop_geometry() %>% 
  dplyr::select(lon, lat)

# Background points
pres_pts <- as_Spatial(sp_df)

set.seed(10)
if(mutually_exclusive_pa) {
  
  backg <- dismo::randomPoints(mask = pred, n=1000, p = pres_pts
                               # prob = T, # use mask as sampling bias grid
  )
  
  # bg <- spatSample(pred, 1000, "random", na.rm = TRUE, xy = TRUE)
  
} else {
  
  backg <- dismo::randomPoints(pred, n=1000) 
}

backg <- backg %>% 
  as_tibble() %>% 
  rename(lon = 'x', lat = 'y')

# Split presence into training and testing 
set.seed(0)
group <- kfold(sp1, 5)
train_1 <- sp1[group != 1, ]
test_1 <- sp1[group == 1, ]

# Split background into training and testing
set.seed(0)
group <- kfold(backg, 5)
train_0 <- backg[group != 1, ]
test_0 <- backg[group == 1, ]

# Extract environmental data 
# Training dataset 
train <- bind_rows(train_1, train_0)
envtrain1 <- raster::extract(pred, train, cellnumbers=T, df=T)

pb_train <- c(rep(1, nrow(train_1)), rep(0, nrow(train_0)))
envtrain1 <- data.frame( cbind(pa = pb_train, envtrain1) ) %>% 
  mutate(across(starts_with(c('biomes', 'ESA', 'usv')), as.factor)) %>% 
  dplyr::select(-ID)

# Remove duplicated cells
if(unq_cells) {
  envtrain1 <- envtrain1 %>% distinct()
}

envtrain <- envtrain1 %>% 
  dplyr::select(-cells) %>% 
  drop_na()

# Testing datasets - get predictors for test presence and background points
testpres <- data.frame( raster::extract(pred, test_1) ) %>%
  mutate(across(starts_with(c('biomes', 'ESA', 'usv')), as.factor))

testbackg <- data.frame( raster::extract(pred, test_0) ) %>%
  mutate(across(starts_with(c('biomes', 'ESA', 'usv')), as.factor))

# Training dataset 
test <- bind_rows(test_1, test_0)
pb_test <- c(rep(1, nrow(test_1)), rep(0, nrow(test_0)))
envtest1 <- raster::extract(pred, test, cellnumbers=T, df=T)
envtest1 <- data.frame( cbind(pa = pb_test, envtest1) ) %>% 
  mutate(across(starts_with(c('biomes', 'ESA', 'usv')), as.factor)) %>% 
  dplyr::select(-ID)
envtest <- envtest1 %>% 
  dplyr::select(-cells) %>% 
  drop_na()

# Set factor levels for test DFs to match training data
vars <- envtrain %>% dplyr::select(where(is.factor)) %>% tbl_vars
for(var in vars){
  levels(testpres[[var]]) <- levels(envtrain[[var]])
  levels(testbackg[[var]]) <- levels(envtrain[[var]])
}

# RF with downsampling ----
library('caret')

# Convert PA column to factor
envtrain <- envtrain %>% 
  mutate(pa = str_c('class_', pa), 
         pa = factor(pa))

envtest <- envtest %>% 
  mutate(pa = str_c('class_', pa), 
         pa = factor(pa))
  
# Get number of points in smallest class
nmin <- envtrain %>% count(pa) %>% select(n) %>% min()

# Train RF model
ctrl <- trainControl(method = "cv",
                     classProbs = TRUE,
                     summaryFunction = twoClassSummary)

set.seed(2)
rfDownsampled <- train(pa ~ ., 
                       data = envtrain,
                       method = "rf",
                       ntree = 1000,
                       tuneLength = 5,
                       metric = "ROC",
                       trControl = ctrl,
                       na.action = na.exclude,
                       # Tell randomForest to sample by strata. Here, 
                       ## that means within each class
                       strata = envtrain$pa,
                       ## Now specify that the number of samples selected
                       ## within each class should be the same
                       sampsize = rep(nmin, 2))

# Compute test set ROC curves 
downProbs <- predict(rfDownsampled, envtest, type = "prob")[,1]
downsampledROC <- pROC::roc(response = envtest$pa, 
                            predictor = downProbs,
                            levels = rev(levels(envtest$pa)))

# Plot curves
plot(downsampledROC, col = rgb(1, 0, 0, .5), lwd = 2)
pROC::auc(downsampledROC)
dismo::evaluate(testpres, testbackg, rfDownsampled)

# Compare to unbalanced RF ----
set.seed(2)
rfUnbalanced <- train(pa ~ ., 
                        data = envtrain,
                        method = "rf",
                        ntree = 1000,
                        tuneLength = 5,
                        metric = "ROC",
                        trControl = ctrl,
                      na.action = na.exclude)

# Compute test set ROC curves 
unbalProbs <- predict(rfUnbalanced, envtest, type = "prob")[,1]
unbalROC <- pROC::roc(response = envtest$pa, 
                predictor = unbalProbs,
                levels = rev(levels(envtest$pa)))

# Plot ROC curves
plot(unbalROC, col = rgb(0, 0, 1, .5), lwd = 2, add = TRUE)
pROC::auc(unbalROC)
dismo::evaluate(testpres, testbackg, rfUnbalanced)

# RF with repeated CV with downsampling ----
training <- bind_rows(mutate(backg, pa = 'class_0'),
                      mutate(sp1, pa = 'class_1'))
train1 <- raster::extract(pred, select(training, -pa), cellnumbers=T, df=T)
training <- data.frame( cbind(pa = training$pa, train1) ) %>% 
  mutate(across(starts_with(c('biomes', 'ESA', 'usv')), as.factor)) %>% 
  dplyr::select(-ID) %>% 
  dplyr::select(-cells) %>% 
  drop_na()

# Convert PA column to factor
training <- training %>% 
  mutate(pa = str_c('class_', pa), 
         pa = factor(pa))

# Get number of points in smallest class
nmin <- training %>% count(pa) %>% select(n) %>% min()

# Train RF model
ctrl <- trainControl(method = "repeatedcv",
                     number = 5,
                     repeats = 10,
                     classProbs = TRUE,
                     summaryFunction = twoClassSummary)

set.seed(2)
rfDownsampled <- train(pa ~ ., 
                       data = training,
                       method = "rf",
                       ntree = 1000,
                       tuneLength = 5,
                       metric = "Sens",
                       trControl = ctrl,
                       na.action = na.exclude,
                       # Tell randomForest to sample by strata. Here, 
                       ## that means within each class
                       strata = training$pa,
                       ## Now specify that the number of samples selected
                       ## within each class should be the same
                       sampsize = rep(nmin, 2))

# Make raster with downsampled model ----
rf <- rfDownsampled

# Filepaths  
sp_nospc <- str_replace(sp_name, ' ', '_')
likelihood_fp <- file.path(pred_dir, 'repeatCV', str_glue(sp_nospc, '.tif'))
dir.create(dirname(likelihood_fp), recursive = TRUE, showWarnings = FALSE)

# Create map and interpolate to fill holes
pr_rf1 <- dismo::predict(predictors, rf, ext=ext)
# pr_rf1 <- raster::focal(pr_rf1, 
#                         w=matrix(1,nrow=3, ncol=3), 
#                         fun=modal, 
#                         NAonly=TRUE, 
#                         na.rm=TRUE) 

# Save likelihood raster
writeRaster(pr_rf1, likelihood_fp, overwrite=T, 
            options=c("dstnodata=-99999"), wopt=list(gdal='COMPRESS=LZW'))
