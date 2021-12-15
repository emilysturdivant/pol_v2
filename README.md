# Pollinators of Mexico

Exploring the dynamic between pollinators and conservation land. This is a new and reduced version of the processing in [polinizadores](https://github.com/emilysturdivant/polinizadores).

## Process overview

1. For each species, create a model that predicts likelihood of species occupancy based on the environment variables (altitude, landuse, WorldClim).
2. Use each model to predict species occupancy across Mexico.
3. Sum the likelihood rasters into a richness raster.

## Inputs:

- species points from GBIF provided by Dr. Quesada 
- predictor variables: all 19 WorldClim bioclimatic variables, elevation, land use/landcover, biomes
- administrative boundaries for context
- polygons of natural protected areas (areas protegidas naturales, ANPs)

## Outputs

### Data 

- species distribution models (SDMs) for all species 
- dataframe with species modeling attributes
- richness layers by species groupings

### Figures

- richness maps
- richness distribution by ANP zone, ecoregion, and pollinator group

## Workflow

### Using the scripts

1. Modify parameters in `00_initialize.R`, such as parameters for the filtering of species observation points and parameters for the random forest model. 
3. If the point data needs to be re-processed, run `prep_Quesada_GBIF_data.R`. The function `add_taxon_info` requires manual input when taxize needs help in identifying the best match. 
4. To rerun the model, run `process_SDM_for_cluster.R`. 
5. To perform spot-check and sum the SDMs into richness layers, run `process_SDM.R`.
6. To perform analysis on the richness results, run `analyze_richness.R`.
