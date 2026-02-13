# IP-35-Joint-Species-Distribution-Modeling-HMSC

This product is built on the Hierarchical Modelling of Species Communities (HMSC) framework, a Bayesian hierarchical joint species distribution modeling (JSDM) approach applied to benthic seaweed and invertebrate communities in the Baltic Sea region (Estonian and Finnish coastal area). HMSC enables the simultaneous modeling of multiple species while accounting for co-occurrence patterns and environmental dependencies. By accommodating the multivariate nature of community data, it is particularly well-suited for benthic seaweed ecosystems, where shared environmental responses and strong interactions among habitat-forming species, seaweeds, and invertebrates are expected. In addition, the Bayesian framework in HMSC provides a robust approach for quantifying uncertainties in parameter estimates and predictions, which is essential given the inherent measurement errors and natural variability in ecological data. HMSC is especially useful for modeling host species linked to habitat-formers, which single-species models struggle with. It also provides a strong community perspective, key for habitat modeling.

# Prerequisites

Software
-
R

Required R packages
-
Loaded in hmsc_script.R:

data.table

Hmsc

terra

corrplot

abind

ggplot2

pROC

Hmisc

Loaded in hmsc_abikoodid.r:

abind (and it assumes data.table + ggplot2 are already loaded by the main script when sourced)

# Input

The script assumes you are using rds files. Check "training_dataset_example.csv", "abiotics_cop_keskmised_example.csv" and "bs1km_vs_cop_example.csv" for examples for the table structures

Training dataset
-


Loaded into sisendandmed and converted to data.table

sisendandmed=readRDS("QuantitativeSamplesBiomassesKeySpeciesWithCopernicusDataAndDepthAndMoreHistory.rds")

Contains abiotic environmental variables and species observations, currently presence and absence data, but can also include cover or biomass.


Important
-
Species dependent variables and environmental covariates independent variables are selected by hard coded column positions (liigid = species, keskkond2 = environmental covariates, logtunnused = variables that need to be log transformed first):

liigid=names(sisendandmed)[c(23,28,34,35,...)]

keskkond2=names(sisendandmed)[c(11,128,132,...)]

logtunnused=names(sisendandmed)[c(11)]

Prediction datasets
-
abiotics_cop_keskmised.rds: abiotic covariates for the prediction grid (cop = native Copernicus product grid, bs1km = our own 1 km resolution Baltic Sea grid); later merged with a mapping table.

andmedvalja=readRDS("abiotics_cop_keskmised.rds")



bs1km_vs_cop.rds: mapping table providing grid coordinates and IDs

abibs=readRDS("bs1km_vs_cop.rds")



# Output

Tables (CSV)
-
- liigid_rruut.csv
  Per-species “R²-like” fit metric on training predictions:
  - For "esinemine": Tjur R² (abi$TjurR2)
  - Otherwise: standard R² (abi$R2)

- R2_testandmetel.csv (only if testime=TRUE)
  Per-species test-set metrics:
  - For occurrence: includes R2, accuracy_metric, discrimination_metric (AUC), calibration_metric, precision_metric
  - For biomass/logbiomass: squared correlation between observed and predicted

- R2_testandmetel_yksik.csv (only if testime=TRUE and yheliigimudelid=TRUE)
Same as above but for single-species models.

Plots (PDF)
-
- variance_partioning.pdf
Barplots per species: variance partitioning across covariate groups and the spatial random effect.
- liik_ja_keskkond.pdf
Species–environment coefficient support plot (plotBeta(..., param="Support", supportLevel=0.95)).

- liik_ja_liik.pdf
Species association heatmap (correlation of residuals / random-effect level associations), with significance filtering by posterior support.

- For each covariate in tunnusenimed = unique(c(keskkond, keskkond2)), several PDFs:
  - *_kesk.pdf: gradient plots with other covariates set to “mean” settings

  - *_soltuv.pdf: gradient plots with other covariates set to “most probable” settings conditional on focal covariate
  
  - *_soltuv_kategooria.pdf: extra set when a categorical covariate exists; stratifies plots by category level

Spatial rasters (GeoTIFF):
- 
For each species name in the prediction output matrix:

- <species_name>.tif
Raster in EPSG:3035 created from predicted values over the 1 km grid. Values are clamped with pmax(pred,0) before rasterisation.

# Methodology

This R script builds a spatial joint species distribution model (JSDM) for benthic species using Hmsc: it loads an RDS dataset, selects target species columns and environmental covariates, filters by year, bins samples into 20 depth classes and does stratified random train/test splitting within each depth bin (with extra downsampling of the shallowest bins to manage imbalance), jitters coordinates slightly to avoid duplicated locations, and optionally log-transforms chosen predictors; it then constructs the response matrix Y as either presence/absence (thresholding biomass > 0), raw biomass, or log(biomass+1), creates predictor data XData and a model formula that can include quadratic (poly(...,2)) terms, and fits an Hmsc model with a spatial random effect implemented via a Gaussian Predictive Process (GPP) using a set of knots (dimension-reduction for spatial autocorrelation), estimating parameters by MCMC sampling (multiple chains/parallel); after fitting it computes variance partitioning, posterior “significant” species–environment effects (Beta support), residual species–species association matrices (Omega correlations), and in-sample fit metrics (R² or Tjur R²), then evaluates out-of-sample performance on the held-out test set by predicting expected values and summarizing metrics (e.g., Tjur-style separation, absolute error, AUC via ROC, calibration-by-bins, and a precision term), optionally repeats the whole workflow as single-species models for comparison, produces gradient-response plots for marginal effects, and finally generates spatial predictions over a gridded abiotic dataset in chunks (to control memory/time) and writes each species’ prediction surface to GeoTIFF rasters.

# Usage instructions

Folder layout
-
Place in one working directory:

hmsc_script.R

hmsc_abikoodid.R 

the required .rds files

Step-by-step run
-

1. Install required packages (once):
   install.packages(c("data.table","Hmsc","terra","corrplot","abind","ggplot2","pROC","Hmisc"))

2. Open hmsc_script.R and update:
- setwd("...") to your directory
- sisendandmed=readRDS("...") to the correct training RDS filename
- the “KASUTAJA SISENDI” block: species (liigid), covariates (keskkond, keskkond2), years, bounding box, response type, etc.

3. Run the script from a clean R session.

4. Review outputs in the working directory (PDFs, CSVs, GeoTIFFs, and tooseis.RData).





