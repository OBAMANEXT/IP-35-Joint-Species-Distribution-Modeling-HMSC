# IP-35-Joint-Species-Distribution-Modeling-HMSC

This product is built on the Hierarchical Modelling of Species Communities (HMSC) framework, a Bayesian hierarchical joint species distribution modeling (JSDM) approach applied to benthic seaweed and invertebrate communities. HMSC enables the simultaneous modeling of multiple species while accounting for co-occurrence patterns and environmental dependencies. By accommodating the multivariate nature of community data, it is particularly well-suited for benthic seaweed ecosystems, where shared environmental responses and strong interactions among habitat-forming species, seaweeds, and invertebrates are expected. In addition, the Bayesian framework in HMSC provides a robust approach for quantifying uncertainties in parameter estimates and predictions, which is essential given the inherent measurement errors and natural variability in ecological data. HMSC is especially useful for modeling host species linked to habitat-formers, which single-species models struggle with. It also provides a strong community perspective, key for habitat modeling.

# Prerequisites

Software:
R

# Input



# Output

Tables (CSV):
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

Plots (PDF):
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
s
