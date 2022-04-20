# TwinAnalysis

This is a package for structural equation modeling and twin analysis using OpenMx package for R. 
The package provides following functionality:
1. It computes standardized parameter estimates and their CIs in a RAM MxModel.
2. It provides functionality to gather the parameter estimates and their CIs in the neat tables and to return the tables in console or in an MS Excel file.
3. It fits the reference models (Saturated and Independence models) to compute fit statistics and compares nested models.
4. It computes descriptive statistics for twin data and cross-twin cross-trait correlations.
5. It provides the definition of a cross-lagged twin model.

The package uses my other package `mlth.data.frame` for the output functionality.

Author: Ivan Voronin, Université Laval, QC, Canada

e-mail: ivan.a.voronin@gmail.com

# To install:

```
# install.packages('remotes')

remotes::install_github('ivanvoronin/mlth.data.frame')
remotes::install_github('ivanvoronin/TwinAnalysis')
```
