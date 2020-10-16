# TwinAnalysis

This is a package to simplify structural equation modeling - and particularly, twin analysis - in R. 
The package is based on OpenMx and provides the following functionality:
1. To define standardized parameters for any RAM MxModel.
2. To define output tables (as MxAlgebra) within any MxModel. To extract output tables, to merge estimates and confidence
intervals, to write the output tables into MS Excel file.
3. To fit reference (Saturated and Independence) models, to compute fit statistics, to compare nested models.
4. To compute twin statistics, to fit twin models (at the moment only cross-lag twin model is implemented).
5. To introduce your favorite twin model.

The package relies mostly on the basic functionality of OpenMx providing the shortcut between model definition 
and publishing the results, especially when similar models used to describe different data. 
It also uses `mlth.data.frame` which is another development of mine, the package for
fancy multi-layered tables in R.

Author: Ivan Voronin, Psychological Institute of Russian Academy of Education, Moscow

e-mail: ivan.a.voronin@gmail.com

# To install:

```
# install.packages('devtools')

library(devtools)

install_github('ivanvoronin/mlth.data.frame')
install_github('ivanvoronin/TwinAnalysis')
```
