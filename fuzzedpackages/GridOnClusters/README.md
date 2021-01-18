The 'GridOnClusters' R package
===============================

[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/GridOnClusters)](https://cran.r-project.org/package=GridOnClusters)
[![CRAN_latest_release_date](https://www.r-pkg.org/badges/last-release/GridOnClusters)](https://cran.r-project.org/package=GridOnClusters)
[![metacran downloads](https://cranlogs.r-pkg.org/badges/GridOnClusters)](https://cran.r-project.org/package=GridOnClusters)
[![metacran downloads](https://cranlogs.r-pkg.org/badges/grand-total/GridOnClusters)](https://cran.r-project.org/package=GridOnClusters)
[![lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)



### Overview

The package offers a method to discretize multivariate continuous data using a grid that captures the joint distribution via preserving clusters in the original data (Wang et al. 2020). Joint grid discretization is applicable as a data transformation step before using other methods to infer association, function, or causality without assuming a parametric model.

### When to use the package

Most available discretization methods process one variable at a time, such as ['Ckmeans.1d.dp'](https://cran.r-project.org/package=Ckmeans.1d.dp). If discretizing each variable independently misses patterns arising from the joint distribution of multiple involved variables, one may benefit from using the joint discretization method in this package.

### To download and install the package

```{r}
install.packages("GridOnClusters")
```

### Examples

See the *Examples* vignette of the package.
