# gaselect R package
This R package implements a genetic algorithm (GA) for variable selection as described in
Kepplinger, D., Filzmoser, P., and Varmuza, K. (2017). Variable selection with genetic algorithms using repeated cross-validation of PLS regression models as fitness measure. https://arxiv.org/abs/1711.06695.

## Installation
To install the latest release from CRAN, run the following R code in the R console:
```r
install.packages('gaselect')
```

The most recent stable version as well as the developing version might not yet be available on
CRAN. These can be directly installed from github using the
[devtools](https://cran.r-project.org/package=devtools) package:
```r
# Install the most recent stable version:
install_github('dakep/gaselect')
# Install the (unstable) develop version:
install_github('dakep/gaselect', ref = 'develop')
```
