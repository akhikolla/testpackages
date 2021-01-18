# Hierarchical Inference Testing

[![CRAN](http://www.r-pkg.org/badges/version/hit)](https://cran.r-project.org/package=hit)

The current built and test status for Linux (Mac)
[![Build Status](https://travis-ci.org/QTCAT/hit.svg)](https://travis-ci.org/QTCAT/hit)
and for Windows 
[![Build status](https://ci.appveyor.com/api/projects/status/kttq4x98q6hra6ct?svg=true)](https://ci.appveyor.com/project/jrklasen/hit)
.

## Description:
Hierarchical inference testing (HIT) for (generalized) linear models with 
correlated covariates. HIT is furthermore applicable to high-dimensional 
settings. For details see:

**Mandozzi, J. and Buehlmann, P. (2015)**. *Hierarchical testing in the 
high-dimensional setting with correlated variables*. Journal of the American 
Statistical Association. [Preprint](https://arxiv.org/abs/1312.5556)

**Klasen, J. R. et al. (2016)**. *A multi-marker association method for genome-wide 
association studies without the need for population structure correction*. Nature 
Communications. [Paper](http://www.nature.com/articles/ncomms13299)

## Installation:
The package can be installed from CRAN,

```R
install.packages("hit")

```

or via 
[`devtools`](https://github.com/hadley/devtools#updating-to-the-latest-version-of-devtools)
, if you haven't `devtools` installed yet you have to do so first.

```R
# install.packages("devtools")
devtools::install_github("QTCAT/hit")
```

## Example:
The `hit`-function example gives an overview of the functionality of the 
package and can be accessed once the package is loaded.

```R
library(hit)
example(hit)
```

--------------------------------------------------------------------------------
[![License](https://img.shields.io/badge/license-GPL%20%28%3E=%202%29-brightgreen.svg)](https://www.gnu.org/licenses/gpl-2.0.html)
&copy; 2016 JR Klasen
