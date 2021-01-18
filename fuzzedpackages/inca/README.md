# inca: an R package for integer calibration

### Authors
[Luca Sartore](mailto://luca.sartore@usda.gov) and [Kelly Toppin](mailto://kelly.toppin@nass.usda.gov)

Maintainer: [Luca Sartore](mailto://drwolf85@gmail.com)

[![](http://www.r-pkg.org/badges/version/inca)](http://www.r-pkg.org/pkg/inca)
[![CRAN RStudio mirror downloads](http://cranlogs.r-pkg.org/badges/inca)](http://www.r-pkg.org/pkg/inca)
[![Mentioned in Awesome Official Statistics ](https://awesome.re/mentioned-badge.svg)](http://www.awesomeofficialstatistics.org)

*This research was supported in part by the U.S. Department of Agriculture, National Agriculture Statistics Service. The findings and conclusions in this publication are those of the authors and should not be construed to represent any official USDA, or US Government determination or policy.*

## Features of the package
Calibration forces the weighted estimates of calibration variables to match known totals. This improves the quality of the design-weighted estimates. It is used to adjust for non-response and/or under-coverage. The commonly used methods of calibration produce non-integer weights. In cases where weighted estimates must be integers, one must "integerize" the calibrated weights. However, this procedure often produces final weights that are very different for the "sample" weights. To counter this problem, the **inca** package provides specific functions for rounding real weights to integers, and performing an integer programming algorithm for calibration problems with integer weights.

For a complete list of exported functions, use `library(help = "inca")` once the **inca** package is installed (see the `inst/INSTALL.md` file for a detailed description of the setup process).

### Example
```R
library(inca)
set.seed(0)
w <- rpois(150, 4)
data <- matrix(rbinom(150000, 1, .3) * rpois(150000, 4), 1000, 150)
y <- data %*% w
w <- runif(150, 0, 7.5)
print(sum(abs(y - data %*% w)))
cw <- intcalibrate(w, ~. + 0, y, lower = 1, upper = 7, sparse = TRUE, data = data)
print(sum(abs(y - data %*% cw)))
barplot(table(cw), main = "Calibrated integer weights")
```

## References
Theberge, A. (1999). Extensions of calibration estimators in survey sampling.  *Journal of the American Statistical Association*, **94**(446), 635-644.

Little, R. J., & Vartivarian, S. (2003). On weighting the rates in non-response weights.

Kish, L. (1992). Weighting for unequal Pi.  *Journal of Official Statistics*, **8**(2), 183.

Rao, J. N. K., & Singh, A. C. (1997). A ridge-shrinkage method for range-restricted weight calibration in survey sampling.  *In Proceedings of the section on survey research methods* (pp. 57-65). American Statistical Association, Washington, DC.

Horvitz, D. G., & Thompson, D. J. (1952). A generalization of sampling without replacement from a finite universe.  *Journal of the American Statistical Association*, **47**(260), 663-685.

Kalton, G., & Flores-Cervantes, I. (2003). Weighting methods.  *Journal of Official Statistics*, **19**(2), 81-98.

Sartore, L., Toppin, K., Young, L., Spiegelman, C. 2019. Developing integer calibration weights for the Census of Agriculture.  *Journal of Agricultural, Biological, and Environmental Statistics*, **24**(1), 26-48.
