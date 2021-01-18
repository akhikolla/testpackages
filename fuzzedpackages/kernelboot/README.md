
# kernelboot

[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/kernelboot)](https://CRAN.R-project.org/package=kernelboot)
[![GitHub Actions CI](https://github.com/twolodzko/kernelboot/workflows/CI/badge.svg)](https://github.com/twolodzko/kernelboot/actions?query=workflow%3ACI)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/twolodzko/kernelboot?branch=master&svg=true)](https://ci.appveyor.com/project/twolodzko/kernelboot)
[![Coverage Status](https://img.shields.io/codecov/c/github/twolodzko/kernelboot/master.svg)](https://codecov.io/github/twolodzko/kernelboot?branch=master)
[![Downloads](https://cranlogs.r-pkg.org/badges/grand-total/kernelboot)](https://CRAN.R-project.org/package=kernelboot)


This package implements random generation procedures for sampling from kernel
densities and smoothed bootstrap, that is an extension of standard bootstrap
procedure, where instead of drawing samples with replacement from the empirical
distribution, they are drawn from kernel density estimate of the distribution.

Three functions are provided to sample from univariate kernel densities (`ruvk`),
multivariate product kernel densities (`rmvk`) and multivariate Gaussian kernel
densities (`rmvg`). The `ruvk` function samples from the kernel densities as 
estimated using the base R `density` function. It offers possibility of sampling
from kernel densities with Gaussian, Epanechnikov, rectangular, triangular, biweight,
cosine, and optcosine kernels. The `rmvk` offers sampling from a multivariate kernel
density constructed from independent univariate kernel densities. It is also possible
to sample from multivariate Gaussian kernel density using the `rmvg` function,
that allows for correlation between the variables.

Smooth bootstrap is possible by using the `kernelboot` function, that draws with
replacement samples from the empirical distribution, enhances them using noise
drawn from the kernel density and evaluates the user-provided statistic on the
samples. This procedure can be thought as an extension of the basic bootstrap
procedure.
