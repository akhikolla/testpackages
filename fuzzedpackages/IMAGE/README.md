
<!-- README.md is generated from README.Rmd. Please edit that file -->

# IMAGE

<!-- badges: start -->

<!-- badges: end -->

mQTL mapping in bisulfite sequencing studies by fitting a binomial mixed
model, incorporating allelic-specific methylation pattern.

## Installation

IMAGE requires the following R packages:

``` r
install.packages(c("foreach", "doParallel", "Matrix", "MASS"))
```

You can install the released version of IMAGE from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("IMAGE")
```

Or install the development version from:

``` r
library(devtools)
install_github("umich-biostatistics/IMAGE")
```

## Usage

The main function is image. You can find the instructions and an example
by ‘?image’.

Example:

``` r
data(ExampleData)
geno <- ExampleData$geno
K <- ExampleData$K
data <- ExampleData$data
res=image(geno,data,K)
```

## Results Reproduced

All the results from all methods used in the PMR-Egger paper can be
reproduced at <https://github.com/3211895/IMAGEreproduce>
