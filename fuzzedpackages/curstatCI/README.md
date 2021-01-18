
<!-- README.md is generated from README.Rmd. Please edit that file -->
curstatCI
=========

The goal of curstatCI is to obtain confidence intervals for the distribution function of a random variable based on current status data. In the current status model, the variable of interest *X* with distribution function *F*<sub>0</sub> is not observed directly. A censoring variable *T* is observed instead together with the indicator *Δ* = (*X* ≤ *T*). *curstatCI* provides functions to estimate the distribution function *F*<sub>0</sub> and to construct pointswise confidence intervals around *F*<sub>0</sub>(*t*) based on an observed sample (*T*<sub>1</sub>, *Δ*<sub>1</sub>),…,(*T*<sub>*n*</sub>, *Δ*<sub>*n*</sub>) of size *n* from the observable random vector (*T*, *Δ*).

Installation
------------

You can install curstatCI from CRAN with:

``` r
# install.packages("curstatCI")
```

The package *curstatCI* requires the library *Rcpp*. To use the functions available in *curstatCI* load:

``` r
load(Rcpp)
load(curstatCI)
```

You can install curstatCI from github with:

``` r
# install.packages("devtools")
devtools::install_github("kimhendrickx/curstatCI")
```

Example
-------

This is a basic example which shows you how to obtain the confidence intervals for the distribution function of the time to infection for the Rubella data set. More information on the data and usage of the package can be found in the vignette "curstatCI":

``` r
library(Rcpp)
library(curstatCI)
set.seed(1)
data(rubella) 
grid <-1:80
bw <-ComputeBW(data=rubella, x=grid)
#> The computations took    1.256   seconds
out<-ComputeConfIntervals(data=rubella,x=grid,alpha=0.05, bw = bw)
#> The program produces the Studentized nonparametric bootstrap confidence intervals for the cdf, using the SMLE.
#> 
#> Number of unique observations:    225
#> Sample size n =  230
#> Number of Studentized Intervals = 80
#> Number of Non-Studentized Intervals = 0
#> The computations took    0.414   seconds

out$MLE
#>          [,1]      [,2]
#>  [1,]  0.0000 0.0000000
#>  [2,]  0.9452 0.2000000
#>  [3,]  1.2301 0.4857143
#>  [4,]  5.6411 0.5000000
#>  [5,]  9.4603 0.5714286
#>  [6,] 12.4548 0.8571429
#>  [7,] 15.4466 0.8666667
#>  [8,] 20.8000 0.8750000
#>  [9,] 23.2219 0.9285714
#> [10,] 25.2932 0.9401709
#> [11,] 77.8027 1.0000000

smle <-  out$SMLE
left<-out$CI[,1]
right<-out$CI[,2]

ConfInt<-cbind(smle, left, right)
head(ConfInt)
#>            smle       left      right
#> [1,] 0.05205647 0.04645370 0.05684488
#> [2,] 0.10387679 0.09062863 0.11518718
#> [3,] 0.15522742 0.13450756 0.17291065
#> [4,] 0.20587997 0.17795836 0.22970949
#> [5,] 0.25561365 0.22085116 0.28528225
#> [6,] 0.30421750 0.26290425 0.33934792
```
