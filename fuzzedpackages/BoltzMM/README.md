
<img src="http://www.r-pkg.org/badges/version-last-release/BoltzMM"></img></a>
[![Downloads from the RStudio CRAN
mirror](http://cranlogs.r-pkg.org/badges/BoltzMM)](https://CRAN.R-project.org/package=BoltzMM)
[![Build
Status](https://travis-ci.org/andrewthomasjones/BoltzMM.svg?branch=master)](https://travis-ci.org/andrewthomasjones/BoltzMM)
[![status](http://joss.theoj.org/papers/23eb189a5e0bdd2b51f668621abcc75a/status.svg)](http://joss.theoj.org/papers/23eb189a5e0bdd2b51f668621abcc75a)

<!-- README.md is generated from README.Rmd. Please edit that file -->

# BoltzMM

The BoltzMM package allows for computation of probability mass functions
of fully-visible Boltzmann machines (FVBMs) via `pfvbm` and `allpfvbm`.
Random data can be generated using `rfvbm`. Maximum pseudolikelihood
estimation of parameters via the MM algorithm can be conducted using
`fitfvbm`. Computation of partial derivatives and Hessians can be
performed via `fvbmpartiald` and `fvbmHessian`. Covariance estimation
and normal standard errors can be computed using `fvbmcov` and
`fvbmstderr`.

## Installation

If `devtools` has already been installed, then the most current build of
`BoltzMM` can be obtained via the
command:

``` r
devtools::install_github('andrewthomasjones/BoltzMM',build_vignettes = TRUE)
```

The latest stable build of `BoltzMM` can be obtain from CRAN via the
command:

``` r
install.packages("BoltzMM", repos='http://cran.us.r-project.org')
```

An archival build of `BoltzMM` is available at
<http://doi.org/10.5281/zenodo.2538256>. Manual installation
instructions can be found within the *R* installation and administration
manual <https://cran.r-project.org/doc/manuals/r-release/R-admin.html>.

## Examples

Compute the probability of every length n=3 binary spin vector under
bvec and Mmat:

``` r
library(BoltzMM)
set.seed(1)

bvec <- c(0,0.5,0.25)
Mmat <- matrix(0.1,3,3) - diag(0.1,3,3)
allpfvbm(bvec,Mmat)
#>           [,1]       [,2]      [,3]      [,4]       [,5]       [,6]
#> [1,] 0.0666189 0.04465599 0.1213876 0.1213876 0.07362527 0.07362527
#>           [,7]      [,8]
#> [1,] 0.2001342 0.2985652
```

Generate num=1000 random strings of n=3 binary spin variables under bvec
and Mmat.

``` r
library(BoltzMM)
set.seed(1)

num <- 1000
bvec <- c(0,0.5,0.25)
Mmat <- matrix(0.1,3,3) - diag(0.1,3,3)
data <- rfvbm(num,bvec,Mmat)

head(data)
#>      [,1] [,2] [,3]
#> [1,]    1    1   -1
#> [2,]   -1   -1    1
#> [3,]   -1    1    1
#> [4,]    1    1    1
#> [5,]   -1    1   -1
#> [6,]    1    1    1
```

Fit a fully visible Boltzmann machine to data, starting from parameters
bvec and Mmat.

``` r
library(BoltzMM)
set.seed(1)

bvec <- c(0,0.5,0.25)
Mmat <- matrix(0.1,3,3) - diag(0.1,3,3)
data <- rfvbm(num,bvec,Mmat)

fitfvbm(data,bvec,Mmat)
#> $pll
#> [1] -1892.661
#> 
#> $bvec
#> [1] 0.02607382 0.46484595 0.27640931
#> 
#> $Mmat
#>           [,1]      [,2]      [,3]
#> [1,] 0.0000000 0.1179001 0.1444486
#> [2,] 0.1179001 0.0000000 0.0351134
#> [3,] 0.1444486 0.0351134 0.0000000
#> 
#> $itt
#> [1] 5
```

Example with real data from
<https://hal.archives-ouvertes.fr/hal-01927188v1>.

``` r
# Load bnstruct library & package
library(bnstruct)
#> Loading required package: bitops
#> Loading required package: Matrix
#> Loading required package: igraph
#> 
#> Attaching package: 'igraph'
#> The following objects are masked from 'package:stats':
#> 
#>     decompose, spectrum
#> The following object is masked from 'package:base':
#> 
#>     union
library(BoltzMM)

# Load data
data(senate)

# Turn data into a matrix
senate_data <- as.matrix(senate)

# Recode Yes as 1, and No as -1
senate_data[senate=="Yes"] <- 1
senate_data[senate=="No"] <- -1

# Conduct imputation
imp_data <- knn.impute(suppressWarnings(matrix(as.numeric(senate_data),
                                        dim(senate_data))),
                       k=1)

# No governement - using as reference level
data_nogov <- imp_data[,-1]


# Initialize parameters
bvec <- rep(0,8)
Mmat <- matrix(0,8,8)
nullmodel<-list(bvec=bvec,Mmat=Mmat)

# Fit a fully visible Boltzmann machine to data, starting from parameters bvec and Mmat.
model <- fitfvbm(data_nogov,bvec,Mmat)
# Compute the sandwich covariance matrix using the data and the model.
covarmat <- fvbmcov(data_nogov,model,fvbmHess)
# Compute the standard errors of the parameter elements according to a normal approximation.
st_errors <- fvbmstderr(data_nogov,covarmat)
# Compute z-scores and p-values under null
test_results<-fvbmtests(data_nogov,model,nullmodel)

test_results
#> $bvec_z
#> [1] -1.3871285  3.0958110  1.8099814 -0.4957960 -0.8230061 -0.1625973
#> [7]  0.5715010  2.4308532
#> 
#> $bvec_p
#> [1] 0.165402598 0.001962754 0.070298676 0.620038322 0.410504513 0.870835547
#> [7] 0.567660056 0.015063316
#> 
#> $Mmat_z
#>            [,1]       [,2]       [,3]       [,4]       [,5]       [,6]
#> [1,]         NA -0.6596632 -1.3152156 -2.0181850  0.1936413  3.8587216
#> [2,] -0.6596632         NA -0.9536614 -2.0795939 -0.4947831  6.1625534
#> [3,] -1.3152156 -0.9536614         NA  1.3258059 -1.0861345  2.6332797
#> [4,] -2.0181850 -2.0795939  1.3258059         NA  3.0115924  0.4798541
#> [5,]  0.1936413 -0.4947831 -1.0861345  3.0115924         NA -1.2899547
#> [6,]  3.8587216  6.1625534  2.6332797  0.4798541 -1.2899547         NA
#> [7,]  0.5671620  0.5877623  5.8430378 -0.9769979  1.5757127 -1.0129338
#> [8,]  0.3126387 -3.4715041 -0.9287578  4.0101249  0.7587521  2.3224311
#>            [,7]       [,8]
#> [1,]  0.5671620  0.3126387
#> [2,]  0.5877623 -3.4715041
#> [3,]  5.8430378 -0.9287578
#> [4,] -0.9769979  4.0101249
#> [5,]  1.5757127  0.7587521
#> [6,] -1.0129338  2.3224311
#> [7,]         NA  2.2571849
#> [8,]  2.2571849         NA
#> 
#> $Mmat_p
#>              [,1]         [,2]         [,3]         [,4]        [,5]
#> [1,]           NA 5.094700e-01 1.884374e-01 4.357200e-02 0.846456748
#> [2,] 0.5094699872           NA 3.402551e-01 3.756280e-02 0.620753227
#> [3,] 0.1884374472 3.402551e-01           NA 1.849040e-01 0.277419500
#> [4,] 0.0435719956 3.756280e-02 1.849040e-01           NA 0.002598813
#> [5,] 0.8464567475 6.207532e-01 2.774195e-01 2.598813e-03          NA
#> [6,] 0.0001139817 7.158116e-10 8.456467e-03 6.313312e-01 0.197066396
#> [7,] 0.5706041434 5.566918e-01 5.125738e-09 3.285702e-01 0.115092039
#> [8,] 0.7545551426 5.175514e-04 3.530146e-01 6.068665e-05 0.448000891
#>              [,6]         [,7]         [,8]
#> [1,] 1.139817e-04 5.706041e-01 7.545551e-01
#> [2,] 7.158116e-10 5.566918e-01 5.175514e-04
#> [3,] 8.456467e-03 5.125738e-09 3.530146e-01
#> [4,] 6.313312e-01 3.285702e-01 6.068665e-05
#> [5,] 1.970664e-01 1.150920e-01 4.480009e-01
#> [6,]           NA 3.110918e-01 2.020973e-02
#> [7,] 3.110918e-01           NA 2.399652e-02
#> [8,] 2.020973e-02 2.399652e-02           NA
```

For more examples, see individual help files.

## Technical references

Please refer to the following sources regarding various facets of the
FVBM models that are implemented in the package.

The FVBM model and the consistency of their maximum pseudolikelihood
estimators (MPLEs) was first considered in
<http://doi.org/10.1162/neco.2006.18.10.2283>. The MM algorithm
implemented in the main function `fitfvbm` was introduced in
<http://doi.org/10.1162/NECO_a_00813>. Here various convergence results
regarding the algorithm is proved. Next, the asymptotic normality
results pertaining to the use of the functions `fvbmstderr` and
`fvbmtests` are proved in <http://doi.org/10.1109/TNNLS.2015.2425898>.
Finally, the `senate` data was introduced and analysed in
<https://hal.archives-ouvertes.fr/hal-01927188v1>.

## Reference to package

If you find this package useful in your work, then please follow the
usual `R` instructions for citing the package in your publications. That
is, follow the instructions from `citation('BoltzMM')`.

``` r
# Citation instructions
citation('BoltzMM')
#> Warning in citation("BoltzMM"): no date field in DESCRIPTION file of
#> package 'BoltzMM'
#> Warning in citation("BoltzMM"): could not determine year for 'BoltzMM' from
#> package DESCRIPTION file
#> 
#> To cite package 'BoltzMM' in publications use:
#> 
#>   Andrew Thomas Jones, Hien Duy Nguyen and Jessica Juanita Bagnall
#>   (NA). BoltzMM: Boltzmann Machines with MM Algorithms. R package
#>   version 0.1.3.
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Manual{,
#>     title = {BoltzMM: Boltzmann Machines with MM Algorithms},
#>     author = {Andrew Thomas Jones and Hien Duy Nguyen and Jessica Juanita Bagnall},
#>     note = {R package version 0.1.3},
#>   }
#> 
#> ATTENTION: This citation information has been auto-generated from
#> the package DESCRIPTION file and may need manual editing, see
#> 'help("citation")'.
```

## Authorship statement

The `BoltzMM` package is co-authored by [Andrew T.
Jones](https://github.com/andrewthomasjones), [Hien D.
Nguyen](https://github.com/hiendn), and Jessica J. Bagnall. The initial
development of the package, in native `R` was conducted by HDN.
Implementation of the core loops of the package in the `C` language was
performed by ATJ. JJB formatted and contributed the `senate` data set as
well as the example analysis on the `senate` data. All three co-authors
contributed to the documentation of the software as well as
troubleshooting and testing.

## Unit testing

Using the package `testthat`, we have conducted the following unit test
for the GitHub build, on the date: 31 January, 2019. The testing files
are contained in the
[tests](https://github.com/andrewthomasjones/BoltzMM/tree/master/tests)
folder of the repository.

## Bug reporting and contributions

Thank you for your interest in `BoltzMM`. If you happen to find any bugs
in the program, then please report them on the Issues page
(<https://github.com/andrewthomasjones/BoltzMM/issues>). Support can
also be sought on this page. Furthermore, if you would like to make a
contribution to the software, then please forward a pull request to the
owner of the repository.
