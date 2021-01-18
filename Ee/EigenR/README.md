EigenR
================

Originally, I entitled this package *Fast Matrix Algebra with ‘Eigen’*,
because I expected it to be faster than R base. But this is not the
case. So I entitled it *Complex Matrix Algebra with ‘Eigen’*, because it
supports some operations on complex matrices which are not supported by
R base: determinant, Cholesky decomposition, and linear least-squares
problems.

``` r
library(EigenR)
library(microbenchmark)
```

## Determinant

``` r
set.seed(666L)
M <- matrix(rnorm(300L*300L, mean = 1), 300L, 300L)
M[sample.int(300L*300L, 300L*270L)] <- 0 # 90% of zeros
Ms <- asSparseMatrix(M)
microbenchmark(
  base          = det(M),
  EigenR        = Eigen_det(M),
  EigenR_sparse = Eigen_det(Ms), # :-(
  times = 200L
)
## Unit: milliseconds
##           expr       min        lq      mean    median        uq      max neval
##           base  2.298912  2.379902  3.852771  2.595700  3.204240 26.37411   200
##         EigenR  3.267192  6.587482  6.920492  6.776801  7.068001 22.67092   200
##  EigenR_sparse 11.928411 12.468541 12.970583 12.652409 12.924786 29.39698   200
##  cld
##  a  
##   b 
##    c
```

Determinants of complex matrices are supported:

``` r
set.seed(666L)
Mr <- matrix(rnorm(100L*100L, mean = 1), 100L, 100L)
Mi <- matrix(rnorm(100L*100L, mean = 1), 100L, 100L)
M <- Mr + 1i * Mi
library(complexplus)
microbenchmark(
  EigenR      = Eigen_det(M), # :-)
  complexplus = Det(M), 
  times = 30L
)
## Unit: milliseconds
##         expr       min        lq      mean    median        uq      max neval
##       EigenR  2.430008  2.852084  2.990067  2.890091  3.207062  3.55512    30
##  complexplus 35.228304 36.151367 37.867418 36.983333 38.792447 45.66115    30
##  cld
##   a 
##    b
```

## Cholesky decomposition

``` r
set.seed(666L)
M <- matrix(rpois(10000L, 25), 100L, 100L)
A <- crossprod(M)
microbenchmark(
  base   = chol(A),
  EigenR = Eigen_chol(A), # :-(
  times = 1000L
)
## Unit: microseconds
##    expr     min      lq     mean   median       uq       max neval cld
##    base 148.137 161.859 429.0553 170.6695 216.3325 24262.700  1000   a
##  EigenR 149.939 330.668 395.4294 372.3930 434.6495  9127.298  1000   a
```

Cholesky decomposition of complex matrices is supported.

## Pivoted Cholesky decomposition

``` r
set.seed(666L)
M <- matrix(rgamma(202L*199L, 10), 202L, 199L)
M <- cbind(M, M[, 1L] + 3*M[, 2L])
A <- crossprod(M)
microbenchmark(
  base   = chol(A, pivot = TRUE),
  EigenR = Eigen_UtDU(A), # :-(
  times = 1000L
)
## Unit: milliseconds
##    expr      min       lq     mean   median       uq       max neval cld
##    base 1.118317 1.175779 1.678367 1.269667 1.471796 17.062454  1000  a 
##  EigenR 1.412606 1.909860 2.074614 1.946364 2.074476  9.196924  1000   b
```

Pivoted Cholesky decomposition of complex matrices is supported.

## Kernel

``` r
set.seed(666L)
M <- matrix(rpois(10000L, 25), 100L, 100L)
A <- cbind(M, M)
At <- t(A)
library(MASS)
microbenchmark(
  MASS       = Null(At),
  EigenR_LU  = Eigen_kernel(A, method = "LU"),  # :-)
  EigenR_COD = Eigen_kernel(A, method = "COD"), 
  times = 100L
)
## Unit: milliseconds
##        expr      min       lq     mean   median       uq        max neval cld
##        MASS 4.658614 4.847165 7.429655 5.005544 5.433269 118.029181   100   b
##   EigenR_LU 2.162524 2.237471 2.403308 2.308238 2.398000   5.455091   100  a 
##  EigenR_COD 5.046402 5.288316 5.691691 5.550173 5.911409   8.938171   100   b
```

## Linear least-squares problems

``` r
set.seed(666L)
n <- 700L; p <- 200L
A <- matrix(rnorm(n * p), n, p)
b <- rnorm(n)
microbenchmark(
  stats   = lm.fit(A, b),
  Eigen_R = Eigen_lsSolve(A, b), # :-(
  times = 20L
)
## Unit: milliseconds
##     expr      min       lq     mean   median       uq      max neval cld
##    stats 14.41073 14.83985 15.23161 15.17467 15.49191 16.64085    20  a 
##  Eigen_R 73.97197 74.77006 75.77076 75.45565 76.82308 78.93335    20   b
```

Complex matrices `A` and `b` are supported.
