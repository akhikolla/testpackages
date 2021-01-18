
GCSM
====

<!-- badges: start -->

[![R-CMD-check](https://github.com/liuyadong/GCSM/workflows/R-CMD-check/badge.svg)](https://github.com/liuyadong/GCSM/actions)
<!-- badges: end -->

The goal of GCSM is to implement the generic composite similarity
measure (GCSM), described in “A generic composite measure of similarity
between geospatial variables” by Liu et al. (2020)
[doi:10.1016/j.ecoinf.2020.101169](https://doi.org/10.1016/j.ecoinf.2020.101169).
This package also provides implements of SSIM and CMSC. Functions are
given to compute composite similarity between vectors (e.g, `gcsm`), on
spatial windows (e.g., `gcsm_sw`) or temporal windows (e.g., `gcsm_tw`).
They are implemented in C++ with
[RcppArmadillo](https://github.com/RcppCore/RcppArmadillo). OpenMP is
used for parallel computing.

Installation
------------

You can install the package from [GitHub](https://github.com/) with:

    # install.packages("devtools")
    devtools::install_github("liuyadong/GCSM")

Examples
--------

Composite similarity between vectors:

    library(GCSM)

    x = runif(9)
    gcsm(x, x)
    #> [1] 1
    cmsc(x, x)
    #> [1] 1

    # mean shift
    gcsm(x, x - 0.2, xmin = 0, xmax = 1, ymin = 0, ymax = 1)
    #> [1] 0.8
    cmsc(x, x - 0.2, xmin = 0, xmax = 1, ymin = 0, ymax = 1)
    #> [1] 0.96
    gcsm(x, x + 0.2, xmin = 0, xmax = 1, ymin = 0, ymax = 1)
    #> [1] 0.8
    cmsc(x, x + 0.2, xmin = 0, xmax = 1, ymin = 0, ymax = 1)
    #> [1] 0.96
    ## dissimilarity
    y = 1 - x # y is the perfect antianalog of x
    gcsm(y, x)
    #> [1] -1
    gcsm(y, x - 0.2, xmin = 0, xmax = 1, ymin = 0, ymax = 1)
    #> [1] -0.8
    gcsm(y, x + 0.2, xmin = 0, xmax = 1, ymin = 0, ymax = 1)
    #> [1] -0.8

    # random noise
    noise = rnorm(9, mean = 0, sd = 0.1)
    gcsm(x, x + noise, xmin = 0, xmax = 1, ymin = 0, ymax = 1)
    #> [1] 0.9201221
    cmsc(x, x + noise, xmin = 0, xmax = 1, ymin = 0, ymax = 1)
    #> [1] 0.9416416
    ## dissimilariry
    gcsm(y, x + noise, xmin = 0, xmax = 1, ymin = 0, ymax = 1)
    #> [1] -0.9201221

Composite similarity on spatial windows:

    x = matrix(runif(36), nrow = 6, ncol = 6)
    gcsm_sw(x, x + 0.2, xmin = 0, xmax = 1, ymin = 0, ymax = 1, ksize = 3)
    #>      [,1] [,2] [,3] [,4] [,5] [,6]
    #> [1,]  0.8  0.8  0.8  0.8  0.8  0.8
    #> [2,]  0.8  0.8  0.8  0.8  0.8  0.8
    #> [3,]  0.8  0.8  0.8  0.8  0.8  0.8
    #> [4,]  0.8  0.8  0.8  0.8  0.8  0.8
    #> [5,]  0.8  0.8  0.8  0.8  0.8  0.8
    #> [6,]  0.8  0.8  0.8  0.8  0.8  0.8
    cmsc_sw(x, x + 0.2, xmin = 0, xmax = 1, ymin = 0, ymax = 1, ksize = 3)
    #>      [,1] [,2] [,3] [,4] [,5] [,6]
    #> [1,] 0.96 0.96 0.96 0.96 0.96 0.96
    #> [2,] 0.96 0.96 0.96 0.96 0.96 0.96
    #> [3,] 0.96 0.96 0.96 0.96 0.96 0.96
    #> [4,] 0.96 0.96 0.96 0.96 0.96 0.96
    #> [5,] 0.96 0.96 0.96 0.96 0.96 0.96
    #> [6,] 0.96 0.96 0.96 0.96 0.96 0.96
    ssim_sw(x, x + 0.2, xmin = 0, xmax = 1, ymin = 0, ymax = 1, ksize = 3)
    #>           [,1]      [,2]      [,3]      [,4]      [,5]      [,6]
    #> [1,] 0.9405428 0.9107526 0.8956004 0.8758824 0.8983908 0.8752976
    #> [2,] 0.9356499 0.9213593 0.9306332 0.9179906 0.9268518 0.9082596
    #> [3,] 0.9266229 0.9361720 0.9497504 0.9331137 0.9312823 0.9243788
    #> [4,] 0.9044219 0.9205696 0.9334963 0.9157745 0.9144879 0.9159464
    #> [5,] 0.9411510 0.9265003 0.9171057 0.9065103 0.9271306 0.9304926
    #> [6,] 0.9580466 0.9272437 0.9095319 0.9179363 0.9485734 0.9454656

Composite similarity on temporal windows:

    x = array(runif(81), dim = c(3, 3, 9))
    gcsm_tw(x, x + 0.2, xmin = 0, xmax = 1, ymin = 0, ymax = 1)
    #>      [,1] [,2] [,3]
    #> [1,]  0.8  0.8  0.8
    #> [2,]  0.8  0.8  0.8
    #> [3,]  0.8  0.8  0.8
    cmsc_tw(x, x + 0.2, xmin = 0, xmax = 1, ymin = 0, ymax = 1)
    #>      [,1] [,2] [,3]
    #> [1,] 0.96 0.96 0.96
    #> [2,] 0.96 0.96 0.96
    #> [3,] 0.96 0.96 0.96
