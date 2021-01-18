
<!-- README.md is generated from README.Rmd. Please edit that file -->
RcppGreedySetCover
==================

A fast implementation of the greedy algorithm for the set cover problem using 'Rcpp'.

Installation
------------

You can install RcppGreedySetCover from github with:

``` r
# install.packages("devtools")
devtools::install_github("matthiaskaeding/RcppGreedySetCover")
```

Example
-------

This is a basic example which shows you how to use the main function:

``` r
# Create some data.
set.seed(333)
X <- data.table::rbindlist(
 lapply(
   seq_len(1e4L),
   function(x) list(element=sample.int(n=1e3L,size=sample.int(50L,1L)))
 ),
 idcol="set"
)

# Input is in long format.
head(X) 
#>    set element
#> 1:   1      85
#> 2:   1     973
#> 3:   1     571
#> 4:   1      21
#> 5:   1     721
#> 6:   1     607

# Run set cover
res <- RcppGreedySetCover::greedySetCover(X,FALSE)
#> 100% covered by 45 sets.

# Result is in long format.
head(res) 
#>   set element
#> 1  19      25
#> 2  19      65
#> 3  19      68
#> 4  19      86
#> 5  19      90
#> 6  19     105

# Check if all elements are covered.
identical(sort(unique(res$element)),sort(unique(X$element)))
#> [1] TRUE
```
