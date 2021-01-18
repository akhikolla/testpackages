## ----setup, echo=FALSE---------------------------------------------------
has_kdtools = require(kdtools, quietly = TRUE)
if (!has_kdtools) cat("kdtools package not installed; results not evaluted\n")
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  eval = has_kdtools
)

## ------------------------------------------------------------------------
library(kdtools)
x = matrix(runif(3e3), nc = 3)
y = matrix_to_tuples(x)
y[1:3, c(1, 3)]

## ------------------------------------------------------------------------
kd_sort(y, inplace = TRUE, parallel = TRUE)

## ------------------------------------------------------------------------
rq = kd_range_query(y, c(0, 0, 0), c(1/4, 1/4, 1/4)); rq
i = kd_nearest_neighbor(y, c(0, 0, 0)); y[i, ]
nns = kd_nearest_neighbors(y, c(0, 0, 0), 100); nns
nni = kd_nn_indices(y, c(0, 0, 0), 10); nni

## ------------------------------------------------------------------------
head(tuples_to_matrix(rq))
head(tuples_to_matrix(nns))

