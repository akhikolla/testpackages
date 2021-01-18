## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

library("samc")
library("raster")

res_data <- samc::ex_res_data
abs_data <- samc::ex_abs_data
occ_data <- samc::ex_occ_data

## ---- error = TRUE------------------------------------------------------------
# Working example
r1 <- ex_res_data
r2 <- ex_res_data
check(r1, r2)


# Remove the NA's in r2 by overwriting all the elements with the number 1.
# check() doesn't check the actual values of the data, but it does check the
# location of NA's
r1 <- ex_res_data
r2 <- ex_res_data
r2[1] <- 1
check(r1, r2)


# Change the dimensions of r2 by subsetting it. check() ensures that the data
# inputs have the same number of rows and columns
r1 <- ex_res_data
r2 <- ex_res_data
r2 <- r2[1:5, 1:5]
check(r1, r2)

## ---- error = TRUE------------------------------------------------------------
# Example: Skipping optional arguments. In this case, the `fidelity` argument is
# optional and the `latlon` argument doesn't apply because our input data are
# matrices, so we skip them. The tr_fun argument, however, is always required, 
# so we pass `function(x) 1/mean(x)` to it. But because we don't specify which
# argument it is, R is trying to find a version of the function that expects it
# as the third argument, but this version does not exist.
samc_obj <- samc(res_data, abs_data, function(x) 1/mean(x))

# Solution
samc_obj <- samc(res_data, abs_data, tr_fun = function(x) 1/mean(x))


# Example: Incorrect input types. In this case, we are attempting to pass a
# single numeric value as absorption data. However, the absorption data must 
# always be in a matrix or RasterLayer object
samc_obj <- samc(res_data, 0.01, tr_fun = function(x) 1/mean(x))

# Solution
samc_obj <- samc(res_data, abs_data, tr_fun = function(x) 1/mean(x))

