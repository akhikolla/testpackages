# Copyright (C) 2018  Paul Fink
#
# This file is part of imptree.
#
# imptree is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# imptree is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with imptree.  If not, see <https://www.gnu.org/licenses/>.

#------------------------ Preparations ------------------------------#

## entropy calculations ##############################################
### no correction
calcEntropyBase <- function(values, nobs) {
  -sum(values[values>0] * log(values[values>0], base = 2))
}
### strobl correction
calcEntropyStrobl <- function(values, nobs) {
  ent <- -sum(values[values>0] * log(values[values>0], base = 2))
  ent + (length(values) - 1) / (2 * nobs)
}

## test template #####################################################
value_tests <- function(pI, values, entropyfun) {
  expect_equal(pI$probint[1,],
               values[[1]])
  expect_equal(pI$probint[2,],
               pmax(0, (values[[1]] - 1) / sum(values[[1]])))
  expect_equal(pI$probint[3,],
               pmin(1, (values[[1]] + 1) / sum(values[[1]])))
  expect_equal(pI$maxEntDist,
               values[[2]])
  expect_equal(pI$maxEntCorr,
               entropyfun(values[[2]], sum(values[[1]])))
  expect_equal(pI$minEntDist,
               values[[3]])
  expect_equal(pI$minEntCorr,
               entropyfun(values[[3]], sum(values[[1]])))
}

# loading library
library(imptree)

#----------------- Test Suite starts here! --------------------------#

## NPI with exact max entroy calculation #############################
context("NPI entropy calculation exact")

ip <- "NPI"

test_that("Exact NPI, krem=0", {
  vec <- c(0,0,1,1)
  pI <- probInterval(table = vec, iptype = ip, correction = "no", 
                     entropymin = TRUE, entropymax = TRUE)
  expected_vals <- list(vec, 
                        rep(0.25, 4), 
                        c(0.5, 0.5, 0, 0))
  
  value_tests(pI, expected_vals, calcEntropyBase)
})

test_that("Exact NPI, krem < k0; h < k1+1", {
  vec <- c(1,0,0,4)
  pI <- probInterval(table = vec, iptype = ip, correction = "no", 
                     entropymin = TRUE, entropymax = TRUE)
  expected_vals <- list(vec, 
                        c(rep(4/30, 3), 0.6),
                        c(0, 0, 0, 1))
  
  value_tests(pI, expected_vals, calcEntropyBase)
})

test_that("Exact NPI, krem < k0; h > k1+1", {
  vec <- c(rep(0,6), 1:4)
  pI <- probInterval(table = vec, iptype = ip, correction = "no", 
                     entropymin = TRUE, entropymax = TRUE)
  expected_vals <- list(vec, 
                        c(rep(0.05, 6), 0.1, 0.1, 0.2, 0.3),
                        c(rep(0, 7), 0.1, 0.4, 0.5))
  
  value_tests(pI, expected_vals, calcEntropyBase)
})

test_that("Exact NPI, krem >= k0", {
  vec <- c(3,3,2,3,4,5,0,0)
  pI <- probInterval(table = vec, iptype = ip, correction = "no", 
                     entropymin = TRUE, entropymax = TRUE)
  expected_vals <- list(vec,
                        c(rep(0.1375,4), 0.15, 0.2, 0.05, 0.05),
                        c(0.2, 0.1, 0.05, 0.1, 0.25, 0.3, 0, 0))
  
  value_tests(pI, expected_vals, calcEntropyBase)
})

## NPI with approximate max entropy calculation ######################

context("NPI entropy calculation approx")

ip <- "NPIapprox"

test_that("Approx NPI, krem < k0", {
  vec <- c(0,0,1,3)
  pI <- probInterval(table = vec, iptype = ip, correction = "no", 
                     entropymin = TRUE, entropymax = TRUE)
  expected_vals <- list(vec, 
                        c(rep(1/6, 3), 0.5), 
                        c(rep(0, 3), 1))
  
  value_tests(pI, expected_vals, calcEntropyBase)
})

test_that("Approx NPI, krem >= k0", {
  vec <- c(3,3,2,3,4,5,0,0)
  pI <- probInterval(vec, ip, correction = "no", 
                     entropymin = TRUE, entropymax = TRUE)
  expected_vals <- list(vec,
                        c(rep(0.1375,4), 0.15, 0.2, 0.05, 0.05),
                        c(0.2, 0.1, 0.05, 0.1, 0.25, 0.3, 0, 0))
  
  value_tests(pI, expected_vals, calcEntropyBase)
})


## NPI entropy correction methods ####################################
context("NPI entropy correction")

test_that("Approx NPI, no correction", {
  vec <- c(3,3,2,3,4,5,0,0)
  pI <- probInterval(table = vec, iptype = ip, correction = "no", 
                     entropymin = TRUE, entropymax = TRUE)
  expected_vals <- list(vec,
                        c(rep(0.1375,4), 0.15, 0.2, 0.05, 0.05),
                        c(0.2, 0.1, 0.05, 0.1, 0.25, 0.3, 0, 0))
  
  value_tests(pI, expected_vals, calcEntropyBase)
})

test_that("Approx NPI, strobl correction", {
  vec <- c(3,3,2,3,4,5,0,0)
  pI <- probInterval(table = vec, iptype = ip, correction = "strobl", 
                     entropymin = TRUE, entropymax = TRUE)
  expected_vals <- list(vec,
                        c(rep(0.1375,4), 0.15, 0.2, 0.05, 0.05),
                        c(0.2, 0.1, 0.05, 0.1, 0.25, 0.3, 0, 0))
  
  value_tests(pI, expected_vals, calcEntropyStrobl)
})

test_that("Approx NPI, abellan correction", {
  vec <- c(3,3,2,3,4,5,0,0)
  
  expect_error(
    probInterval(table = vec, iptype = ip, correction = "abellan", 
                 entropymin = TRUE, entropymax = TRUE), 
               sprintf("'correction' should be one of %s",
                       paste(dQuote(c("no","strobl")),
                             collapse = ", ")))
})

