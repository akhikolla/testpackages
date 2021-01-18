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
calcEntropyBase <- function(values, nobs, s) {
  -sum(values[values>0] * log(values[values>0], base = 2))
}
### abellan correction
calcEntropyAbellan <- function(values, nobs, s) {
  calcEntropyBase(values, nobs, s) + 
    (s * log(length(values), base = 2)) / (nobs + s)
}
### strobl correction
calcEntropyStrobl <- function(values, nobs, s) {
  calcEntropyBase(values, nobs, s) + 
    (length(values) + 1) / (2 * nobs + s)
}

## test template #####################################################
value_tests <- function(pI, values, s, entropyfun) {
  expect_equal(pI$probint[1,],
               values[[1]])
  expect_equal(pI$probint[2,],
               values[[1]]  / (sum(values[[1]]) + s))
  expect_equal(pI$probint[3,],
               (values[[1]] + s) / (sum(values[[1]]) + s))
  expect_equal(pI$maxEntDist,
               values[[2]])
  expect_equal(pI$maxEntCorr,
               entropyfun(values[[2]], sum(values[[1]]), s))
  expect_equal(pI$minEntDist,
               values[[3]])
  expect_equal(pI$minEntCorr,
               entropyfun(values[[3]], sum(values[[1]]), s))
}

# loading library
library(imptree)

#----------------- Test Suite starts here! --------------------------#
ip <- "IDM"

## IDM entropy calculation ###########################################
context("IDM entropy calculation")

test_that("IDM, nonsense 's'", {
  vec <- c(0,0,1,1)
  
  expect_error(probInterval(vec, ip, correction = "no", s = -2, 
                            entropymin = TRUE, entropymax = TRUE),
               "value of 's' \\(-2.000000\\) must be strictly positive")
})

test_that("IDM, direct", {
  vec <- c(0,0,1,1)
  pI <- probInterval(vec, ip, correction = "no", s = 2, 
                     entropymin = TRUE, entropymax = TRUE)
  expected_vals <- list(vec, 
                        rep(0.25, 4), 
                        c(0, 0, 0.75, 0.25))
  
  value_tests(pI, expected_vals, 2, calcEntropyBase)
})

test_that("IDM, more iterations", {
  vec <- c(2,3,1,2,0)
  pI <- probInterval(vec, ip, correction = "no", s = 2, 
                     entropymin = TRUE, entropymax = TRUE)
  expected_vals <- list(vec,
                        c(0.2, 0.3, 0.15, 0.2, 0.15), 
                        c(0.2, 0.5, 0.1, 0.2, 0))
  
  value_tests(pI, expected_vals, 2, calcEntropyBase)
})


## IDM entropy correction methods ####################################
context("IDM entropy correction")

test_that("IDM, no correction", {
  vec <- c(2,3,1,2,0)
  pI <- probInterval(vec, ip, correction = "no", s = 2, 
                     entropymin = TRUE, entropymax = TRUE)
  expected_vals <- list(vec,
                        c(0.2, 0.3, 0.15, 0.2, 0.15), 
                        c(0.2, 0.5, 0.1, 0.2, 0))
  
  value_tests(pI, expected_vals, 2, calcEntropyBase)
})

test_that("IDM, strobl correction", {
  vec <- c(2,3,1,2,0)
  pI <- probInterval(vec, ip, correction = "strobl", s = 2, 
                     entropymin = TRUE, entropymax = TRUE)
  expected_vals <- list(vec,
                        c(0.2, 0.3, 0.15, 0.2, 0.15), 
                        c(0.2, 0.5, 0.1, 0.2, 0))
  
  value_tests(pI, expected_vals, 2, calcEntropyStrobl)
})

test_that("IDM, abellan correction", {
  vec <- c(2,3,1,2,0)
  pI <- probInterval(vec, ip, correction = "abellan", s = 2, 
                     entropymin = TRUE, entropymax = TRUE)
  expected_vals <- list(vec,
                        c(0.2, 0.3, 0.15, 0.2, 0.15), 
                        c(0.2, 0.5, 0.1, 0.2, 0))
  
  value_tests(pI, expected_vals, 2, calcEntropyAbellan)
})

test_that("IDM, nonsense correction", {
  vec <- c(2,3,1,2,0)
  
  expect_error(probInterval(vec, ip, correction = "bla", s = 2, 
                            entropymin = TRUE, entropymax = TRUE))
})
