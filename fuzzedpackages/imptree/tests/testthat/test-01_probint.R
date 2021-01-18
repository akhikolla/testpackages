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

## test template #####################################################
structure_tests <- function(pI, hasmin, hasmax, nc) {
  expect_equal(length(pI), 1 + 2 * hasmax + 2 * hasmin)
  listnames <- c("probint", "maxEntDist", "maxEntCorr", "minEntDist",
                 "minEntCorr")[c(1, c(2,3) * hasmax, c(4,5) * hasmin)]

  expect_named(pI,  listnames)
  expect_equivalent(sapply(pI, length), 
                    c(3 * nc, rep(c(nc, 1), times = hasmin + hasmax)))
                                          
  expect_equal(rownames(pI$probint), c("Frequency", "Lower", "Upper"))
  expect_null(colnames(pI$probint))
  expect_equal(dim(pI$probint), c(3, nc))
}


# loading library
library(imptree)

#----------------- Test Suite starts here! --------------------------#

## some vector and iptype
ip <- "NPI"
vec <- c(0,0,1,1)

## Structure of probInterval #########################################
context("ProbInterval structure")

test_that("Probint, full", {
  pI <- probInterval(table = vec, iptype = ip, correction = "no",
                     entropymin = TRUE, entropymax = TRUE)
  structure_tests(pI, TRUE, TRUE, length(vec))
})

test_that("Probint, only max entropy", {
  pI <- probInterval(table = vec, iptype = ip, correction = "no",
                     entropymin = FALSE, entropymax = TRUE)
  structure_tests(pI, FALSE, TRUE, length(vec))
})

test_that("Probint, only min entropy", {
  pI <- probInterval(table = vec, iptype = ip, correction = "no",
                     entropymin = TRUE, entropymax = FALSE)
  structure_tests(pI, TRUE, FALSE, length(vec))
})

test_that("Probint, no entropy", {
  pI <- probInterval(table = vec, iptype = ip, correction = "no",
                     entropymin = FALSE, entropymax = FALSE)
  structure_tests(pI, FALSE, FALSE, length(vec))
})
