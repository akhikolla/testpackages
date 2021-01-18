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

# loading library and data
library(imptree)
data("carEvaluation")

#----------------- Test Suite starts here! --------------------------#

## imprecise tree (stub)
ip <- imptree(acceptance~., data = carEvaluation, 
              method = "NPI", 
              method.param = list(splitmetric = "globalmax"), 
              control = list(depth = 1, minbucket = 1))

## Imptree (stub) ####################################################
context("Imptree (stub)")

test_that("imptree", {
  expect_s3_class(ip,"imptree")
  expect_length(ip, 4)
  expect_named(ip, c("call", "tree", "traindata", "formula"))
  expect_true(imptree:::hasRoot_cpp(ip$tree))
    
})

test_that("imptree.summary", {
  summimp <- summary(ip)
  expect_is(summimp, "summary.imptree")
  expect_s3_class(summimp, "summary.imptree")
  expect_length(summimp, 6)
  expect_named(summimp, c("call", "utility", "dominance",
                          "sizes", "acc", "meta"))
  
  expect_equal(summimp$utility, 0.65)
  
  expect_equal(summimp$dominance, "strong")
  
  expect_equal(summimp$sizes$obs, 1728)
  expect_equal(summimp$sizes$iobs, 0L)
  
  expect_equal(summimp$meta,
               c("depth" = 1, "nleaves" = 3, "nnodes" = 4))
  
  expect_equal(summimp$acc, c("Determinacy" = 1,
                              "Average indeterminate size" = NA_real_,
                              "Single-Set Accuracy" = 1210/1728,
                              "Set-Accuracy" = NA_real_,
                              "Discounted Accuracy" = 1210/1728,
                              "0.65 utility based Accuracy" = 1210/1728))
})

test_that("root node", {
  expect_message(rn <- node_imptree(ip, NULL),
                 "extracting probability information from root node")
  expect_is(rn, "node_imptree")
  expect_true(is.list(rn))
  expect_length(rn, 6L)
  expect_named(rn, c("probint", "depth", "splitter","children",
                     "traindataIdx", "ipmodel"))
  expect_equal(rn$depth, 0L)
  expect_equal(rn$splitter, "safety")
  expect_equal(rn$children, 3L)
  expect_length(rn$traindataIdx, 1728)
  expect_equal(range(rn$traindataIdx), c(1,1728))
  expect_length(rn$ipmodel, 1)
  expect_equal(rn$ipmodel$iptype, "NPI")
  expect_is(rn$probint, "matrix")
  expect_equal(dim(rn$probint), c(3,4))
  expect_equal(dimnames(rn$probint), 
               list(c("Frequency", "Lower", "Upper"), 
                    c("acc", "good", "unacc", "vgood")))
  expect_equivalent(rn$probint[1,], c(384, 69, 1210, 65))
  expect_true(all(rn$probint[2:3,] >= 0))
  expect_true(all(rn$probint[2:3,] <= 1))
})

test_that("child (index 1) below root", {
  cn <- node_imptree(ip, 3)
  expect_is(cn, "node_imptree")
  expect_true(is.list(cn))
  expect_length(cn, 6L)
  expect_named(cn, c("probint", "depth", "splitter","children",
                     "traindataIdx", "ipmodel"))
  expect_equal(cn$depth, 1L)
  expect_equal(cn$splitter, NA_character_)
  expect_equal(cn$children, 0L)
  
  expect_length(cn$traindataIdx, 576)
  expect_equal(range(cn$traindataIdx), c(2,1727))
  expect_length(cn$ipmodel, 1)
  expect_equal(cn$ipmodel$iptype, "NPI")
  expect_is(cn$probint, "matrix")
  expect_equal(dim(cn$probint), c(3,4))
  expect_equal(dimnames(cn$probint), 
               list(c("Frequency", "Lower", "Upper"), 
                    c("acc", "good", "unacc", "vgood")))
  expect_equivalent(cn$probint[1,], c(180, 39, 357, 0))
  expect_true(all(cn$probint[2:3,] >= 0))
  expect_true(all(cn$probint[2:3,] <= 1))
})

test_that("Out of range child query", {
  expect_error(node_imptree(ip, 4), 
               "Queried index (4) > child size (3)", fixed = TRUE)
})
