library(testthat)
library(AnaCoDa)

context("CovarianceMatrix")

test_that("general covariance matrix functions", {
  expect_equal(testCovarianceMatrix(), 0)
})

covM <- new(CovarianceMatrix)

# TODO: Implement the following below.
#test_that("set Covariance Matrix", {
#  expect_equal(g$checkIndex(2, 1, 10), TRUE)
#  expect_equal(g$checkIndex(5, 1, 10), TRUE)
#  expect_equal(g$checkIndex(20, 1, 30), TRUE)
#  expect_equal(g$checkIndex(5, 4, 6), TRUE)
#  
#  #Checking invalid cases
#  expect_equal(g$checkIndex(11, 3, 10), FALSE)
#  expect_equal(g$checkIndex(5, 6, 11), FALSE)
#})
