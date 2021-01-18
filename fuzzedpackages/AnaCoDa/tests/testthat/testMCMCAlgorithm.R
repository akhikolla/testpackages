library(testthat)
library(AnaCoDa)

context("MCMCAlgorithm")

test_that("general MCMCAlgorithm functions", {
  expect_equal(testMCMCAlgorithm(), 0)
})

samples <- 1000
thinning <- 1
adaptiveWidth <- 100
mcmc <- new(MCMCAlgorithm, samples, thinning, adaptiveWidth, TRUE, TRUE, TRUE)

test_that("get Samples", {
  expect_equal(mcmc$getSamples(), samples)
})

test_that("get Thinning", {
  expect_equal(mcmc$getThinning(), thinning)
})

test_that("get Adaptive Width", {
  expect_equal(mcmc$getAdaptiveWidth(), adaptiveWidth)
})

test_that("set Samples", {
  mcmc$setSamples(10)
  expect_equal(mcmc$getSamples(), 10)
})

test_that("set Thinning", {
  mcmc$setThinning(10)
  expect_equal(mcmc$getThinning(), 10)
})

test_that("set Adaptive Width", {
  mcmc$setAdaptiveWidth(10)
  expect_equal(mcmc$getAdaptiveWidth(), 10)
})

test_that("set Log Posterior Trace", {
  vect <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)
  mcmc$setLogPosteriorTrace(vect)
  expect_equal(mcmc$getLogPosteriorTrace(), vect)
})
