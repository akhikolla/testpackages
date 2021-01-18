context("test-ecf_mod")

test_that("multivariate sample at multiple vectors", {
  nn <- 500
  smp <- cbind(rnorm(nn), rnorm(nn))
  x <- matrix(c(0, 0, 0, 1), ncol = 2)
  out <- ecf_mod(x, smp)
  expect_true(length(out) == 2)
})


test_that("multivariate sample at single vector", {
  nn <- 500
  smp <- cbind(rnorm(nn), rnorm(nn))
  out <- ecf_mod(0:1, smp)
  expect_true(length(out) == 1)
})

test_that("univariate sample at multiple points", {
  out <- ecf_mod(seq(-2, 2, length.out = 100), rnorm(10))
  expect_true(length(out) == 100)
})


test_that("univariate sample at single point", {
  out <- ecf_mod(0, rnorm(10))
  expect_true(length(out) == 1)
})

test_that("positive number", {
  out <- ecf_mod(1, rnorm(10))
  expect_true(out >= 0)
})

test_that("error when input is incorrect", {
  expect_error(ecf_mod(0, "drumpf"))
})

test_that("error when dimensions are different", {
  expect_error(ecf_mod(matrix(1:12, 3, 4), matrix(1:12, 4, 3)))
})

test_that("error when an entry is an 3D-array", {
  expect_error(ecf_mod(array(1:24, 2:4), 1:3))
})
