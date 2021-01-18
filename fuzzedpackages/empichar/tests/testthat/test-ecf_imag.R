context("test-ecf_imag")

test_that("multivariate sample at multiple vectors", {
  nn <- 500
  smp <- cbind(rnorm(nn), rnorm(nn))
  x <- matrix(c(0, 0, 0, 1), ncol = 2)
  out <- ecf_imag(x, smp)
  expect_true(length(out) == 2)
})


test_that("multivariate sample at single vector", {
  nn <- 500
  smp <- cbind(rnorm(nn), rnorm(nn))
  out <- ecf_imag(0:1, smp)
  expect_true(length(out) == 1)
})

test_that("univariate sample at multiple points", {
  out <- ecf_imag(seq(-2, 2, length.out = 100), rnorm(10))
  expect_true(length(out) == 100)
})


test_that("univariate sample at single point", {
  out <- ecf_imag(0, rnorm(10))
  expect_true(length(out) == 1)
})

test_that("real number", {
  out <- ecf_imag(0, rnorm(10))
  expect_true(is.numeric(out))
})

test_that("error when input is incorrect", {
  expect_error(ecf_imag(0, "drumpf"))
})

test_that("error when dimensions are different", {
  expect_error(ecf_imag(matrix(1:12, 3, 4), matrix(1:12, 4, 3)))
})

test_that("error when an entry is an 3D-array", {
  expect_error(ecf_imag(array(1:24, 2:4), 1:3))
})
