context("MOK")

test_that("MOK works for (1, 2, 3) and k=1", {
  expect_equal(MOK(1:3, 1), sqrt(3)*2/3)
})

test_that("MOK works for (1, 2, 3, 4) and k=1", {
  expect_equal(MOK(1:4, 1), 2*2/3)
})

test_that("MOK works for rnorm(50) and k=1", {
  set.seed(1)
  X <- rnorm(50)
  expect_equal(MOK(X, 1), sqrt(50)*0.1771429, tolerance = 1e-6)
})

test_that("MOK works for rnorm(50) and k=2", {
  set.seed(1)
  X <- rnorm(50)
  expect_equal(MOK(X, 2), sqrt(50)*0.1928571, tolerance = 1e-6)
})
