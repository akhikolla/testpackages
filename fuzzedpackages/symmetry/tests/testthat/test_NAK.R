context("NAK")

test_that("NAK works for (1, 2, 3) and k=2", {
  expect_equal(NAK(1:3, 2), sqrt(3)*2/3)
})

test_that("NAK works for (1, 2, 3, 4) and k=2", {
  expect_equal(NAK(1:4, 2), 2*2/3)
})

test_that("NAK works for rnorm(50) and k=2", {
  set.seed(1)
  X <- rnorm(50)
  expect_equal(NAK(X, 2), sqrt(50)*0.1771429, tolerance = 1e-6)
})
