context("K2")

test_that("K2 works for (1, 2)", {
  expect_equal(K2(1:2), 2)
})

test_that("K2 works for rnorm(50)", {
  set.seed(1)
  X <- rnorm(50)
  expect_equal(K2(X), 50*0.044, tolerance = 1e-5)
})

