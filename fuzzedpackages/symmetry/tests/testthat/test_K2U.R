context("K2U")

test_that("K2U works for (1, 2)", {
  expect_equal(K2U(1:2), 2*1)
})

test_that("K2U works for rnorm(50)", {
  set.seed(1)
  X <- rnorm(50)
  expect_equal(K2U(X), 50*0.03265306, tolerance = 1e-5)
})

