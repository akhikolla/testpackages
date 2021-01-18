context("BHK")

test_that("BHK works for (1, 2, 3)", {
  expect_equal(BHK(1:3), sqrt(3)*1/3, tolerance = 1e-8)
})

test_that("BHK works for (1, 2, 3, 4)", {
  expect_equal(BHK(1:4), 2*1/3, tolerance = 1e-8)
})

test_that("BHK works for rnorm(50)", {
  set.seed(1)
  X <- rnorm(50)
  expect_equal(BHK(X), 0.6262946, tolerance = 1e-6)
})
