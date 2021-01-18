context("BHI")

test_that("BHI works for (1, 2, 3)", {
  expect_equal(BHI(1:3), sqrt(3)*2/9, tolerance = 1e-8)
})

test_that("BHI works for (1, 2, 3, 4)", {
  expect_equal(BHI(1:4), 2*5/24, tolerance = 1e-8)
})

test_that("BHI works for rnorm(50)", {
  set.seed(1)
  X <- rnorm(50)
  expect_equal(BHI(X), sqrt(50)*0.0363346, tolerance = 1e-6)
})
