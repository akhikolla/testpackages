context("NAI")

test_that("NAI works for (1, 2, 3) and k=2", {
  expect_equal(NAI(1:3, 2), sqrt(3)*4/9)
})

test_that("NAI works for rnorm(50) and k=2", {
  set.seed(1)
  X <- rnorm(50)
  expect_equal(NAI(X, 2), sqrt(50)*0.07266939, tolerance = 1e-7)
})
