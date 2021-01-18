context("MOI")

test_that("MOI works for (1, 2, 3) and k=1", {
  expect_equal(MOI(1:3, k=1), sqrt(3)*4/9)
})

test_that("MOI works for rnorm(50) and k=1", {
  set.seed(1)
  X <- rnorm(50)
  expect_equal(MOI(X, k=1), sqrt(50)*0.07266939, tolerance = 1e-7)
})

test_that("MOI works for rnorm(50) and k=2", {
  set.seed(1)
  X <- rnorm(50)
  expect_equal(MOI(X, k=2), sqrt(50)*0.08693539, tolerance = 1e-7)
})
