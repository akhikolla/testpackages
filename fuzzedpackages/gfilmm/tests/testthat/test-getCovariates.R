test_that("getCovariates", {
  dat <- data.frame(
    x = 1:3,
    f = c("a", "a", "b"),
    g = c("x", "y", "y")
  )
  #
  cvrts <- getCovariates(dat, ~ 1)
  expect_length(cvrts$continuous, 0L)
  expect_length(cvrts$categorical, 0L)
  #
  cvrts <- getCovariates(dat, ~ x)
  expect_identical(cvrts$continuous, "x")
  expect_identical(unname(cvrts$categorical), list())
  #
  cvrts <- getCovariates(dat, ~ f)
  expect_identical(cvrts$continuous, character(0L))
  expect_identical(cvrts$categorical, list(f = c("a", "b")))
  #
  cvrts <- getCovariates(dat, ~ x*f)
  expect_identical(cvrts$continuous, "x")
  expect_identical(cvrts$categorical, list(f = c("a", "b")))
  #
  cvrts <- getCovariates(dat, ~ 0)
  expect_identical(cvrts$continuous, character(0L))
  expect_identical(unname(cvrts$categorical), list())
  #
  cvrts <- getCovariates(dat, ~ f*g)
  expect_identical(cvrts$continuous, character(0L))
  expect_identical(cvrts$categorical, list(f = c("a", "b"), g = c("x", "y")))
  #
  cvrts <- getCovariates(dat, ~ f*x)
  expect_identical(cvrts$continuous, "x")
  expect_identical(cvrts$categorical, list(f = c("a", "b")))
})
