context("test-obsquantiles.R")

test_that("obsquantiles works", {
  set.seed(123)
  tt <- seq(0, 5, .01)
  pars <- c(.8, 2, .5, .5, .5, # condition 1
           .8, 3, .5, .5, .5)  # condition 3
  pdfND <- dbeta(tt, 10, 30)
  data <- simData(n = 1e3, pars = pars, tt = tt, pdfND = pdfND)
  probs <- seq(0, 1, .2)
  q <- obsQuantiles(data, probs = probs)
  expect_equal(
    object = q,
    expected = structure(
      c(0.16, 0.24, 0.3, 0.33, 0.42, 0.87, 0.14, 0.28, 0.322,
        0.38, 0.46, 0.9, 0.17, 0.26, 0.3, 0.33, 0.42, 0.63, 0.13, 0.26,
        0.31, 0.36, 0.43, 0.94), .Dim = c(6L, 4L),
      .Dimnames = list(c("0%", "20%", "40%", "60%", "80%", "100%"),
                       c("lower.1", "upper.1", "lower.2", "upper.2"))),
    tol = 1e-1

  )
})
