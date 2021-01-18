context("GEV functions 2")

## Example shape parameters

# Positive
xi1 <- 0.1
low <- -1 / xi1

# Zero
xi2 <- 0

# Negative
xi3 <- -1e-7
up <- -1 / xi3

## Example input vectors
# For testing pdf, log_pdf and cdf
xvec <- c(-Inf, 0, Inf, NA)
x1 <- c(low, xvec)
x2 <- xvec
x3 <- c(up, xvec)
# For testing quantile
pvec <- c(0, 0.25, 0.5, 0.75, 1, NA)

test_that("rgev, length of output", {
  expect_length(rgev(n = 1, shape = xi1), 1)
  expect_length(rgev(n = 100, shape = xi1), 100)
  expect_length(rgev(n = 0, shape = xi1), 0)
  expect_error(rgev(n = -2, shape = xi1))

  expect_length(rgev(n = 1, shape = xi2), 1)
  expect_length(rgev(n = 100, shape = xi2), 100)
  expect_length(rgev(n = 0, shape = xi2), 0)
  expect_error(rgev(n = -2, shape = xi2))

  expect_length(rgev(n = 1, shape = xi3), 1)
  expect_length(rgev(n = 100, shape = xi3), 100)
  expect_length(rgev(n = 0, shape = xi3), 0)
  expect_error(rgev(n = -2, shape = xi3))
})

test_that("dgev, values and length of output", {
  p <- pvec[2:4]
  qq <- qgev(p, shape = xi1)
  expect_equal(dgev(x1, shape = xi1), c(0, 0, exp(-1), 0, NA))
  expect_equal(dgev(qq, shape = xi1), (-log(p)) ^ (1 + xi1) * p)
  expect_length(dgev(seq_len(0), shape = xi1), 0)
  expect_length(dgev(seq_len(1), shape = xi1), 1)
  expect_length(dgev(seq_len(10), shape = xi1), 10)

  qq <- qgev(p, shape = xi2)
  expect_equal(dgev(x1, shape = xi2), c(0, 0, exp(-1), 0, NA))
  expect_equal(dgev(qq, shape = xi2), (-log(p)) ^ (1 + xi2) * p)
  expect_length(dgev(seq_len(0), shape = xi2), 0)
  expect_length(dgev(seq_len(1), shape = xi2), 1)
  expect_length(dgev(seq_len(10), shape = xi2), 10)

  qq <- qgev(p, shape = xi3)
  expect_equal(dgev(x1, shape = xi3), c(0, 0, exp(-1), 0, NA))
  expect_equal(dgev(qq, shape = xi3), (-log(p)) ^ (1 + xi3) * p)
  expect_length(dgev(seq_len(0), shape = xi3), 0)
  expect_length(dgev(seq_len(1), shape = xi3), 1)
  expect_length(dgev(seq_len(10), shape = xi3), 10)
})

test_that("pgev, values and length of output", {
  expect_equal(pgev(x1, shape = xi1), c(0, 0, exp(-1), 1, NA))
  expect_length(pgev(seq_len(0), shape = xi1), 0)
  expect_length(pgev(seq_len(1), shape = xi1), 1)
  expect_length(pgev(seq_len(10), shape = xi1), 10)

  expect_equal(pgev(x2, shape = xi2), c(0, exp(-1), 1, NA))
  expect_length(pgev(seq_len(0), shape = xi2), 0)
  expect_length(pgev(seq_len(1), shape = xi2), 1)
  expect_length(pgev(seq_len(10), shape = xi2), 10)

  expect_equal(pgev(x3, shape = xi3), c(1, 0, exp(-1), 1, NA))
  expect_length(pgev(seq_len(0), shape = xi3), 0)
  expect_length(pgev(seq_len(1), shape = xi3), 1)
  expect_length(pgev(seq_len(10), shape = xi3), 10)
})

test_that("qgev, values and length of output", {
  q1 <- ((-log(pvec[2:4])) ^ (-xi1) - 1) / xi1
  expect_equal(qgev(pvec, shape = xi1), c(low, q1, Inf, NA))
  expect_length(qgev(seq_len(0), shape = xi1), 0)
  expect_length(qgev(c(0, 1), shape = xi1), 2)
  expect_length(qgev(seq_len(10) / 10, shape = xi1), 10)

  q2 <- -log(-log(pvec[2:4]))
  expect_equal(qgev(pvec, shape = xi2), c(-Inf, q2, Inf, NA))
  expect_length(qgev(seq_len(0), shape = xi2), 0)
  expect_length(qgev(c(0, 1), shape = xi2), 2)
  expect_length(qgev(seq_len(10) / 10, shape = xi2), 10)

  q3 <- ((-log(pvec[2:4])) ^ (-xi3) - 1) / xi3
  expect_equal(qgev(pvec, shape = xi3), c(-Inf, q3, up, NA))
  expect_length(qgev(seq_len(0), shape = xi3), 0)
  expect_length(qgev(c(0, 1), shape = xi3), 2)
  expect_length(qgev(seq_len(10) / 10, shape = xi3), 10)
})

test_that("pgev and qgev are consistent", {
  expect_equal(pgev(qgev(pvec, shape = xi1), shape = xi1), pvec)
  expect_equal(pgev(qgev(pvec, shape = xi2), shape = xi2), pvec)
  expect_equal(pgev(qgev(pvec, shape = xi3), shape = xi3), pvec)
})
