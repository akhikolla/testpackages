context("Jack expansions")

test_that("Schur expansion", {
  library(jack)
  genpoch <- function(a, kappa, alpha){
    prod(sapply(seq_along(kappa), function(i){
      prod(a - (i-1L)/alpha + seq_len(kappa[i]) - 1)
    }))
  }
  m <- 2
  alpha <- 1
  a <- c(2, 3)
  b <- c(4, 1i)
  x <- c(0.3i, 0.7)
  o1 <-
    genpoch(a[1], c(0), alpha) * genpoch(a[2], c(0), alpha) /
    genpoch(b[1], c(0), alpha) / genpoch(b[2], c(0), alpha) *
    Schur(x, c(0)) +
    genpoch(a[1], c(1), alpha) * genpoch(a[2], c(1), alpha) /
    genpoch(b[1], c(1), alpha) / genpoch(b[2], c(1), alpha) *
    Schur(x, c(1)) +
    genpoch(a[1], c(1,1), alpha) * genpoch(a[2], c(1,1), alpha) /
    genpoch(b[1], c(1,1), alpha) / genpoch(b[2], c(1,1), alpha) *
    Schur(x, c(1,1)) / 2 +
    genpoch(a[1], c(2), alpha) * genpoch(a[2], c(2), alpha) /
    genpoch(b[1], c(2), alpha) / genpoch(b[2], c(2), alpha) *
    Schur(x, c(2)) / 2
  o2 <- hypergeomPFQ(m, a, b, x, alpha)
  expect_equal(o1, o2)
  #
  a <- c(2, 3)
  b <- c(4, 5)
  x <- c(0.3, 0.7)
  o1 <-
    genpoch(a[1], c(0), alpha) * genpoch(a[2], c(0), alpha) /
    genpoch(b[1], c(0), alpha) / genpoch(b[2], c(0), alpha) *
    Schur(x, c(0)) +
    genpoch(a[1], c(1), alpha) * genpoch(a[2], c(1), alpha) /
    genpoch(b[1], c(1), alpha) / genpoch(b[2], c(1), alpha) *
    Schur(x, c(1)) +
    genpoch(a[1], c(1,1), alpha) * genpoch(a[2], c(1,1), alpha) /
    genpoch(b[1], c(1,1), alpha) / genpoch(b[2], c(1,1), alpha) *
    Schur(x, c(1,1)) / 2 +
    genpoch(a[1], c(2), alpha) * genpoch(a[2], c(2), alpha) /
    genpoch(b[1], c(2), alpha) / genpoch(b[2], c(2), alpha) *
    Schur(x, c(2)) / 2
  o2 <- hypergeomPFQ(m, a, b, x, alpha)
  expect_equal(o1, o2)
})

test_that("Zonal expansion", {
  library(jack)
  genpoch <- function(a, kappa, alpha){
    prod(sapply(seq_along(kappa), function(i){
      prod(a - (i-1L)/alpha + seq_len(kappa[i]) - 1)
    }))
  }
  m <- 2
  alpha <- 2
  a <- c(2, 3)
  b <- c(4, 1i)
  x <- c(0.3i, 0.7)
  o1 <-
    genpoch(a[1], c(0), alpha) * genpoch(a[2], c(0), alpha) /
    genpoch(b[1], c(0), alpha) / genpoch(b[2], c(0), alpha) *
    Zonal(x, c(0)) +
    genpoch(a[1], c(1), alpha) * genpoch(a[2], c(1), alpha) /
    genpoch(b[1], c(1), alpha) / genpoch(b[2], c(1), alpha) *
    Zonal(x, c(1)) +
    genpoch(a[1], c(1,1), alpha) * genpoch(a[2], c(1,1), alpha) /
    genpoch(b[1], c(1,1), alpha) / genpoch(b[2], c(1,1), alpha) *
    Zonal(x, c(1,1)) / 2 +
    genpoch(a[1], c(2), alpha) * genpoch(a[2], c(2), alpha) /
    genpoch(b[1], c(2), alpha) / genpoch(b[2], c(2), alpha) *
    Zonal(x, c(2)) / 2
  o2 <- hypergeomPFQ(m, a, b, x, alpha)
  expect_equal(o1, o2)
  #
  a <- c(2, 3)
  b <- c(4, 5)
  x <- c(0.3, 0.7)
  o1 <-
    genpoch(a[1], c(0), alpha) * genpoch(a[2], c(0), alpha) /
    genpoch(b[1], c(0), alpha) / genpoch(b[2], c(0), alpha) *
    Zonal(x, c(0)) +
    genpoch(a[1], c(1), alpha) * genpoch(a[2], c(1), alpha) /
    genpoch(b[1], c(1), alpha) / genpoch(b[2], c(1), alpha) *
    Zonal(x, c(1)) +
    genpoch(a[1], c(1,1), alpha) * genpoch(a[2], c(1,1), alpha) /
    genpoch(b[1], c(1,1), alpha) / genpoch(b[2], c(1,1), alpha) *
    Zonal(x, c(1,1)) / 2 +
    genpoch(a[1], c(2), alpha) * genpoch(a[2], c(2), alpha) /
    genpoch(b[1], c(2), alpha) / genpoch(b[2], c(2), alpha) *
    Zonal(x, c(2)) / 2
  o2 <- hypergeomPFQ(m, a, b, x, alpha)
  expect_equal(o1, o2)
})

test_that("ZonalQ expansion", {
  library(jack)
  genpoch <- function(a, kappa, alpha){
    prod(sapply(seq_along(kappa), function(i){
      prod(a - (i-1L)/alpha + seq_len(kappa[i]) - 1)
    }))
  }
  m <- 2
  alpha <- 1/2
  a <- c(2, 3)
  b <- c(4, 1i)
  x <- c(0.3i, 0.7)
  o1 <-
    genpoch(a[1], c(0), alpha) * genpoch(a[2], c(0), alpha) /
    genpoch(b[1], c(0), alpha) / genpoch(b[2], c(0), alpha) *
    ZonalQ(x, c(0)) +
    genpoch(a[1], c(1), alpha) * genpoch(a[2], c(1), alpha) /
    genpoch(b[1], c(1), alpha) / genpoch(b[2], c(1), alpha) *
    ZonalQ(x, c(1)) +
    genpoch(a[1], c(1,1), alpha) * genpoch(a[2], c(1,1), alpha) /
    genpoch(b[1], c(1,1), alpha) / genpoch(b[2], c(1,1), alpha) *
    ZonalQ(x, c(1,1)) / 2 +
    genpoch(a[1], c(2), alpha) * genpoch(a[2], c(2), alpha) /
    genpoch(b[1], c(2), alpha) / genpoch(b[2], c(2), alpha) *
    ZonalQ(x, c(2)) / 2
  o2 <- hypergeomPFQ(m, a, b, x, alpha)
  expect_equal(o1, o2)
  #
  a <- c(2, 3)
  b <- c(4, 5)
  x <- c(0.3, 0.7)
  o1 <-
    genpoch(a[1], c(0), alpha) * genpoch(a[2], c(0), alpha) /
    genpoch(b[1], c(0), alpha) / genpoch(b[2], c(0), alpha) *
    ZonalQ(x, c(0)) +
    genpoch(a[1], c(1), alpha) * genpoch(a[2], c(1), alpha) /
    genpoch(b[1], c(1), alpha) / genpoch(b[2], c(1), alpha) *
    ZonalQ(x, c(1)) +
    genpoch(a[1], c(1,1), alpha) * genpoch(a[2], c(1,1), alpha) /
    genpoch(b[1], c(1,1), alpha) / genpoch(b[2], c(1,1), alpha) *
    ZonalQ(x, c(1,1)) / 2 +
    genpoch(a[1], c(2), alpha) * genpoch(a[2], c(2), alpha) /
    genpoch(b[1], c(2), alpha) / genpoch(b[2], c(2), alpha) *
    ZonalQ(x, c(2)) / 2
  o2 <- hypergeomPFQ(m, a, b, x, alpha)
  expect_equal(o1, o2)
})

test_that("MSF expansion", {
  library(jack)
  genpoch <- function(a, kappa, alpha){
    prod(sapply(seq_along(kappa), function(i){
      prod(a - (i-1L)/alpha + seq_len(kappa[i]) - 1)
    }))
  }
  m <- 2
  alpha <- Inf
  a <- c(2, 3)
  b <- c(4, 1i)
  x <- c(0.3i, 0.7)
  o1 <-
    genpoch(a[1], c(0), alpha) * genpoch(a[2], c(0), alpha) /
    genpoch(b[1], c(0), alpha) / genpoch(b[2], c(0), alpha) *
    MSF(x, c(0)) +
    genpoch(a[1], c(1), alpha) * genpoch(a[2], c(1), alpha) /
    genpoch(b[1], c(1), alpha) / genpoch(b[2], c(1), alpha) *
    MSF(x, c(1)) +
    genpoch(a[1], c(1,1), alpha) * genpoch(a[2], c(1,1), alpha) /
    genpoch(b[1], c(1,1), alpha) / genpoch(b[2], c(1,1), alpha) *
    MSF(x, c(1,1)) / 2 * 2 +
    genpoch(a[1], c(2), alpha) * genpoch(a[2], c(2), alpha) /
    genpoch(b[1], c(2), alpha) / genpoch(b[2], c(2), alpha) *
    MSF(x, c(2)) / 2
  o2 <- hypergeomPFQ(m, a, b, x, alpha = 10000000)
  expect_equal(o1, o2, tolerance = 1e-7)
  #
  a <- c(2, 3)
  b <- c(4, 1)
  x <- c(0.3, 0.7)
  o1 <-
    genpoch(a[1], c(0), alpha) * genpoch(a[2], c(0), alpha) /
    genpoch(b[1], c(0), alpha) / genpoch(b[2], c(0), alpha) *
    MSF(x, c(0)) +
    genpoch(a[1], c(1), alpha) * genpoch(a[2], c(1), alpha) /
    genpoch(b[1], c(1), alpha) / genpoch(b[2], c(1), alpha) *
    MSF(x, c(1)) +
    genpoch(a[1], c(1,1), alpha) * genpoch(a[2], c(1,1), alpha) /
    genpoch(b[1], c(1,1), alpha) / genpoch(b[2], c(1,1), alpha) *
    MSF(x, c(1,1)) / 2 * 2 +
    genpoch(a[1], c(2), alpha) * genpoch(a[2], c(2), alpha) /
    genpoch(b[1], c(2), alpha) / genpoch(b[2], c(2), alpha) *
    MSF(x, c(2)) / 2
  o2 <- hypergeomPFQ(m, a, b, x, alpha = 10000000)
  expect_equal(o1, o2, tolerance = 1e-7)
})
