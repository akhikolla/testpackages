library(testthat)
library(FunChisq)

if(0) {
  context("Testing comparative chi-squared tests")

  test_that("Testing the comparative functional chi-squared test", {
    x <- matrix(c(4,0,4,0,4,0,1,0,1), 3)
    y <- t(x)
    z <- matrix(c(1,0,1,4,0,4,0,4,0), 3)
    data <- list(x,y,z)
    expect_equivalent(signif(cp.fun.chisq.test(data)$p.value, 8), 0.00018762119)
    expect_equivalent(signif(cp.fun.chisq.test(data, method="nfchisq")$p.value, 8),
                      1.0052639e-07)
  })


  test_that("Testing the comparative chi-squared test", {

    x <- list()

    x[[1]] <- matrix(c(0,0,0,
                       0,0,0,
                       0,0,0), nrow=3)
    x[[2]] <- x[[1]]
    x[[3]] <- x[[1]]

    h <- cp.chisq.test(x)
    expect_equivalent(signif(h$p.value, 8), 1)
    expect_equivalent(signif(h$statistic, 8), 0)
    expect_equivalent(h$parameter, 0)

    h <- cp.chisq.test(x, method="nchisq")
    expect_equivalent(signif(h$p.value, 8), 1)

    x <- list()

    x[[1]] <- matrix(c(4,0,0,
                       0,4,0,
                       0,0,4), nrow=3)
    x[[2]] <- x[[1]]
    x[[3]] <- x[[1]]
    h <- cp.chisq.test(x)
    expect_equivalent(signif(h$p.value, 8), 1)
    expect_equivalent(signif(h$statistic, 8), 0)
    expect_equivalent(h$parameter, 8)

    h <- cp.chisq.test(x, method="nchisq")
    expect_equivalent(signif(h$p.value, 8), 0.97724987)

    x <- list()

    x[[1]] <- matrix(c(4,0,0,
                       0,4,0,
                       0,0,4), nrow=3)

    x[[2]] <- matrix(c(0,4,4,
                       4,0,4,
                       4,4,0), nrow=3)

    h <- cp.chisq.test(x)
    expect_equivalent(signif(h$p.value, 8), 2.8936962e-07)
    expect_equivalent(signif(h$statistic, 8), 36)
    expect_equivalent(h$parameter, 4)

    h <- cp.chisq.test(x, method="nchisq")
    expect_equivalent(signif(h$p.value, 8), 0)

    x <- matrix(c(4,0,4,0,4,0,1,0,1), 3)
    y <- t(x)
    z <- matrix(c(1,0,1,4,0,4,0,4,0), 3)
    data <- list(x,y,z)
    h <- cp.chisq.test(data)
    expect_equivalent(signif(h$p.value, 8), 1.3542453e-06)
    expect_equivalent(signif(h$statistic, 8), 42)
    expect_equivalent(h$parameter, 8)

    h <- cp.chisq.test(data, method="nchisq")
    expect_equivalent(signif(h$p.value, 8), 9.4795348e-18)
  })
}
