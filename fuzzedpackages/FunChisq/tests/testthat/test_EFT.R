library(testthat)
library(FunChisq)

context("Testing the exact functional test")

test_that("Testing the EFT-DQP Dynamic and Quadratic Programming", {

  exact.functional.test <- EFTDQP

  x1 <- matrix(c(12, 26, 18, 0, 8, 12), nrow=2, byrow = TRUE)
  expect_equal(signif(exact.functional.test(x1), 8), 0.042556227)
  expect_equal(signif(exact.functional.test(t(x1)), 8), 0.027271581)

  x2 <- matrix(c(0,0,0,0,0,0,0,0,0), nrow=3, byrow = TRUE)
  expect_equal(exact.functional.test(x2), 1)

  x3 <- matrix(c(4,0,4,0,4,0,1,0,1), 3)
  expect_equal(signif(exact.functional.test(x3), 8), 0.002997003)
  expect_equal(signif(exact.functional.test(t(x3)), 8), 0.0065490065)

  x4 <- matrix(rep(10,25), nrow=5)
  expect_equal(exact.functional.test(x4), 1)

  x5 <- matrix(c(4,0,0,0,4,0,0,0,4), nrow=3, byrow = TRUE)
  expect_equal(signif(exact.functional.test(x5), 8), 0.00017316017)

  x6 <- matrix(c(2,0,0,2), nrow=2, byrow = TRUE)
  expect_equivalent(signif(exact.functional.test(x6), 8),
                    signif(stats::fisher.test(x6)$p.value, 8))


  x7 <- matrix(c(2,2,2,2), nrow=2, byrow = TRUE)
  expect_equivalent(signif(exact.functional.test(x7), 8),
                    signif(stats::fisher.test(x7)$p.value, 8))

  x8 <- matrix(c(0,10,15,20,5,0,25,0,0), nrow=3, byrow = TRUE)
  expect_equivalent(signif(exact.functional.test(x8), 8),
                    signif(fun.chisq.test(x8)$p.value, 8))

  x9 <- matrix(c(1,1,1,1,1,1,1,1,1), nrow=3, byrow = TRUE)
  expect_equal(exact.functional.test(x9), 1)

  x10 = matrix(c(
    6,    0,    8,    9,    6,
    0,   11,    0,    0,    0,
    0,    0,    0,    0,    0
  ),
  nrow=3, byrow=TRUE
  )
  expect_equal(exact.functional.test(x10), 4.325631e-10)

  x11 = matrix(c(
    0, 0, 0, 10, 9,
    0, 0, 0, 0, 0,
    8, 3, 10, 0, 0
  ),
  nrow=3, byrow=TRUE
  )
  expect_equal(exact.functional.test(x11), 1.523433e-11)

  x12 = matrix(c(
    8,    0,   11,    0,   11,
    0,    4,    0,    0,    0,
    0,    0,    0,    6,    0), nrow=3, byrow=TRUE
  )
  expect_equal(exact.functional.test(x12), 5.617703e-12)

  x13 = matrix( c(
    0, 2, 0,
    0, 2, 0,
    0, 2, 0 ), nrow = 3, byrow =T)
  expect_equal(exact.functional.test(x13), 1)
  expect_equal(exact.functional.test(t(x13)), 1)

  x14 = matrix( c(
    1, 2, 1,
    2, 2, 3,
    1, 1, 0
  ), nrow = 3, byrow =T)
  expect_equal(exact.functional.test(x14), 1)
  expect_equal(exact.functional.test(t(x14)), 0.8321678322)

  x15 = matrix( c(
    1, 0, 0,
    8, 2, 0,
    0, 1, 1
  ), nrow = 3, byrow =T)
  expect_equal(exact.functional.test(x15), 0.1818181818)
  expect_equal(exact.functional.test(t(x15)), 0.1818181818)

  x16 = matrix( c(0,    0,    0,    7,
                  7,    0,    0,    0,
                  0,    0,    0,    20,
                  0,    0,    0,    6
  ), nrow = 4, byrow =T)
  expect_equal(exact.functional.test(x16), 1.072756491e-07)
  expect_equal(exact.functional.test(t(x16)), 1.072756491e-07)

})

test_that("Testing the EFT-DP Dynamic Programming", {

  exact.functional.test <- EFTDP

  x1 <- matrix(c(12, 26, 18, 0, 8, 12), nrow=2, byrow = TRUE)
  expect_equal(signif(exact.functional.test(x1), 8), 0.042556227)
  expect_equal(signif(exact.functional.test(t(x1)), 8), 0.027271581)

  x2 <- matrix(c(0,0,0,0,0,0,0,0,0), nrow=3, byrow = TRUE)
  expect_equal(exact.functional.test(x2), 1)

  x3 <- matrix(c(4,0,4,0,4,0,1,0,1), 3)
  expect_equal(signif(exact.functional.test(x3), 8), 0.002997003)
  expect_equal(signif(exact.functional.test(t(x3)), 8), 0.0065490065)

  if (0) {
    x4 <- matrix(rep(10,25), nrow=5)
    expect_equal(exact.functional.test(x4), 1)
  }

  x5 <- matrix(c(4,0,0,0,4,0,0,0,4), nrow=3, byrow = TRUE)
  expect_equal(signif(exact.functional.test(x5), 8), 0.00017316017)

  x6 <- matrix(c(2,0,0,2), nrow=2, byrow = TRUE)
  expect_equivalent(signif(exact.functional.test(x6), 8),
                    signif(stats::fisher.test(x6)$p.value, 8))


  x7 <- matrix(c(2,2,2,2), nrow=2, byrow = TRUE)
  expect_equivalent(signif(exact.functional.test(x7), 8),
                    signif(stats::fisher.test(x7)$p.value, 8))

  x8 <- matrix(c(0,10,15,20,5,0,25,0,0), nrow=3, byrow = TRUE)
  expect_equivalent(signif(exact.functional.test(x8), 8),
                    signif(fun.chisq.test(x8)$p.value, 8))

  x9 <- matrix(c(1,1,1,1,1,1,1,1,1), nrow=3, byrow = TRUE)
  expect_equal(exact.functional.test(x9), 1)

  x10 = matrix(c(
    6,    0,    8,    9,    6,
    0,   11,    0,    0,    0,
    0,    0,    0,    0,    0
  ),
  nrow=3, byrow=TRUE
  )
  expect_equal(exact.functional.test(x10), 4.325631e-10)

  x11 = matrix(c(
    0, 0, 0, 10, 9,
    0, 0, 0, 0, 0,
    8, 3, 10, 0, 0
  ),
  nrow=3, byrow=TRUE
  )
  expect_equal(exact.functional.test(x11), 1.523433e-11)

  x12 = matrix(c(
    8,    0,   11,    0,   11,
    0,    4,    0,    0,    0,
    0,    0,    0,    6,    0), nrow=3, byrow=TRUE
  )
  expect_equal(exact.functional.test(x12), 5.617703e-12)

  x13 = matrix( c(
    0, 2, 0,
    0, 2, 0,
    0, 2, 0 ), nrow = 3, byrow =T)
  expect_equal(exact.functional.test(x13), 1)
  expect_equal(exact.functional.test(t(x13)), 1)

  x14 = matrix( c(
    1, 2, 1,
    2, 2, 3,
    1, 1, 0
  ), nrow = 3, byrow =T)
  expect_equal(exact.functional.test(x14), 1)
  expect_equal(exact.functional.test(t(x14)), 0.8321678322)

  x15 = matrix( c(
    1, 0, 0,
    8, 2, 0,
    0, 1, 1
  ), nrow = 3, byrow =T)
  expect_equal(exact.functional.test(x15), 0.1818181818)
  expect_equal(exact.functional.test(t(x15)), 0.1818181818)

  x16 = matrix( c(0,    0,    0,    7,
                  7,    0,    0,    0,
                  0,    0,    0,    20,
                  0,    0,    0,    6
  ), nrow = 4, byrow =T)
  expect_equal(exact.functional.test(x16), 1.072756491e-07)
  expect_equal(exact.functional.test(t(x16)), 1.072756491e-07)

})

test_that("Testing the EFT_QP Quadratic Programming", {

  exact.functional.test <- ExactFunctionalTest

  x1 <- matrix(c(12, 26, 18, 0, 8, 12), nrow=2, byrow = TRUE)
  expect_equal(signif(exact.functional.test(x1, TRUE), 8), 0.042556227)
  expect_equal(signif(exact.functional.test(t(x1), TRUE), 8), 0.027271581)
  expect_equal(signif(exact.functional.test(x1, TRUE), 8),
               signif(exact.functional.test(x1, FALSE), 8))

  x2 <- matrix(c(0,0,0,0,0,0,0,0,0), nrow=3, byrow = TRUE)
  expect_equal(exact.functional.test(x2, TRUE), 1)
  expect_equal(signif(exact.functional.test(x2, TRUE), 8),
               signif(exact.functional.test(x2, FALSE), 8))

  x3 <- matrix(c(4,0,4,0,4,0,1,0,1), 3)
  expect_equal(signif(exact.functional.test(x3, TRUE), 8), 0.002997003)
  expect_equal(signif(exact.functional.test(t(x3), TRUE), 8), 0.0065490065)
  expect_equal(signif(exact.functional.test(x3, TRUE), 8),
               signif(exact.functional.test(x3, FALSE), 8))

  if(0) { # This test case causes hang on windows. To be fixed.
    x4 <- matrix(rep(10,25), nrow=5)
    expect_equal(exact.functional.test(x4, TRUE), 1)
  }

  x5 <- matrix(c(4,0,0,0,4,0,0,0,4), nrow=3, byrow = TRUE)
  expect_equal(signif(exact.functional.test(x5, TRUE), 8), 0.00017316017)
  expect_equal(signif(exact.functional.test(x5, TRUE), 8),
               signif(exact.functional.test(x5, FALSE), 8))

  x6 <- matrix(c(2,0,0,2), nrow=2, byrow = TRUE)
  expect_equivalent(signif(exact.functional.test(x6, TRUE), 8),
                    signif(stats::fisher.test(x6)$p.value, 8))
  expect_equal(signif(exact.functional.test(x6, TRUE), 8),
               signif(exact.functional.test(x6, FALSE), 8))

  x7 <- matrix(c(2,2,2,2), nrow=2, byrow = TRUE)
  expect_equivalent(signif(exact.functional.test(x7, TRUE), 8),
                    signif(stats::fisher.test(x7)$p.value, 8))
  expect_equal(signif(exact.functional.test(x7, TRUE), 8),
               signif(exact.functional.test(x7, FALSE), 8))

  x8 <- matrix(c(0,10,15,20,5,0,25,0,0), nrow=3, byrow = TRUE)
  expect_equivalent(signif(exact.functional.test(x8, TRUE), 8),
                    signif(fun.chisq.test(x8)$p.value, 8))
  expect_equal(signif(exact.functional.test(x8, TRUE), 8),
               signif(exact.functional.test(x8, FALSE), 8))

  x9 <- matrix(c(1,1,1,1,1,1,1,1,1), nrow=3, byrow = TRUE)
  expect_equal(exact.functional.test(x9, TRUE), 1)
  expect_equal(signif(exact.functional.test(x9, TRUE), 8),
               signif(exact.functional.test(x9, FALSE), 8))
})
