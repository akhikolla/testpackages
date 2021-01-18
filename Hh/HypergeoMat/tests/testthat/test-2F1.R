context("2F1")

test_that("2F1 for scalar x", {
  obtained <- hypergeomPFQ(m = 25, a = c(1,2), b = c(3), x = 0.5)
  expected <- gsl::hyperg_2F1(1, 2, 3, 0.5)
  expect_equal(obtained, expected)
})

test_that("Some values for 2F1", {
  obtained <- hypergeomPFQ(m = 10, a = c(1,2), b = c(3), x = c(0.2,0.5))
  expect_equal(obtained, 1.79412894456143)
  obtained <- hypergeomPFQ(m = 10, a = c(1i,2), b = c(3i), x = c(0.2,0.5))
  expect_equal(obtained, 1.677558924-0.183004016i)
  obtained <- hypergeomPFQ(m = 10, a = c(1,2), b = c(3), x = c(0.2i,0.5))
  expect_equal(obtained, 1.513810425+0.20576184i)
  obtained <- hypergeomPFQ(m = 10, a = c(1,2i), b = c(3), x = c(0.2i,0.5))
  expect_equal(obtained, 0.7733140719+0.3092059749i)
})

test_that("Gauss formula", {
  a <- 1
  b <- 2
  c <- 9
  o1 <- mvgamma(c,3)*mvgamma(c-a-b,3)/mvgamma(c-a,3)/mvgamma(c-b,3)
  o2 <- hypergeomPFQ(100, c(a,b), c, c(1,1,1))
  expect_equal(o1, o2, tolerance = 1e-6)
  #
  a <- 1i
  o1 <- mvgamma(c,3)*mvgamma(c-a-b,3)/mvgamma(c-a,3)/mvgamma(c-b,3)
  o2 <- hypergeomPFQ(100, c(a,b), c, c(1,1,1))
  expect_equal(o1, o2, tolerance = 1e-6)
})

