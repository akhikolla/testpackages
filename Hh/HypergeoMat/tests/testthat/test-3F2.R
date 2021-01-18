context("3F2")

test_that("Kummer relation", {
  a <- c(1, 2, 3)
  b <- c(9, 10)
  c <- sum(b)-sum(a)
  p <- 4
  o1 <- mvgamma(b[2],p)*mvgamma(c,p)/mvgamma(b[2]-a[3],p)/mvgamma(c+a[3],p) *
    hypergeomPFQ(m=100, c(b[1]-a[1], b[1]-a[2], a[3]), c(b[1], c+a[3]), diag(p))
  o2 <- hypergeomPFQ(m=15, a, b, diag(p))
  expect_equal(o1, o2, tolerance = 1e-3)
  #
  a <- c(1, 2, 3i)
  b <- c(9i, 10)
  c <- sum(b)-sum(a)
  p <- 3
  o1 <- mvgamma(b[2],p)*mvgamma(c,p)/mvgamma(b[2]-a[3],p)/mvgamma(c+a[3],p) *
    hypergeomPFQ(m=100, c(b[1]-a[1], b[1]-a[2], a[3]), c(b[1], c+a[3]), diag(p))
  o2 <- hypergeomPFQ(m=15, a, b, diag(p))
  expect_equal(o1, o2, tolerance = 1e-5)
})
