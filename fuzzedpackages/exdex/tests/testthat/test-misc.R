context("log0const")

# Test that log0const(), a quick but less transparent version of
# log0const_slow() gives the correct answers

x <- c(runif(100), 0, 0, 0)
const <- -log(10)

test1 <- log0const_slow(x, const)
test2 <- log0const(x, const)

test_that("log0const slow = log0const fast", {
  testthat::expect_equal(test1, test2)
})

context("Empirical c.d.f.")

x <- newlyn
temp <- disjoint_maxima(newlyn)
y <- temp$y
x <- temp$x

test0 <- stats::ecdf(x)(y)
test1 <- ecdf1(x, y, lenx = length(x))
test2 <- ecdf2(x, y)
test3 <- ecdf3(x, y)

test_that("ecdf1 = stats::ecdf", {
  testthat::expect_equal(test0, test1)
})
test_that("ecdf2 = stats::ecdf", {
  testthat::expect_equal(test0, test2)
})
test_that("ecdf3 = stats::ecdf", {
  testthat::expect_equal(test0, test3)
})


