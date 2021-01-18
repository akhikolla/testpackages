context("testing psqn-bfgs")

test_that("we get the right result with the Rosenbrock Banana function", {
  fn <- function(x) {
    x1 <- x[1]
    x2 <- x[2]
    100 * (x2 - x1 * x1)^2 + (1 - x1)^2
  }
  gr_psqn <- function(x) {
    x1 <- x[1]
    x2 <- x[2]
    out <- c(-400 * x1 * (x2 - x1 * x1) - 2 * (1 - x1),
             200 *      (x2 - x1 * x1))
    attr(out, "value") <- 100 * (x2 - x1 * x1)^2 + (1 - x1)^2
    out
  }

  res <- psqn_bfgs(c(-1.2, 1), fn, gr_psqn)
  expect_equal(res$par, c(1, 1))
  expect_equal(res$value, 0)
  expect_true(res$convergence)
  expect_equal(res$info, 0L)
})
