context("choose_uk")

# Check that calling choose_uk() with vector arguments u and k gives
# the same results as calling kgaps repeatedly with scalar arguments

u <- stats::quantile(newlyn, probs = c(0.1, 0.90))
k_vals <- 1:3
cres <- choose_uk(newlyn, u = u, k = k_vals)
n_u <- length(u)
comp <- function(i, j) {
  return((i - 1) * n_u + j)
}
for (i in 1:length(k_vals)) {
  for (j in 1:length(u)) {
    res <- kgaps(newlyn, u[j], k_vals[i])
    temp <- cres$theta[[comp(i, j)]]
    # The calls will be different1
    temp$call <- res$call <- NULL
    test_that("choose_k agrees with kgaps", {
      testthat::expect_equal(temp, res)
    })
  }
}

# =============================== plot.choose_uk ==============================

# Check that plot.choose_uk works

# S&P 500 index
u <- quantile(sp500, probs = seq(0.1, 0.9, by = 0.1))
imt_theta <- choose_uk(sp500, u = u, k = 1:5)

ukplot <- plot(imt_theta)
test_that("plot.choose_uk works", {
  testthat::expect_identical(ukplot, NULL)
})
ukplot <- plot(imt_theta, y = "theta", ylim = c(0, 1), xlab = "my xlab", lwd = 2,
               col = 1:5)
test_that("plot.choose_b works, user plot args", {
  testthat::expect_identical(ukplot, NULL)
})

# One run parameter K, many thresholds u
u <- quantile(sp500, probs = seq(0.1, 0.9, by = 0.1))
imt_theta <- choose_uk(sp500, u = u, k = 1)
ukplot <- plot(imt_theta)
test_that("plot.choose_uk works", {
  testthat::expect_identical(ukplot, NULL)
})
ukplot <- plot(imt_theta, y = "theta", ylim = c(0, 1), xlab = "my xlab", lwd = 2,
               col = 1:5)
test_that("plot.choose_b works, user plot args", {
  testthat::expect_identical(ukplot, NULL)
})

# One threshold u, many run parameters K
u <- quantile(sp500, probs = 0.9)
imt_theta <- choose_uk(sp500, u = u, k = 1:5)
test_that("plot.choose_uk works", {
  testthat::expect_identical(ukplot, NULL)
})
ukplot <- plot(imt_theta, y = "theta", ylim = c(0, 1), xlab = "my xlab", lwd = 2,
               col = 1:5)
test_that("plot.choose_b works, user plot args", {
  testthat::expect_identical(ukplot, NULL)
})

