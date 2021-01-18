context("iwls")

# Check that iwls() works (on the newlyn data) over a range of thresholds

u <- stats::quantile(newlyn, probs = seq(0.1, 0.9, by = 0.1))

for (i in 1:length(u)) {
  res <- iwls(newlyn, u[i])
  test_that("iwls works", {
    testthat::expect_identical(res$conv, 0)
  })
}
