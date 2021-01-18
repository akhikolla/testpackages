context("choose_b")

# =============================== plot.choose_b ===============================

# Check that plot.choose_b works

# S&P 500 index: vastly reduced number of block sizes, for speed
b_vals <- seq(from = 150, to = 350, by = 100)
res500 <- choose_b(sp500, b_vals)

bplot <- plot(res500)
test_that("plot.choose_b works", {
  testthat::expect_identical(bplot, NULL)
})

bplot <- plot(res500, ylim = c(0, 1), xlab = "my xlab", lwd = 2, col = "blue")
test_that("plot.choose_b works, user plot args", {
  testthat::expect_identical(bplot, NULL)
})
