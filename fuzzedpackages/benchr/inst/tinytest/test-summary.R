library(benchr)

b <- benchr::benchmark(1 + 1, 2 + 2)
s <- summary(b)

expect_equal(class(s), c("summaryBenchmark", "data.frame"))
expect_equal(dim(s), c(2L, 10L))
expect_equal(names(s), c("expr", "n.eval", "min", "lw.qu", "median", "mean", "up.qu", "max", "total", "relative"))

expect_equal(class(s$expr), "factor")
expect_true(all(sapply(s[-1], is.numeric)))
expect_equal(levels(s$expr), c("1 + 1", "2 + 2"))

fun <- function(x, FUN, ...) sapply(split(x$time, x$expr), FUN, ...)
expect_equivalent(s$min, fun(b, min, na.rm = TRUE))
expect_equivalent(s$max, fun(b, max, na.rm = TRUE))
expect_equivalent(s$mean, fun(b, mean, na.rm = TRUE))
expect_equivalent(s$median, fun(b, median, na.rm = TRUE))
expect_equivalent(s$total, fun(b, sum, na.rm = TRUE))
expect_equivalent(s$lw.qu, fun(b, quantile, na.rm = TRUE, prob = 0.25))
expect_equivalent(s$up.qu, fun(b, quantile, na.rm = TRUE, prob = 0.75))
expect_equal(sum(!is.na(b$time)), sum(s$n.eval))

expect_error(summary(b, relative = c("min", "max")))
expect_error(summary(b, relative = "z"))
expect_error(summary(b, relative = 100))

expect_true(is.null(mean(benchr::benchmark(1 + 1))$relative))
expect_true(is.null(mean(b, relative = NULL)$relative))

s1 <- summary(b, relative = "median")
expect_equal(s1$relative, signif(s1$median / min(s1$median), 3))

s2 <- summary(b, relative = "mean")
expect_equal(s2$relative, signif(s2$mean / min(s2$mean), 3))

expect_equal(summary(b, relative = "min"), summary(b, relative = 3))
