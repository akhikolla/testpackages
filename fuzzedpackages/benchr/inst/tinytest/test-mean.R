library(benchr)

b <- benchr::benchmark(1 + 1, 2 + 2)
m <- mean(b)

expect_equal(class(m), c("summaryBenchmark", "data.frame"))
expect_equal(dim(m), c(2L, 7L))
expect_equal(names(m), c("expr", "n.eval", "mean", "trimmed", "lw.ci", "up.ci", "relative"))

expect_equal(class(m$expr), "factor")
expect_equal(levels(m$expr), c("1 + 1", "2 + 2"))
expect_true(all(sapply(m[-1], is.numeric)))

expect_equal(sum(!is.na(b$time)), sum(m$n.eval))

expect_equal(mean(b, relative = "mean"), mean(b, relative = 3))

expect_error(mean(b, relative = "z"))
expect_error(mean(b, relative = 100))

expect_true(is.null(mean(benchr::benchmark(1 + 1))$relative))
expect_true(is.null(mean(b, relative = NULL)$relative))

m1 <- mean(b, relative = "mean")
expect_equal(m1$relative, signif(m1$mean / min(m1$mean), 3))

m2 <- mean(b, relative = "trimmed", trim = 0.05)
expect_equal(m2$relative, signif(m2$trimmed / min(m2$trimmed), 3))

fun <- function(x, FUN, ...) sapply(split(x$time, x$expr), FUN, ...)
expect_equivalent(m$mean, fun(b, mean, na.rm = TRUE))
expect_equivalent(m$trimmed, fun(b, mean, na.rm = TRUE, trim = 0.05))

conf.ints <- unlist(fun(b, wilcox.test, conf.int = TRUE, conf.level = 0.95)[c(8, 17)])
expect_equivalent(m$lw.ci, conf.ints[c(1, 3)])
expect_equivalent(m$up.ci, conf.ints[c(2, 4)])

m3 <- mean(benchr::benchmark(1 + 1, times = 5))
expect_equal(class(m3$lw.ci), "numeric")
expect_equal(class(m3$up.ci), "numeric")
expect_true(m3$lw.ci > 0)
expect_true(m3$up.ci >= m3$lw.ci)

m4 <- mean(benchr::benchmark(1 + 1, times = 1))
expect_true(is.null(m4$lw.ci))
expect_true(is.null(m4$up.ci))
