library(benchr)

b <- benchr::benchmark(1 + 1, 2 + 2)

expect_equal(class(b), c("benchmark", "data.frame"))
expect_equal(names(b), c("expr", "time"))
expect_equal(dim(b), c(200L, 2L))

expect_equal(class(b$expr), "factor")
expect_equal(class(b$time), "numeric")
expect_equal(levels(b$expr), c("1 + 1", "2 + 2"))

expect_true(!is.null(attr(b, "error")))
expect_equal(class(attr(b, "error")), "numeric")
expect_true(attr(b, "error") > 0)

expect_true(!is.null(attr(b, "precision")))
expect_equal(class(attr(b, "precision")), "numeric")
expect_true(attr(b, "precision") >= 0)
expect_true(attr(b, "precision") <= 1)

expect_true(!is.null(attr(b, "times")))
expect_equal(class(attr(b, "times")), "integer")
expect_equal(attr(b, "times"), 100L)

expect_true(!is.null(attr(b, "units")))
expect_identical(class(attr(b, "units")), "character")
expect_identical(attr(b, "units"), "s")

expect_true(!is.null(attr(b, "gc")))
expect_identical(attr(b, "gc"), FALSE)

b_gc <- benchr::benchmark(1 + 1, times = 1L, gcDuring = TRUE)
expect_true(!is.null(attr(b_gc, "gc")))
expect_identical(attr(b_gc, "gc"), TRUE)

b_r <- benchr::benchmark(1 + 1, times = 1L, order = "random")
expect_true(!is.null(attr(b_r, "order")))
expect_identical(attr(b_r, "order"), "random")

b_b <- benchr::benchmark(1 + 1, times = 1L, order = "block")
expect_true(!is.null(attr(b_b, "order")))
expect_identical(attr(b_b, "order"), "block")

b_i <- benchr::benchmark(1 + 1, times = 1L, order = "inorder")
expect_true(!is.null(attr(b_i, "order")))
expect_identical(attr(b_i, "order"), "inorder")

expect_error(benchmark(), "No expressions to benchmark")
