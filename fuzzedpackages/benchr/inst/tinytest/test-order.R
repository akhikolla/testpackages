library(benchr)

set.seed(0)
l <- list(a = 1, b = 2, c = 3, d = 4)
b1 <- benchr:::make_order(l, 1, "block")
i1 <- benchr:::make_order(l, 1, "inorder")
r1 <- benchr:::make_order(l, 1, "random")
b5 <- benchr:::make_order(l, 5, "block")
i5 <- benchr:::make_order(l, 5, "inorder")
r5 <- benchr:::make_order(l, 5, "random")

expect_equal(class(b1), "factor")
expect_equal(levels(b1), names(l))

expect_equal(length(b1), length(l))
expect_equal(length(i1), length(l))
expect_equal(length(r1), length(l))

expect_equal(length(b5), length(l) * 5)
expect_equal(length(i5), length(l) * 5)
expect_equal(length(r5), length(l) * 5)

expect_equal(b5, as.factor(rep(names(l), each = 5)))
expect_equal(i5, as.factor(rep(names(l), times = 5)))

expect_false(identical(benchr:::make_order(l, 5, "random"), benchr:::make_order(l, 5, "random")))
