library(benchr)

expect_equal(class(benchr:::dots(1, 2, 3)), "list")
expect_equal(benchr:::dots(1, 2, {3}), list(1, 2, quote({3})))
expect_equal(benchr:::dots("a", "b", {"c"}), list("a", "b", quote({"c"})))
expect_equal(benchr:::dots(1:2, 1:3), list(quote(1:2), quote(1:3)))
expect_equal(benchr:::dots(a, b), list(as.symbol("a"), as.symbol("b")))
expect_equal(
  benchr:::dots(a = mean(1), b = mean(2)),
  list(a = quote(mean(1)), b = quote(mean(2)))
)

expect_equal(class(benchr:::named_dots(NULL)), "list")
expect_equal(names(benchr:::named_dots(1, 2, a, b = 1:5, {3})), c("1", "2", "a", "b", "{ 3 }"))
expect_equal(names(benchr:::named_dots(NULL)), c("NULL"))
expect_equal(names(benchr:::named_dots({{NULL}})), c("{ { NULL } }"))
expect_equal(names(benchr:::named_dots({{NULL; NULL}})), c("{ { NULL; NULL } }"))
expect_equal(names(benchr:::named_dots({f <- function() {NULL; NULL}})), c("{ f <- function() { NULL; NULL } }"))
