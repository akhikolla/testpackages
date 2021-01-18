library(benchr)

p <- benchr:::timer_precision()
e <- benchr:::timer_error()

expect_equal(class(p), "numeric")
expect_equal(length(p), 1L)
expect_true(p < 0.001)

expect_equal(class(e), "numeric")
expect_equal(length(e), 1L)
expect_true(e > 0)

if (.Platform$OS.type != "windows") {
  expect_true(benchr:::do_timing(quote(Sys.sleep(0.1)), .GlobalEnv) >= 0.1)
  expect_true(benchr:::do_timing(quote(Sys.sleep(0.2)), .GlobalEnv) >= 0.2)
  expect_true(benchr:::do_timing(quote(Sys.sleep(0.3)), .GlobalEnv) >= 0.3)

  expect_true(benchr:::do_timing(quote(Sys.sleep(0.1)), .GlobalEnv) >= e)
}
