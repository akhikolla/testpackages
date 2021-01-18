library(testthat)
library(breathteststan)
Sys.unsetenv("R_TESTS") # https://github.com/r-lib/testthat/issues/603
test_check("breathteststan")
