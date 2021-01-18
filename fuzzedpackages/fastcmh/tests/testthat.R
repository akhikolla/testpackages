library(testthat)
library(fastcmh)

Sys.setenv("R_TESTS" = "")
test_check("fastcmh")
