library(testthat)
library(dplyr)

Sys.setenv("R_TESTS" = "")

test_check("incgraph")

