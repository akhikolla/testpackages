## This file follows the structure of aaa.R in the free group package.

## Define some checker functions, and call them at the end.  They
## should all return TRUE if the package works, and stop with error if
## a test is failed.  Function checker1() has one argument, checker2()
## two, and checker3() has three.  

test_that("Test suite aab.R",{

  x <- clifford(list(1,2,1:4),1:3)
  expect_true(all(grades(grade(x,4)) == 4))

  expect_output(print(as.clifford(0)))
  expect_output(print(as.clifford(1)))
  expect_output(print(rcliff()))

  options("separate" = TRUE)
  expect_output(print(rcliff()))
  options("separate" = FALSE)


})
