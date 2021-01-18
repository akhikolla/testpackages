
library(testthat)
library(synthACS)

#----------------------------------------------------------
context("CRAN -- adding constraints")
#----------------------------------------------------------

test_that("constraint errors work", {
  # geography
  diamonds <- data.frame(
    carat= rexp(100),
    cut= factor(sample(c("A", "B", "C"), size= 100, replace= TRUE)),
    x= runif(100, min= 0, max= 10),
    y= runif(100, min= 0, max= 10),
    x= runif(100, min= 0, max= 10)
  )
  let <- letters
  
  # errors:
  expect_error(add_constraint(attr_name= "age", attr_totals= a, micro_data= diamonds))
})
