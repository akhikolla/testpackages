#Testing CondFunChisq.R
#Created by: Sajal Kumar
#Date : 10/16/2018

library(testthat)
library(FunChisq)

context("Testing cond.fun.chisq.test()")

test_that("Testing the conditional functional chi-square test", {

  # Test 1 : Perfect conditional functional dependency
  x1 = rep(c(0,0,1,1),400)
  x2 = rep(c(0,1,0,1),400)
  y = rep(c(0,1,1,0),400)

  data <- data.frame(x=x1,y=y,z=x2)

  test_res = cond.fun.chisq.test(x="x",y="y",z="z",data=data)

  expect_equivalent(test_res$p.value, 0)
  expect_equivalent(test_res$statistic, 1600)
  expect_equivalent(test_res$estimate, 1)


  # Test 2 : More natural example
  x1 = c(1,1,1,1,1,1,1,2,2,2,2,2,2,2,3,3,3,3,3,3,3)
  x2 = c(1,1,1,1,1,2,3,1,2,3,3,3,3,3,1,2,2,2,2,2,3)
  y = c(1,1,2,2,2,1,1,2,1,1,2,2,2,2,2,1,1,1,2,2,2)

  data <- data.frame(x=x1,y=y,z=x2)

  test_res_xy_z = cond.fun.chisq.test(x="x",y="y",z="z",data=data)
  test_res_xz_y = cond.fun.chisq.test(x="x",y="z",z="y",data=data)
  test_res_zy_x = cond.fun.chisq.test(x="z",y="y",z="x",data=data)

  expect_equivalent(signif(test_res_xy_z$p.value, 4), 0.6304)
  expect_equivalent(signif(test_res_xy_z$statistic, 2), 4.3)
  expect_equivalent(signif(test_res_xy_z$estimate, 4), 0.4971)

  expect_equivalent(signif(test_res_xz_y$p.value, 4), 0.01546)
  expect_equivalent(test_res_xz_y$statistic, 15.7)
  expect_equivalent(signif(test_res_xz_y$estimate, 4), 0.6386)

  expect_equivalent(signif(test_res_zy_x$p.value, 4), 0.3566)
  expect_equivalent(signif(test_res_zy_x$statistic, 3), 6.63)
  expect_equivalent(signif(test_res_zy_x$estimate, 4), 0.5778)


  # Test 3 : Conditional independency
  x1 = c(rep(1,10), rep(2,10), rep(3,10))
  x2 = c(rep(2,10), rep(1,10), rep(2,10))
  y = c(rep(1,10), rep(2,10), rep(1,10))

  test_res = cond.fun.chisq.test(x=x1,y=x2,z=y)

  expect_equivalent(test_res$p.value, 1)
  expect_equivalent(test_res$statistic, 0)
  expect_equivalent(test_res$estimate, 0)

})
