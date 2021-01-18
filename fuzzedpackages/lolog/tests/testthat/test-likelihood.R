context("Latent Order Likeihood")

#library(lolog)
library(testthat)
library(ergm)

test_that("lolog", {
  data(sampson)
  lol <- createLatentOrderLikelihood(samplike ~ edges())
  net <- lol$generateNetwork()$network
  
  lol <- createLatentOrderLikelihood(samplike ~ edges() | 1:net$size())
  net2 <- lol$generateNetwork()$network
  o2 <- net2[["__order__"]]
  expect_true(is.integer(o2))
  expect_equivalent(as.integer(o2), as.integer(0:17))
  
  lol <- createLatentOrderLikelihood(samplike ~ edges() | c(1, 1, 3:18))
  net3 <- lol$generateNetwork()$network
  o3 <- as.integer(net3[["__order__"]])
  op1 <- as.integer(0:17)
  op2 <- as.integer(c(1, 0, 2:17))
  expect_true(all(o3 == op1) | all(o3 == op2))
  
})
