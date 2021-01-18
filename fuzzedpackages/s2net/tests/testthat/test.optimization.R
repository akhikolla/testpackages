context("Optimization")
library(s2net)

lm_test = function(train){
  obj = s2netR(train, s2Params(0))
  
  lm_fit = lm.fit(x = train$xL, y = train$yL)
  true_beta = unname(lm_fit$coefficients)
  found_beta = as.vector(obj$beta)
  expect_equal(found_beta, true_beta, tolerance = .01)
  
  lm_error = mean((train$xL%*%lm_fit$coefficients - train$yL)^2)
  error = mean((obj$intercept + train$xL%*%obj$beta - train$yL)^2)
  expect_equal(error, lm_error, tolerance = 0.00001)
}

test_that("s2net Fista solves the supervised problem correctly",{
  data("auto_mpg")
  train = s2Data(auto_mpg$P2$xL, auto_mpg$P2$yL, preprocess = TRUE)
  lm_test(train)
  
  train = s2Data(auto_mpg$P1$xL[,-1], auto_mpg$P1$yL, preprocess = TRUE)
  lm_test(train)
})
