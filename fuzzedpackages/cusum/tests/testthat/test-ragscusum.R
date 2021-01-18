context("ragscusum")

test_that("RAGSCUSUM results", {
  set.seed(201511)
  y <- rbinom(10,1,.2)
  weight_f <- c(rep(0.5,5), rep(0.6,5))
  weight_s <- c(rep(-.2,5), rep(-.1,5))
  blocks <- c(1,1,1,2,2,3,3,3,3,3)
  outcomes <- cbind(y, weight_f, weight_s, blocks)
  
  observed <- ragscusum(input_ra_outcomes = outcomes,
                        limit = 3,
                        quantiles = (0.5),
                        max_num_shuffles = 10000,
                        seed = 2910)
  exp_sig <- rep(0, 10)
  expect_equal(observed[,1], exp_sig)
  
  exp_mean <- c(0.0, 0.0, 0.0, 0.0, 0.0, 0.1, 0.2, 0.3, 0.4, 0.4)
  expect_equal(round(observed[,2],1), exp_mean)
  
  exp_med <- c(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.4, 0.4 ,0.4)
  expect_equal(round(observed[,3],1), exp_med)
  
  
})


