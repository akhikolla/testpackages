context("gscusum")

test_that("GSCUSUM results", {
  set.seed(201511)
  y <- rbinom(10,1,.2)
  blocks <- c(1,1,1,2,2,3,3,3,3,3)
  outcomes <- cbind(y, blocks)
  
  observed <- gscusum(input_outcomes = outcomes,
                      failure_probability = .2,
                      odds_multiplier = 2,
                      limit = 3,
                      quantiles = (0.5),
                      max_num_shuffles = 1000,
                      seed = 2910)
  exp_sig <- rep(0, 10)
  expect_equal(observed[,1], exp_sig)
  
  exp_mean <- c(0.0, 0.0, 0.0, 0.0, 0.0, 0.1, 0.2, 0.2, 0.2, 0.2)
  expect_equal(round(observed[,2],1), exp_mean)
  
  exp_med <- c(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1, 0.1, 0.1)
  expect_equal(round(observed[,3],1), exp_med)
})

