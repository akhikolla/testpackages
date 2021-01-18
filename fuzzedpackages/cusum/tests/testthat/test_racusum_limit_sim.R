context("racusum_limit_sim")

test_that("Output of Control limit simulation", {
  
  expected_results <- 2.581
  observed_results <- round(racusum_limit_sim(patient_risks = c(rep(0.1,10), rep(0.15, 10), rep(0.2,10)),
                    odds_multiplier = 2,
                    n_simulation = 1000,
                    alpha = 0.05,
                    seed = 2910), 3)
  
  expect_equal(observed_results, expected_results)
})


test_that("Error if OM = 1",
          expect_that(racusum_limit_sim(patient_risks = c(rep(0.1,10), rep(0.15, 10), rep(0.2,10)),
                                        odds_multiplier = 1,
                                        n_simulation = 1000,
                                        alpha = 0.05,
                                        seed = 2910), 
                      throws_error())
)

test_that("Output of Control limit simulation improv", {
  
  expected_results <- -1.747
  observed_results <- round(racusum_limit_sim(patient_risks = c(rep(0.1,10), rep(0.15, 10), rep(0.2,10)),
                                              odds_multiplier = 0.5,
                                              n_simulation = 1000,
                                              alpha = 0.05,
                                              seed = 2910), 3)
  
  expect_equal(observed_results, expected_results)
})
