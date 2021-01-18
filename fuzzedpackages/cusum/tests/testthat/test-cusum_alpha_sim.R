context("cusum_alpha_sim")


test_that("Output of alpha simulation", {
  expected_results <- 0.038
  works <- round(
    cusum_alpha_sim(
      failure_probability = 0.05,
      n_patients = 100,
      odds_multiplier = 2,
      n_simulation = 10000,
      limit = 2.96,
      seed = 2046
    ),
    3
  )
  expect_equal(works, expected_results)
})

test_that("Output of alpha simulation improv", {
  expected_results <- 0.026
  works <- round(
    cusum_alpha_sim(
      failure_probability = 0.1,
      n_patients = 100,
      odds_multiplier = 0.5,
      n_simulation = 10000,
      limit = - 2.96,
      seed = 2046
    ),
    3
  )
  expect_equal(works, expected_results)
})
test_that("Warning if CL is not same direction as OM 1",
          expect_that(cusum_alpha_sim(
            failure_probability = 0.1,
            n_patients = 100,
            odds_multiplier = 0.5,
            n_simulation = 1000,
            limit =  2.96,
            seed = 2046
          ), 
                      gives_warning())
)

test_that("Warning if CL is not same direction as OM 2",
          expect_that(cusum_alpha_sim(
            failure_probability = 0.1,
            n_patients = 100,
            odds_multiplier = 2,
            n_simulation = 1000,
            limit = - 2.96,
            seed = 2046
          ), 
          gives_warning())
)
test_that("Error if OM = 1",
          expect_that(cusum_alpha_sim(
            failure_probability = 0.1,
            n_patients = 100,
            odds_multiplier = 1,
            n_simulation = 1000,
            limit = 2,
            seed = 2046), 
            throws_error())
)

test_that("Warning for recoding failure_prob",
          expect_that(cusum_alpha_sim(
            failure_probability = .95,
            n_patients = 100,
            odds_multiplier = 2,          
            n_simulation = 1000,
            limit = 2.96,            
            seed = 2046
            ),
            gives_warning()
          )
)
