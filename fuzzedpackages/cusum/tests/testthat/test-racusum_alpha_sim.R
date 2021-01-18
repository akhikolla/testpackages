context("racusum_alpha_sim")

risks <- c(0.001, 0.01, 0.1, 0.002, 0.02, 0.2)

set.seed(2046)
patient_risks <- sample(x = risks, size = 100, replace = TRUE)


test_that("Output of alpha simulation", {
  expected_results <- 0.025
  works <- round(racusum_alpha_sim(patient_risks,
                                   odds_multiplier = 2,
                                   n_simulation = 1000,
                                   limit = 2.96,
                                   seed = 2046
                                   ),
                 3
                 )
  expect_equal(works, expected_results)
})

test_that("Warning if CL is not same direction as OM 1",
          expect_that(racusum_alpha_sim(patient_risks,
                                        odds_multiplier = .5,
                                        n_simulation = 1000,
                                        limit = 2.96,       
                                        seed = 2046
                                        ), 
                      gives_warning())
)

test_that("Warning if CL is not same direction as OM 2",
          expect_that(racusum_alpha_sim(patient_risks,
                                        odds_multiplier = 2,
                                        n_simulation = 1000,
                                        limit = -2.96,       
                                        seed = 2046
          ),
                      gives_warning())
)

test_that("Error if OM = 1",
          expect_that(racusum_alpha_sim(patient_risks,
                                        odds_multiplier =  1,
                                        n_simulation = 1000,
                                        limit = 2.96,       
                                        seed = 2046
          ), 
                      throws_error())
)

test_that("Output of alpha simulation improv", {
  expected_results <- 0.047
  works <- round(racusum_alpha_sim(patient_risks,
                                   odds_multiplier = .5,
                                   n_simulation = 1000,
                                   limit = -2.,
                                   seed = 2046
  ),
  3
  )
  expect_equal(works, expected_results)
})