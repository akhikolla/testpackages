context("racusum")

risks <- c(0.001, 0.01, 0.1, 0.002, 0.02, 0.2)

set.seed(2046)
patient_risks <- sample(x = risks, size = 10, replace = TRUE)

set.seed(2046)
outcomes <- as.logical(rbinom(10, 1, patient_risks))

test_that("CUSUM", {
  expected_results <- c(0, 0, 0, 0, 0, 0, 0.51, 1.02, 1.53, 1.53)
  works <-
    racusum(patient_risks,
            patient_outcomes = outcomes,
            limit = 2.96,
            odds_multiplier = 2,
            reset = TRUE
    )
  result <- round(works$ct, 2)
  expect_equal(result, expected_results)
})

test_that("Error if number of patients and patient risks don't match",
          expect_that(racusum(patient_risks,
                              patient_outcomes = rbinom(100,1,.1),
                              limit = 2.96,
                              odds_multiplier = 2,
                              reset = TRUE
          ), 
          throws_error()
          )
)

test_that("Warning if CL is not same direction as OM 1",
          expect_that(racusum(patient_risks,
                              patient_outcomes = outcomes,
                              limit = 2.96,
                              odds_multiplier = 0.5,
                              reset = TRUE
          ), 
          gives_warning()
          )
)

test_that("Warning if CL is not same direction as OM 2",
          expect_that(racusum(patient_risks,
                              patient_outcomes = outcomes,
                              limit = -2.96,
                              odds_multiplier = 2,
                              reset = TRUE
          ), 
          gives_warning()
          )
)


test_that("Warning if limit_method is passed",
          expect_that(racusum(patient_risks,
                              patient_outcomes = outcomes,
                              limit = 2.96,
                              odds_multiplier = 2,
                              reset = TRUE,
                              limit_method = "constant"
          ), 
          gives_warning()
          )
)