context("cusum_limit_sim")

test_that("Output of Control limit simulation", {
  expected_results <- 2.85
  works <- round(
    cusum_limit_sim(
      failure_probability = 0.05,
      n_patients = 100,
      odds_multiplier = 2,
      n_simulation = 1000,
      alpha = 0.05,
      seed = 2046
    ),
    2
  )
  expect_equal(works, expected_results)
})

test_that("Output of Control limit simulation - Improvement", {
  expected_results <- -1.85
  works <- round(
    cusum_limit_sim(
      failure_probability = 0.05,
      n_patients = 100,
      odds_multiplier = 0.5,
      n_simulation = 1000,
      alpha = 0.05,
      seed = 2046
    ),
    2
  )
  expect_equal(works, expected_results)
})


test_that("Warning for recoding of failure probability",
          expect_that(cusum_limit_sim(
            failure_probability = 0.95,
            n_patients = 100,
            odds_multiplier = 0.5,
            n_simulation = 1000,
            alpha = 0.05,
            seed = 2046
          ), 
          gives_warning()
          )
)

test_that("Error if OM = 1",
          expect_that(cusum_limit_sim(
            failure_probability = 0.05,
            n_patients = 100,
            odds_multiplier = 1,
            n_simulation = 1000,
            alpha = 0.05,
            seed = 2046
          ), 
          throws_error()
          )
)
