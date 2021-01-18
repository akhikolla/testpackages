context("cusum_limit_exact")

test_that("Exact calculation of CL",{
  x <- unname(round(cusum_limit_exact(failure_probability = 0.1,
                                      n_patients  = 10,
                                      odds_multiplier = 2,
                                      alpha = 0.05),5))
  expect_equal(x,
               1.41227)
}
  )

test_that("Error if OM = 1",
          expect_that(cusum_limit_exact(failure_probability = 0.1,
                                        n_patients  = 8,
                                        odds_multiplier = 1,
                                        alpha = 0.05), 
                      throws_error())
)

test_that("Message for large sample size",
          expect_message(cusum_limit_exact(failure_probability = 0.1,
                                        n_patients  = 13,
                                        odds_multiplier = 2,
                                        alpha = 0.05)
                      )
          )

test_that("Error for very large sample size",
          expect_that(cusum_limit_exact(failure_probability = 0.1,
                                        n_patients  = 150,
                                        odds_multiplier = 1,
                                        alpha = 0.05),
                      throws_error())
)