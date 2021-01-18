context("test-check")

test_that("error if decline_time unit is not 'day'", {

   expect_error(decline_time(c(1:7300), unit = "days")   )

})

test_that("error if decline_time unit is not 'month'", {

   expect_error(decline_time(c(1:7300), unit = "months")   )

})

test_that("error if decline_time unit is not 'year'", {

   expect_error(decline_time(c(1:7300), unit = "years")   )

})

test_that("error if decline_time unit is not 'date'", {

   expect_error(decline_time(c(1:7300), unit = "Date")   )

})

test_that("error if exponential decline_param b value is not 'zero'", {

   expect_error(decline_param(input_unit = "Field", output_unit = "Field", fluid = "oil",
                                model = "exponential", qi = 1000, Di = 0.0015, b = 1, q_abnd = NULL))

})

test_that("error if harmonic decline_param b value is not 'one'", {

   expect_error(decline_param(input_unit = "Field", output_unit = "Field", fluid = "oil",
                              model = "harmonic", qi = 1000, Di = 0.0015, b = 0, q_abnd = NULL))

})

