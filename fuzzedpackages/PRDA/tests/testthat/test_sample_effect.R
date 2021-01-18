#################################
####    Test sample effect   ####
#################################

library(PRDA)

#----    input checks    ----

context("Evaluate functions to sample effect")

error_truncation <- "Argument 'tl' has to be greater than argument 'tu'"
message_truncation <- "Truncation could require long computational time"
error_tol <- "Truncation requires too long computational time, consider possible misspecification"



#----    sample_effect    ----
test_that("evaluate sample_effect", {
  expect_equal(with_seed(2020, sample_effect(FUN = function(x) rnorm(x),
                                             B_effect = 100)$effect_samples),
               with_seed(2020, rnorm(110)[11:110]))

  expect_error(sample_effect(FUN = function(x,y) rnorm(x,y), B_effect = 100))
  expect_error(sample_effect(FUN = function(x) x, B_effect = 100))
  expect_error(sample_effect(FUN = "ciao", B_effect = 100))
  expect_error(sample_effect(FUN = function(x) sample(c(1,"2"), x, replace = TRUE), B_effect = 100))

  expect_message(sample_effect(FUN = function(x) rnorm(x), B_effect = 100, tl=0), message_truncation)
  expect_error(sample_effect(FUN = function(x) rnorm(x), B_effect = 100, tl=0, tu=-1), error_truncation)

  expect_error(sample_effect(function(x)runif(x), B_effect = 100, tl=2), error_tol)
  })

#----
