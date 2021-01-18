#############################
####    Test utilities   ####
#############################

library(PRDA)


context("Evaluate utility functions")



#----    round_arg    ----

test_that("evaluate round_arg", {
  list_A <- list(a = 1:6/3,
                 b = "ciao")
  list_B <- list_A
  list_B$a <- round(list_A$a)

  expect_equal(round_arg(list_A, 0),
               list_B)
})

#----    list2data    ----

test_that("evaluate list2data", {
  list_A <- list(a = 7:12,
                 b = 1:6,
                 c = 5:10)
  df_A <- data.frame(a = 7:12,
                     b = 1:6)

  expect_equal(list2data(list_A, select = c("a","b")),
               df_A)
})


#----    sign_effect    ----

test_that("evaluate list2data", {
  expect_match(sign_effect(.3, "greater"),
                "0.3")
  expect_match(sign_effect(.3, "two_sided"),
               "± 0.3")
})

#----


