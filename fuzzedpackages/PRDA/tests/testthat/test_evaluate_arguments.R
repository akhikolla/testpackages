######################################
####    Test evaluate arguments   ####
######################################

library(PRDA)

#----    input checks    ----

context("Evaluate arguments")

groups <- with_seed(2020, list(x = rnorm(15, .3),
                            y = rnorm(15)))
diff <- groups$x- groups$y
mx <- mean(groups$x)
my <- mean(groups$y)
nx <- length(groups$x)
ny <- length(groups$y)
ny2 <- ny+2
sig_level <- .10
mu <- .5
mu2 <-  -.5




#----    eval_test_method    ----

test_that("evaluate the correct test method", {

  # Redefine function to avoid specify arguments each the times
  test_eval_test_method <- function(effect_type = "cohen_d", effect_target = .3,
                                    test_method, sample_n1 = nx, sample_n2 = ny,
                                    sig_level = .05, alternative = "two_sided",
                                    ratio_sd = 1){
    with_seed(2020, eval_test_method(
      effect_type, effect_target, test_method, sample_n1, sample_n2,
      alternative, sig_level, ratio_sd))
  }

  # Cohen's d
  expect_equal(test_eval_test_method(test_method = "one_sample", sample_n2 = NULL),
               t.test(groups$x))
  expect_equal(test_eval_test_method(test_method = "two_sample", sample_n2 = ny),
               t.test(groups$x, groups$y, var.equal = TRUE))
  expect_equal(test_eval_test_method(test_method = "welch", sample_n2 = ny),
               t.test(groups$x, groups$y))

  # paired
  groups <- with_seed(2020, list(x = rnorm(15, .3*sqrt(2), 1),
                                y = rnorm(15, 0, 1)))
  expect_equal(test_eval_test_method(test_method = "paired", sample_n2 = ny),
               t.test(groups$x, groups$y, paired = TRUE))

  # welch and ratio_n2
  groups <- with_seed(2020, list(x = rnorm(15, .3*sqrt(5/2), 2),
                                 y = rnorm(15, 0, 1)))
  expect_equal(test_eval_test_method(
    test_method = "welch", sample_n2 = ny, ratio_sd = 2),
               t.test(groups$x, groups$y))

  #Correlation
  groups <- with_seed(2020, sample_obs_cor(nx, .3))
  expect_equal(test_eval_test_method(
    effect_type = "correlation", test_method = "pearson", sample_n2 = ny),
               cor.test(groups$x, groups$y))
})



#----    eval_effect_size    ----

test_that("evaluate the correct effect size", {
  expect_equal(eval_effect_size(effect_type = "correlation", effect_size = .3)$effect_samples, .3)
  expect_equal(eval_effect_size(effect_type = "correlation", effect_size = .3)$effect_summary[["Mean"]], .3)
  expect_equal(with_seed(2020,eval_effect_size(effect_type = "cohen_d", effect_size = function(x) rnorm(x),B_effect = 100)$effect_samples),
               with_seed(2020, rnorm(110)[11:110]))

  error_effect_corr <- "If 'effect_type = correlation', argument 'effect_size' must be between -1 and 1"
  expect_error(eval_effect_size(effect_type = "correlation", effect_size = 2),
               error_effect_corr)
  message_trunc_cor <- "If 'effect_type = correlation', effect_size distribution is truncated between -1 and 1\n"
  message_trunc <- "Truncation could require long computational time\n"
  res <- evaluate_promise(eval_effect_size(effect_type = "correlation", effect_size = function(x) rnorm(x),B_effect = 100))
  expect_equal(res$messages[1], message_trunc_cor)
  expect_equal(res$messages[2], message_trunc)
})


#----    eval_rgn_function    ----

test_that("evaluate the correct number of sampled values", {
  expect_true(eval_rgn_function(FUN = function(x) rnorm(x), n = 10))
  expect_false(eval_rgn_function(FUN = function(x) x, n = 10))
  expect_false(eval_rgn_function(FUN = function(x) sample(c(1,"2"),
                                       x, replace = TRUE), n = 10))

})


#----    eval_effect_type    ----

test_that("evaluate the correct effect type", {
  expect_match(eval_effect_type(test_method = "pearson"), "correlation")
  expect_match(eval_effect_type(test_method = "two_sample"), "cohen_d")
  expect_match(eval_effect_type(test_method = "welch"), "cohen_d")
  expect_match(eval_effect_type(test_method = "paired"), "cohen_d")
  expect_match(eval_effect_type(test_method = "one_sample"), "cohen_d")
})

#----
