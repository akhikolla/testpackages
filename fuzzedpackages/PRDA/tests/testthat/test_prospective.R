###############################
####    Test prospective   ####
###############################

library(PRDA)

#----    input checks    ----

context("prospective inputs specification")




test_that("inputs are correctly specified", {

  # Redefine function to avoid specify arguments each the times
  test_prospective <- function(effect_size = .3, power = .8, ratio_n = 1,
                         test_method = "pearson", alternative = "two_sided",
                         sig_level = .05, ratio_sd = 1, B = 10, tl = -Inf,
                         tu = Inf, B_effect = 10, sample_range = c(2, 1000),
                         eval_power = "mean", tol = .01,
                         display_message = FALSE, seed = 2020){
    with_seed(seed = seed,
              code = prospective(effect_size, power, ratio_n, test_method,
                                 alternative, sig_level, ratio_sd, B, tl, tu, B_effect,
                                 sample_range, eval_power, tol, display_message))
  }

  #----    Arguments    ----

  # effect_size
  effect_size_text <- "Argument 'effect_size' has to be a single numeric value or a function"
  expect_error(test_prospective(effect_size = Inf), effect_size_text)
  expect_error(test_prospective(effect_size = "ciao"), effect_size_text)
  expect_error(test_prospective(effect_size = c(1,2)), effect_size_text)

  # power
  power_text <- "Argument 'power' has to be a single value between 0 and 1"
  expect_error(test_prospective(power = Inf), power_text)
  expect_error(test_prospective(power = "ciao"), power_text)
  expect_error(test_prospective(power = c(.5,.6)), power_text)
  expect_error(test_prospective(power = 1), power_text)
  expect_error(test_prospective(power = -1), power_text)

  # ratio_n
  ratio_n_text <- "If specified, argument 'ratio_n' has to be a single positive value"
  expect_error(test_prospective(ratio_n = Inf), ratio_n_text)
  expect_error(test_prospective(ratio_n = "ciao"), ratio_n_text)
  expect_error(test_prospective(ratio_n = c(.5,.6)), ratio_n_text)
  expect_error(test_prospective(ratio_n = -.5), ratio_n_text)

  # sig_level
  sig_level_text <- "Argument 'sig_level' has to be a single value between 0 and 1"
  expect_error(test_prospective(sig_level = 1), sig_level_text)
  expect_error(test_prospective(sig_level = Inf), sig_level_text)
  expect_error(test_prospective(sig_level = 0), sig_level_text)
  expect_error(test_prospective(sig_level = "ciao"), sig_level_text)
  expect_error(test_prospective(sig_level = c(.1,.2)), sig_level_text)

  # ratio_sd
  ratio_sd_text <- "Argument 'ratio_sd' has to be a single finite number grater than 0"
  expect_error(test_prospective(ratio_sd = -1), ratio_sd_text)
  expect_error(test_prospective(ratio_sd = Inf), ratio_sd_text)
  expect_error(test_prospective(ratio_sd = 0), ratio_sd_text)
  expect_error(test_prospective(ratio_sd = "ciao"), ratio_sd_text)
  expect_error(test_prospective(ratio_sd = c(.1,.2)), ratio_sd_text)

  # B
  B_text <- "Argument 'B' has to be a single integer value grater than 1"
  expect_error(test_prospective(B = 1), B_text)
  expect_error(test_prospective(B = Inf), B_text)
  expect_error(test_prospective(B = "ciao"), B_text)
  expect_error(test_prospective(B = c(10,20)), B_text)

  # tl
  tl_text <- "Argument 'tl' has to be a single numeric value"
  expect_error(test_prospective(tl = "ciao"), tl_text)
  expect_error(test_prospective(tl = c(10,20)), tl_text)

  # tu
  tu_text <- "Argument 'tu' has to be a single numeric value"
  expect_error(test_prospective(tu = "ciao"), tu_text)
  expect_error(test_prospective(tu = c(10,20)), tu_text)

  # B_effect
  B_effect_text <- "Argument 'B_effect' has to be a single integer value grater than 1"
  expect_error(test_prospective(B_effect = 1), B_effect_text)
  expect_error(test_prospective(B_effect = Inf), B_effect_text)
  expect_error(test_prospective(B_effect = "ciao"), B_effect_text)
  expect_error(test_prospective(B_effect = c(10,20)), B_effect_text)

  # sample_range
  sample_range_text1 <- "Argument 'sample_range' has to be a length-2 numeric vector"
  sample_range_text2 <- "Argument 'sample_range' minimum has to be grater than 1 and less than sample range maximum"
  expect_error(test_prospective(sample_range = 1), sample_range_text1)
  expect_error(test_prospective(sample_range = c(10,20,20)), sample_range_text1)
  expect_error(test_prospective(sample_range = c(10,"20")), sample_range_text1)
  expect_error(test_prospective(sample_range = c(10,Inf)), sample_range_text1)
  expect_error(test_prospective(sample_range = c(10,8)), sample_range_text2)
  expect_error(test_prospective(sample_range = c(1,8)), sample_range_text2)

  # tol
  tol_text <- "Argument 'tol' has to be a single value between 0 and 1"
  expect_error(test_prospective(tol = 1), tol_text)
  expect_error(test_prospective(tol = Inf), tol_text)
  expect_error(test_prospective(tol = 0), tol_text)
  expect_error(test_prospective(tol = "ciao"), tol_text)
  expect_error(test_prospective(tol = c(.1,.2)), tol_text)

  # display_message
  display_text <- "Argument 'display_message' has to be logical"
  expect_error(test_prospective(display_message = 1), display_text)
  expect_error(test_prospective(display_message = Inf), display_text)
  expect_error(test_prospective(display_message = NULL), display_text)
  expect_error(test_prospective(display_message = "ciao"), display_text)
  expect_error(test_prospective(display_message = c(.1,.2)), display_text)


  #----    Other cases    ----

  # coherence effect_type and test_method
  # coherence_corr <- "If  'effect_type = correlation', argument 'test_method' has to be 'pearson'"
  # expect_error(test_prospective(effect_type = "correlation", test_method = "paired"),
  #              coherence_corr)
  # coherence_cohen <- "No appropriate 'test_method' for 'effect_type = cohen_d'"
  # expect_error(test_prospective(effect_type = "cohen_d", test_method = "pearson"),
  #              coherence_cohen)

  # correlation and ratio_n
  correlation_text <- "If 'test_method = pearson', argument 'ratio_n' is set to NULL"
  expect_warning(test_prospective(test_method = "pearson", ratio_n = 2),
                 correlation_text)

  # one_sample and ratio_n
  one_sample_text <- "If 'test_method = one_sample', argument 'ratio_n' must be set to NULL"
  expect_error(test_prospective(ratio_n = 3, test_method = "one_sample"),
               one_sample_text)

  # paired and ratio_n
  paired_text <- "If 'test_method = paired', argument 'ratio_n' has to be 1"
  expect_error(test_prospective(test_method = "paired", ratio_n = 2),
                 paired_text)
  expect_error(test_prospective(test_method = "paired", ratio_n = NULL),
               paired_text)

  # two_sample or welch and ratio_n
  t_test_text <- "Argument 'ratio_n' is required for the specified 'test_method'"
  expect_error(test_prospective(ratio_n = NULL, test_method = "two_sample"),
               t_test_text)
  expect_error(test_prospective(ratio_n = NULL, test_method = "welch"),
               t_test_text)

  # sample_range
  sample_range <- "Actual power = 0.16 with n = 100\n  try to increase maximum of sample_range > 100."
  expect_error(test_prospective(effect_size = .1, test_method = "two_sample",
                                sample_range =  c(5,100), B=100), sample_range)

  # tol issue
  tol_text <- "Required power according to tolerance value can not be obtained.\nIncrease tolerance value."
  expect_message(test_prospective(effect_size = .35, B=100, sample_range = c(40,120),
                                  tol=.001), tol_text)

  # welch and ratio_sd
  t_test_ratio_text1 <- "Argument 'ratio_sd' is required only for 'test_method = welch'"
  t_test_ratio_text2 <- "Argument 'ratio_sd' can not be 1 for 'test_method = welch'\n  Consider 'test_method = two_sample' instead"
  expect_error(test_prospective(ratio_sd = 1.5, test_method = "two_sample"),
               t_test_ratio_text1)
  expect_error(test_prospective(ratio_sd = 1, test_method = "welch"),
               t_test_ratio_text2)
})


#----    obtain same results    ----

test_that("same results as previous run", {
  expect_known_value(with_seed(seed = 2020,
                               prospective(effect_size = .3, power = .8, ratio_n = 1, B = 100, display_message = FALSE)$effect_info),
                     file = "test_cache/effect_info_single_pro_cor")
  expect_known_value(with_seed(seed = 2020,
                               prospective(effect_size = .3, power = .8, ratio_n = 1, B = 100, display_message = FALSE,
                                 test_method = "two_sample")$effect_info),
                     file = "test_cache/effect_info_single_pro_cohen")
  expect_known_value(with_seed(seed = 2020,
                               prospective(effect_size = .3, power = .8, ratio_n = 1, B = 100, display_message = FALSE)$prospective_res),
                     file="test_cache/res_corr_single_pro")
  expect_known_value(with_seed(seed = 2020,
                               prospective(effect_size = .3, power = .8, ratio_n = NULL, test_method = "one_sample", B = 100, display_message = FALSE)$prospective_res),
                     file="test_cache/res_one_sample_single_pro")
  expect_known_value(with_seed(seed = 2020,
                               prospective(effect_size = .3, power = .8, test_method = "paired", B = 100, display_message = FALSE)$prospective_res),
                     file="test_cache/res_cohen_paired_pro")
  expect_known_value(with_seed(seed = 2020,
                               prospective(effect_size = .3, power = .8, ratio_n = 1,  test_method = "welch", ratio_sd = 1.5, B = 100, display_message = FALSE)$prospective_res),
                     file="test_cache/res_cohen_single_pro")
  expect_known_value(with_seed(seed = 2020,
                               prospective(effect_size = .3, power = .8, ratio_n = 2,  test_method = "two_sample", B = 100, display_message = FALSE)$prospective_res),
                     file="test_cache/res_cohen_ratio_n_pro")


  expect_known_value(with_seed(seed = 2020,
                               prospective(effect_size = function(x) rnorm(x),  test_method = "welch", ratio_sd = 1.5, power = .8, ratio_n = 1, B = 100, B_effect = 10, display_message = FALSE)$effect_info),
                     file = "test_cache/effect_info_dist_pro")
  expect_known_value(with_seed(seed = 2020,
                               prospective(effect_size = function(x) rnorm(x), power = .8, ratio_n = 1,
                                 B = 100, B_effect = 10, eval_power = "mean", display_message = FALSE)$prospective_res),
                     file = "test_cache/res_corr_dist_pro")
  expect_known_value(with_seed(seed = 2020,
                               prospective(effect_size = function(x) rnorm(x), power = .8, ratio_n = NULL,  test_method = "one_sample",
                                 B = 100, B_effect = 10, eval_power = "mean", display_message = FALSE)$prospective_res),
                     file = "test_cache/res_one_sample_dist_pro")
  expect_known_value(with_seed(seed = 2020,
                               prospective(effect_size = function(x) rnorm(x), power = .8, ratio_n = 1,  test_method = "welch", ratio_sd = 1.5,
                                 B = 100, B_effect = 10, eval_power = "mean", display_message = FALSE)$prospective_res),
                     file = "test_cache/res_cohen_dist_pro")

})


#----
