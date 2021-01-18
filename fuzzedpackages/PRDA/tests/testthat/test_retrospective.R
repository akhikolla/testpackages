#################################
####    Test retrospective   ####
#################################

library(PRDA)

#----    input checks    ----

context("retrospective inputs specification")

test_that("inputs are correctly specified", {

  # Redefine function to avoid specify arguments each the times
  test_retrospective <- function(effect_size = .3, sample_n1 = 20,
                          sample_n2 = NULL, test_method = "pearson",
                          alternative = "two_sided", sig_level = .05,
                          ratio_sd = 1, B = 10, tl = -Inf, tu = Inf,
                          B_effect = 10, display_message = TRUE, seed = 2020){
    with_seed(seed = seed,
              retrospective(effect_size, sample_n1, sample_n2, test_method,
                  alternative, sig_level, ratio_sd, B, tl, tu, B_effect,
                  display_message))
  }

  #----    Arguments    ----

  # sample_n1
  sample_n1_text <- "Argument 'sample_n1' has to be a single integer value grater than 1"
  expect_error(test_retrospective(sample_n1 = 1), sample_n1_text)
  expect_error(test_retrospective(sample_n1 = Inf), sample_n1_text)
  expect_error(test_retrospective(sample_n1 = -3), sample_n1_text)
  expect_error(test_retrospective(sample_n1 = "ciao"), sample_n1_text)
  expect_error(test_retrospective(sample_n1 = c(10,20)), sample_n1_text)

  # effect_size
  effect_size_text <- "Argument 'effect_size' has to be a single numeric value or a function"
  expect_error(test_retrospective(effect_size = Inf), effect_size_text)
  expect_error(test_retrospective(effect_size = "ciao"), effect_size_text)
  expect_error(test_retrospective(effect_size = c(.5,.7)), effect_size_text)

  # sample_n2
  sample_n2_text <- "If specified, argument 'sample_n2' has to be a single integer value grater than 1"
  expect_error(test_retrospective(sample_n2 = 1), sample_n2_text)
  expect_error(test_retrospective(sample_n2 = Inf), sample_n2_text)
  expect_error(test_retrospective(sample_n2 = -3), sample_n2_text)
  expect_error(test_retrospective(sample_n2 = "ciao"), sample_n2_text)
  expect_error(test_retrospective(sample_n2 = c(10,20)), sample_n2_text)

  # sig_level
  sig_level_text <- "Argument 'sig_level' has to be a single value between 0 and 1"
  expect_error(test_retrospective(sig_level = 1), sig_level_text)
  expect_error(test_retrospective(sig_level = Inf), sig_level_text)
  expect_error(test_retrospective(sig_level = 0), sig_level_text)
  expect_error(test_retrospective(sig_level = "ciao"), sig_level_text)
  expect_error(test_retrospective(sig_level = c(.1,.2)), sig_level_text)

  # ratio_sd
  ratio_sd_text <- "Argument 'ratio_sd' has to be a single finite number grater than 0"
  expect_error(test_retrospective(ratio_sd = -1), ratio_sd_text)
  expect_error(test_retrospective(ratio_sd = Inf), ratio_sd_text)
  expect_error(test_retrospective(ratio_sd = 0), ratio_sd_text)
  expect_error(test_retrospective(ratio_sd = "ciao"), ratio_sd_text)
  expect_error(test_retrospective(ratio_sd = c(.1,.2)), ratio_sd_text)

  # B
  B_text <- "Argument 'B' has to be a single integer value grater than 1"
  expect_error(test_retrospective(B = 1), B_text)
  expect_error(test_retrospective(B = Inf), B_text)
  expect_error(test_retrospective(B = "ciao"), B_text)
  expect_error(test_retrospective(B = c(10,20)), B_text)

  # tl
  tl_text <- "Argument 'tl' has to be a single numeric value"
  expect_error(test_retrospective(tl = "ciao"), tl_text)
  expect_error(test_retrospective(tl = c(10,20)), tl_text)

  #tu
  tu_text <- "Argument 'tu' has to be a single numeric value"
  expect_error(test_retrospective(tu = "ciao"), tu_text)
  expect_error(test_retrospective(tu = c(10,20)), tu_text)

  #B_effect
  B_effect_text <- "Argument 'B_effect' has to be a single integer value grater than 1"
  expect_error(test_retrospective(B_effect = 1), B_effect_text)
  expect_error(test_retrospective(B_effect = Inf), B_effect_text)
  expect_error(test_retrospective(B_effect = "ciao"), B_effect_text)
  expect_error(test_retrospective(B_effect = c(10,20)), B_effect_text)

  # display_message
  display_text <- "Argument 'display_message' has to be logical"
  expect_error(test_retrospective(display_message = 1), display_text)
  expect_error(test_retrospective(display_message = Inf), display_text)
  expect_error(test_retrospective(display_message = NULL), display_text)
  expect_error(test_retrospective(display_message = "ciao"), display_text)
  expect_error(test_retrospective(display_message = c(.1,.2)), display_text)

  #----    Other cases    ----

  # # coherence effect_type and test_method
  # coherence_corr <- "If  'effect_type = correlation', argument 'test_method' has to be 'pearson'"
  # expect_error(test_retrospective(effect_type = "correlation", test_method = "paired"),
  #              coherence_corr)
  # coherence_cohen <- "No appropriate 'test_method' for 'effect_type = cohen_d'"
  # expect_error(test_retrospective(effect_type = "cohen_d", test_method = "pearson"),
  #              coherence_cohen)

  # correlation and sample_n2
  correlation_text <- "If 'test_method = pearson', argument 'ratio_n' is set to NULL"
  expect_warning(test_retrospective(sample_n2 = 30),
                 correlation_text)

  # one_sample and sample_n2
  one_sample_text <- "If 'test_method = one_sample', argument 'sample_n2' must be set to NULL"
  expect_error(test_retrospective(sample_n2 = 30,  test_method = "one_sample"),
               one_sample_text)

  # paired and sample_n2
  paired_text <- "If 'test_method = paired', arguments 'sample_n1' and 'sample_n2' must be equal"
  expect_error(test_retrospective(sample_n1 = 20, sample_n2 = 30,  test_method = "paired"),
               paired_text)
  expect_error(test_retrospective(sample_n1 = 20, sample_n2 = NULL,  test_method = "paired"),
               paired_text)

  # two_sample or welch and sample_n2
  t_test_text <- "Argument 'sample_n2' is required for the specified 'test_method'"
  expect_error(test_retrospective(sample_n2 = NULL,  test_method = "two_sample"),
               t_test_text)
  expect_error(test_retrospective(sample_n2 = NULL,  test_method = "welch"),
               t_test_text)

  # welch and ratio_sd
  t_test_ratio_text1 <- "Argument 'ratio_sd' is required only for 'test_method = welch'"
  t_test_ratio_text2 <- "Argument 'ratio_sd' can not be 1 for 'test_method = welch'\n  Consider 'test_method = two_sample' instead"
  expect_error(test_retrospective(sample_n2 = 20, ratio_sd = 1.5,  test_method = "two_sample"),
               t_test_ratio_text1)
  expect_error(test_retrospective(sample_n2 = 20, ratio_sd = 1,  test_method = "welch"),
               t_test_ratio_text2)
})


#----    obtain same results    ----

test_that("same results as previous run", {
  expect_known_value(with_seed(seed = 2020,
                               retrospective(sample_n1 = 10, effect_size = .3)$effect_info),
                            file = "test_cache/effect_info_single_cor")
  expect_known_value(with_seed(seed = 2020,
                               retrospective(sample_n1 = 10, effect_size = .3,  test_method = "one_sample")$effect_info),
                     file = "test_cache/effect_info_single_cohen")
  expect_known_value(with_seed(seed = 2020,
                               retrospective(sample_n1 = 10, effect_size = .3)$retrospective_res),
                            file = "test_cache/res_corr_single")
  expect_known_value(with_seed(seed = 2020,
                               retrospective(sample_n1 = 10, effect_size = .3,  test_method = "one_sample")$retrospective_res),
                     file="test_cache/res_cohen_single")
  expect_known_value(with_seed(seed = 2020,
                               retrospective(sample_n1 = 10, sample_n2 = 10, effect_size = .3,  test_method = "paired")$retrospective_res),
                     file="test_cache/res_cohen_paired")
  expect_known_value(with_seed(seed = 2020,
                               retrospective(sample_n1 = 10, sample_n2 = 20, effect_size = .3,  test_method = "two_sample")$retrospective_res),
                     file="test_cache/res_cohen_two_sample")
  expect_known_value(with_seed(seed = 2020,
                               retrospective(sample_n1 = 10, sample_n2 = 10, effect_size = .3,  test_method = "welch",
                                   ratio_sd = 1.5)$retrospective_res), file="test_cache/res_cohen_welch")


  expect_known_value(with_seed(seed = 2020,
                               retrospective(sample_n1 = 10, effect_size = function(x) rnorm(x), B = 10, B_effect = 10)$effect_info),
                     file = "test_cache/effect_info_dist")
  expect_known_value(with_seed(seed = 2020,
                               retrospective(sample_n1 = 10, effect_size = function(x) rnorm(x), B=100, B_effect = 10)$retrospective_res),
                     file = "test_cache/res_corr_dist")
  expect_known_value(with_seed(seed = 2020,
                               retrospective(sample_n1 = 10, effect_size = function(x) rnorm(x),  test_method = "one_sample", B=100, B_effect = 10)$retrospective_res),
                     file = "test_cache/res_cohen_dist")
  expect_known_value(with_seed(seed = 2020,
                               retrospective(sample_n1 = 10, sample_n2 = 10, effect_size = function(x) rnorm(x),  test_method = "welch",
                                   ratio_sd = 1.5, B=100, B_effect = 10)$retrospective_res), file = "test_cache/res_welch_dist")

  })


#----
