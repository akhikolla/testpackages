###########################################
####    Test methods design analysis   ####
###########################################

library(PRDA)


#----    print    ----

context("Evaluate methods for design analysis")

test_that("evaluate print.design_analysis", {
  # Retrospective
  verify_output(path = "test_cache/print_retro_cohen.txt",
                with_seed(seed = 2020,
                          retrospective(effect_size = function(x) rnorm(x, .3, .1),
                              sample_n1 = 20, sample_n2 = 20, test_method = "paired",
                              tl = .15, tu = .45, B_effect = 10, B = 50)))
  verify_output(path = "test_cache/print_retro_cor.txt",
                with_seed(seed = 2020,
                          retrospective(effect_size = .3, sample_n1 = 20, B = 50)))

  # Prospective
  verify_output(path = "test_cache/print_pro_cor.txt",
                with_seed(seed = 2020,
                          prospective(effect_size = function(x) rnorm(x, .3, .1), power = .8,
                            tl = .15, tu = .45, B_effect = 10, B = 50)))
  verify_output(path = "test_cache/print_pro_cohen.txt",
                with_seed(seed = 2020,
                          prospective(effect_size = .3, power = .8, ratio_n = 2,
                            test_method = "two_sample", B = 50)))
})


#----    summary    ----

test_that("evaluate summary.design_analysis", {
  # Retrospective
  verify_output(path = "test_cache/summary_retro_cor.txt",
                with_seed(seed = 2020,
                          summary(retrospective(effect_size = .3, sample_n1 = 20, B = 50))))
  verify_output(path = "test_cache/summary_retro_cohen.txt",
                with_seed(seed = 2020,
                          summary(retrospective(effect_size = function(x) rnorm(x, .3, .1),
                              sample_n1 = 20, sample_n2 = 20, test_method = "paired",
                              tl = .15, tu = .45, B_effect = 10, B = 50))))

  # Prospective
  verify_output(path = "test_cache/summary_pro_cor.txt",
                with_seed(seed = 2020,
                          summary(prospective(effect_size = function(x) rnorm(x, .3, .1),
                            power = .8, tl = .15, tu = .45, B_effect = 10, B = 50))))
  verify_output(path = "test_cache/summary_pro_cohen.txt",
                with_seed(seed = 2020,
                          summary(prospective(effect_size = .3, power = .8, ratio_n = 2,
                            test_method = "two_sample", B = 50))))
})

#----


