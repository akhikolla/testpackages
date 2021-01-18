context("Test basic Stan models")

# only for debugging
if (FALSE) {
  library(rstan)
  library(testthat)
  library(gastempt)
  library(assertthat)
  library(gastempt)
}

test_that("stanmodels exist", {
  expect_true(exists("stanmodels"), "No stanmodels found")
  mn = stan_model_names()
  expect_gt(nrow(mn), 3)
  expect_equal(names(mn), c("model_name", "description"))
})

test_that("stanmodels$linexp_gastro_1b exists", {
  mod = stanmodels$linexp_gastro_1b
  expect_s4_class(mod, "stanmodel")
  expect_identical(mod@model_name, "linexp_gastro_1b")
  #message(str(mod))
})

gastempt_data = function(){
  simulate_gastempt(n_records = 6, seed = 471)$stan_data
}


run_precompiled_model = function(model, iter = 500){
  cat("\nTesting ", model,"\n") # Generate output for travis
  mod = gastempt:::stanmodels[[model]]
  testthat::expect_s4_class(mod, "stanmodel")
  testthat::expect_identical(mod@model_name, model)
  data = gastempt_data()
  data$lkj = 2
  data$student_df = 5
  cap = capture_output({
    fit = suppressWarnings(
      rstan::sampling(mod, data = data, chains = 2, iter = iter,
                      refresh = -1, verbose = FALSE))
  })
  expect_is(fit, "stanfit")
}

test_that("Running precompiled models linexp _1x directly returns valid result", {
  iter = ifelse(identical(Sys.getenv("NOT_CRAN"), "true"), 500, 100)
  run_precompiled_model("linexp_gastro_1b", iter = iter)
})


test_that("Running precompiled powexp models directly returns valid result", {
  run_precompiled_model("powexp_gastro_1b")
})

test_that("Running precompiled powexp models directly returns valid result", {
  skip_on_cran()
  run_precompiled_model("powexp_gastro_2c")
})

test_that("Running stan_gastempt fit with default parameters returns valid result", {
  skip_on_cran()
  d = simulate_gastempt(n_records = 6, seed = 471)
  v0_d = d$rec$v0
  chains = 2 # Problems with more chains on travis
  options(mc.cores = min(parallel::detectCores(), chains))
  ret = stan_gastempt(d$data, chains = chains,  iter = 500, refresh = -1)
  expect_is(ret, "stan_gastempt")
  expect_is(ret$plot, "ggplot")
  expect_is(ret$plot, "ggplot")
  expect_equal(ret$coef, coef(ret))
  expect_s4_class(ret$fit, "stanfit")
  # residual standard deviation
  v0_f = ret$coef$v0
  expect_lt(sqrt(var(v0_d - v0_f)), 8)
  expect_true(all(c("sigma", "mu_kappa", "sigma_kappa", "lp") %in%
                    names(attributes(ret$coef))))
})


test_that("Running stan_gastempt with powexp returns valid result", {
  skip_on_cran()
  options(mc.cores = 1)
  d = simulate_gastempt(n_records = 6, seed = 471, model = powexp,
                        beta_mean = 2.5, missing = 0.3)
  v0_d = d$rec$v0
  ret = stan_gastempt(d$data, model_name = "powexp_gastro_2c",
                      refresh = -1, iter = 500)
  expect_is(ret, "stan_gastempt")
  expect_is(ret$plot, "ggplot")
  expect_s4_class(ret$fit, "stanfit")
  # residual standard deviation
  v0_f = ret$coef$v0
  expect_lt(sqrt(var(v0_d - v0_f)), 20)
  expect_true(all(c("sigma", "mu_beta", "sigma_beta", "lp") %in%
                    names(attributes(ret$coef))))
})

test_that("Running stan_gastempt fit with non-default parameters returns valid result", {
  skip_on_travis()
  skip_on_cran()
  d = simulate_gastempt(n_records = 6, seed = 471)
  v0_d = d$rec$v0
  ret = stan_gastempt(d$data, model_name = "linexp_gastro_2b",
                      refresh = -1, chains = 2, init_r = 0.3)
  expect_is(ret, "stan_gastempt")
  v0_f = ret$coef$v0
  expect_lt(sqrt(var(v0_d - v0_f)), 8)
})

test_that("Running stan_gastempt with many missing data returns valid result", {
  skip_on_cran()
  d = simulate_gastempt(n_records = 6, missing = 0.3, seed = 471)
  v0_d = d$rec$v0
  ret = stan_gastempt(d$data, model_name = "linexp_gastro_1b", refresh = -1)
  expect_is(ret, "stan_gastempt")
  expect_is(ret$plot, "ggplot")
  expect_s4_class(ret$fit, "stanfit")
  # residual standard deviation
  v0_f = ret$coef$v0
  expect_lt(sqrt(var(v0_d - v0_f)), 8)
  expect_true(all(c("sigma", "mu_kappa", "sigma_kappa", "lp") %in%
                    names(attributes(ret$coef))))
})

test_that("Direct use of sample model returns valid results", {
  skip("Slow. Only use on errors in other Stan functions.")
  data = gastempt_data()
  stan_model = "../../inst/stan/linexp_gastro_2b.stan"
  print(dir("../../inst/stan/"))
  expect_true(file.exists(stan_model))
  rstan_options(auto_write = TRUE)
  iter = 500
  cap = capture_output({
    mr_b =  suppressWarnings(
      stan(stan_model, data = data, chains = 4, iter = iter,
           seed = 4711, refresh = FALSE))
  })
  expect_is(mr_b, "stanfit")


  stan_model = "../../inst/stan/linexp_gastro_1b.stan"
  expect_true(file.exists(stan_model))
  mr_d = stan(stan_model, data = data, chains = 4, iter = iter,
              seed = 4711, refresh = FALSE)

  m_d = get_posterior_mean(mr_d)[,5];
  m_c = get_posterior_mean(mr_c)[,5];
  m_b = get_posterior_mean(mr_b)[,5];
  cbind(m_b, m_c, m_c)

  sum(get_elapsed_time(mr_b))
  sum(get_elapsed_time(mr_c))
  sum(get_elapsed_time(mr_d))
})

