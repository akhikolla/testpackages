
context("blrm_exnex tests")


set.seed(123144)

eps <- 1E-4

## reference runs
single_agent  <- run_example("single_agent")
combo2  <- run_example("combo2")
combo3  <- run_example("combo3")

suppressPackageStartupMessages(library(dplyr))

## perform some basic checks on the shape of outputs
check_model_basic  <- function(fit, envir) {
  skip_on_cran()
  data <- fit$data
  num_obs  <- nrow(data)
  ss  <- summary(fit, interval_prob=c(0,0.16,0.33,1))

  expect_equal(nrow(ss), num_obs)
  expect_equal(ncol(ss), 8)
  expect_equal(sum(is.na(ss)), 0)
  ss2  <- summary(fit, interval_prob=c(0,1))
  expect_equal(sum((ss2[,"(0,1]"] - 1)<eps), nrow(data))

  ## check that the median and 50% interval prob align
  s_med  <- ss[,"50%"]
  for(i in 1:length(s_med)) {
    i <- 1
    q_med  <- s_med[i]
    s_q  <- summary(fit, interval_prob=c(0,q_med,1))
    expect_true(abs(s_q[i,6] - 0.5) < eps)
    expect_true(abs(s_q[i,7] - 0.5) < eps)
  }

  iter  <- getOption("OncoBayes2.MC.iter")
  warmup  <- getOption("OncoBayes2.MC.warmup")
  chains  <- getOption("OncoBayes2.MC.chains")
  num_sim <- chains * (iter - warmup)

  plin  <- posterior_linpred(fit)
  expect_equal(nrow(plin), num_sim)
  expect_equal(ncol(plin), nrow(data))
  expect_equal(sum(is.na(ss)), 0)

  envir$fit  <- fit
  envir$data1  <- data[1,,drop=FALSE]

  suppressMessages(capture.output(single_obs_fit <- with(envir, update(fit, data=data1))))
  ss1  <- summary(single_obs_fit)
  expect_equal(nrow(ss1), 1)
  expect_equal(ncol(ss1), 5)
  expect_equal(sum(is.na(ss1)), 0)

  plin1  <- posterior_linpred(single_obs_fit)
  expect_equal(nrow(plin1), num_sim)
  expect_equal(ncol(plin1), 1)
}


test_that("blrm_exnex data handling consistency single-agent",
          check_model_basic(single_agent$blrmfit, single_agent))

test_that("blrm_exnex data handling consistency combo2",
          check_model_basic(combo2$blrmfit, combo2))

test_that("blrm_exnex data handling consistency combo3",
          check_model_basic(combo3$blrmfit, combo3))

test_that("interval probabilites are consistent", {
  skip_on_cran()
  combo2_sens  <- combo2
  sens_low <- hist_combo2 %>%
    mutate(num_toxicities=0)
  combo2_sens$sens_low  <- sens_low
  outp <- suppressMessages(capture.output(sens_fit_low <- with(combo2_sens, update(blrmfit, data=sens_low))))
  s_low  <- summary(sens_fit_low, interval_prob=c(0,0.9,1))
  expect_equal(sum((s_low[,"(0,0.9]"] - 1)<eps), nrow(sens_low))
  expect_equal(sum((s_low[,"(0.9,1]"])<eps), nrow(sens_low))

  sens_high <- hist_combo2 %>%
    mutate(num_toxicities=num_patients)
  combo2_sens$sens_high <- sens_high
  outp <- suppressMessages(capture.output(sens_fit_high <- with(combo2_sens, update(blrmfit, data=sens_high))))
  s_high  <- summary(sens_fit_high, interval_prob=c(0,0.1,1))
  expect_equal(sum((s_high[,"(0.1,1]"] - 1)<eps), nrow(sens_low))
  expect_equal(sum((s_high[,"(0,0.1]"])<eps), nrow(sens_low))
})


## TODO: Add proper runs which are compared against published results,
## e.g. from the book chapter

test_that("interval probabilites are not NaN", {
  ss1  <- summary(combo2$blrmfit, interval_prob=c(0.1, 0.4))
  expect_true(all(!is.na(ss1[,6])))
  ss2  <- summary(combo2$blrmfit, transform=FALSE, interval_prob=logit(c(0.1, 0.4)))
  expect_true(all(!is.na(ss2[,6])))
})

test_that("correctness of numerical stable log1m_exp_max0", {
    b <- -1E-22
    expect_true(abs(OncoBayes2:::log1m_exp_max0(b - 1E-5) - log1p(-exp(b - 1E-5))) < 1E-10)
    expect_true(abs(OncoBayes2:::log1m_exp_max0(b - 1E-3) - log1p(-exp(b - 1E-3))) < 1E-10)
    expect_true(abs(OncoBayes2:::log1m_exp_max0(b - 1E-1) - log1p(-exp(b - 1E-1))) < 1E-10)
    expect_true(abs(OncoBayes2:::log1m_exp_max0(b - 1E-0) - log1p(-exp(b - 1E-0))) < 1E-10)
})

test_that("all expected posterior quantiles are returned", {

  prob <- 0.95
  f <- function(){
    summary(combo2$blrmfit, prob = prob)
  }
  expect_silent(f())
  ss1  <- f()
  p <- unique(c((1 - prob)/ 2, (1 - (1 - prob) / 2), 0.5))
  plab <- as.character(100 * p)
  expect_true(all(
    sapply(plab, function(lab){
      any(grepl(names(ss1), pattern = lab))
    })
  ))

  prob <- c(0.5, 0.95)
  expect_silent(f())
  ss2 <- f()
  p <- unique(c((1 - prob)/ 2, (1 - (1 - prob) / 2), 0.5))
  plab <- as.character(100 * p)
  expect_true(all(
    sapply(plab, function(lab){
      any(grepl(names(ss2), pattern = lab))
    })
  ))

})


query_fit <- function(fit) {
  capture_output(print(fit))
  prior_summary(fit)
  s1 <- summary(fit)
  s2 <- summary(fit, interval_prob=c(0,0.16,0.33,1.0))
  expect_true(ncol(s1) == ncol(s2)-3)
  pl <- posterior_linpred(fit)
  pi <- posterior_interval(fit)
  pp <- predictive_interval(fit)
}

test_that("blrmfit methods work for single agent models", {
  ## this test is successfull if methods run without errors as we do
  ## not have running examples given they are too costly.
  query_fit(single_agent$blrmfit)
})

test_that("blrmfit methods work for combo2 models", {
  ## this test is successfull if methods run without errors as we do
  ## not have running examples given they are too costly.
  query_fit(combo2$blrmfit)
})

test_that("blrmfit methods work for combo3 models", {
  ## this test is successfull if methods run without errors as we do
  ## not have running examples given they are too costly.
  query_fit(combo3$blrmfit)
})


test_that("blrm_exnex accepts single-stratum data sets with general prior definition", {

    num_comp <- 1 # one investigational drug
    num_inter <- 0 # no drug-drug interactions need to be modeled
    num_groups <- nlevels(hist_SA$group_id) # no stratification needed
    num_strata <- 1 # no stratification needed

    dref <- 50

    hist_SA_alt  <- mutate(hist_SA, stratum=factor(1))

    ## in case a single stratum is used in the data, the priors for tau should accept
    blrmfit <- blrm_exnex(
        cbind(num_toxicities, num_patients - num_toxicities) ~
            1 + log(drug_A / dref) |
            0 |
            stratum/group_id,
        data = hist_SA_alt,
        prior_EX_mu_mean_comp = matrix(
            c(logit(1/2), # mean of intercept on logit scale
              log(1)),    # mean of log-slope on logit scale
            nrow = num_comp,
            ncol = 2
        ),
        prior_EX_mu_sd_comp = matrix(
            c(2,  # sd of intercept
              1), # sd of log-slope
            nrow = num_comp,
            ncol = 2
        ),
        prior_EX_tau_mean_comp = array(0, c(1,num_comp,2)),
        prior_EX_tau_sd_comp   = array(1, c(num_comp,2)),
        prior_EX_prob_comp = matrix(1, nrow = num_comp, ncol = 1),
        prior_tau_dist = 0,
        prior_PD = FALSE,
        cores=1,
        iter=10,
        warmup=5
    )

    expect_true(nrow(summary(blrmfit)) == nrow(hist_SA_alt))

})


test_that("update.blrmfit grows the data set", {
    single_agent_new  <- single_agent
    single_agent_new$new_cohort_SA <- data.frame(group_id="trial_A", num_patients=4, num_toxicities=2, drug_A=50)
    single_agent_new$new_blrmfit_1 <- with(single_agent_new, update(blrmfit, add_data=new_cohort_SA))
    expect_true(nrow(summary(single_agent_new$new_blrmfit_1)) == nrow(hist_SA)+1)

    ## ensure that the data accumulates
    new_blrmfit_2 <- with(single_agent_new, update(new_blrmfit_1, add_data=new_cohort_SA))
    expect_true(nrow(summary(new_blrmfit_2)) == nrow(hist_SA)+2)

    combo2_new  <- combo2
    combo2_new$new_codata <- data.frame(group_id=c("IIT", "trial_AB"),
                                        drug_A=c(8, 8),
                                        drug_B=c(800, 900),
                                        drug_C=c(10, 20),
                                        num_patients=10, num_toxicities=2, stringsAsFactors = TRUE)

    ## this one will fail due to a factor levels mismatch
    expect_error(new_blrmfit_3 <- with(combo2_new, update(blrmfit, add_data=new_codata)),
                 "Mismatch in factor level defintion of grouping", fixed=TRUE)

    combo2_new$new_codata  <-  mutate(combo2_new$new_codata,
                                      group_id=as.character(group_id))
    set.seed(123144)
    new_blrmfit_3 <- with(combo2_new, update(blrmfit, add_data=new_codata))
    expect_true(nrow(summary(new_blrmfit_3)) == nrow(codata_combo2)+2)
    ## Test that adding dummy data does not change results in the other rows
    set.seed(123144)
    combo2_new_with_dummy <- combo2
    combo2_new_with_dummy$new_codata <- add_row(combo2_new_with_dummy$new_codata, group_id = factor("IIT"), drug_A = 1, drug_B = 1, drug_C = 1, num_patients = 0, num_toxicities = 0)

    new_blrmfit_3_with_dummy <- with(combo2_new_with_dummy, update(blrmfit, add_data=new_codata))
    expect_equal(nrow(summary(new_blrmfit_3)) + 1, nrow(summary(new_blrmfit_3_with_dummy)))

    ## test if the log-likelihood is the same for parameter-vector 0
    ## on unconstrained space
    num_pars  <- rstan::get_num_upars(new_blrmfit_3$stanfit)
    theta_uconst <- rep(0.0, num_pars)
    log_prob_group <- rstan::log_prob(new_blrmfit_3$stanfit, theta_uconst, gradient=FALSE)
    log_prob_group_and_dummy  <- rstan::log_prob(new_blrmfit_3_with_dummy$stanfit, theta_uconst, gradient=FALSE)
    expect_equal(log_prob_group, log_prob_group_and_dummy)


    ## Same for empty group
    set.seed(123144)
    new_blrmfit_with_empty_group <- with(combo2_new, update(blrmfit, data=blrmfit$data))
    set.seed(123144)
    new_blrmfit_with_empty_group_and_dummy <- with(combo2_new, update(blrmfit, data=add_row(blrmfit$data, group_id = factor("IIT"), drug_A = 1, drug_B = 1, num_patients = 0, num_toxicities = 0)))
    expect_equal(nrow(summary(new_blrmfit_with_empty_group)) + 1, nrow(summary(new_blrmfit_with_empty_group_and_dummy)))
    ##summary(new_blrmfit_with_empty_group) - summary(new_blrmfit_with_empty_group_and_dummy)[1:nrow(summary(new_blrmfit_with_empty_group)),]
    ## change level order and test

    ## test if the log-likelihood is the same for parameter-vector 0
    ## on unconstrained space
    num_pars  <- rstan::get_num_upars(new_blrmfit_with_empty_group$stanfit)
    theta_uconst <- rep(0.0, num_pars)
    log_prob_prob_ref  <- rstan::log_prob(combo2_new$blrmfit$stanfit, theta_uconst, gradient=FALSE)
    log_prob_empty_group <- rstan::log_prob(new_blrmfit_with_empty_group$stanfit, theta_uconst, gradient=FALSE)
    log_prob_empty_group_and_dummy  <- rstan::log_prob(new_blrmfit_with_empty_group_and_dummy$stanfit, theta_uconst, gradient=FALSE)
    expect_equal(log_prob_prob_ref, log_prob_empty_group)
    expect_equal(log_prob_prob_ref, log_prob_empty_group_and_dummy)

    combo2_new$flaky_new_codata  <-  mutate(combo2_new$new_codata,
                                            group_id=NULL)

    expect_error(new_blrmfit_4 <- with(combo2_new, update(blrmfit, add_data=flaky_new_codata)),
                 "Assertion on 'grouping and/or stratum columns'.*")

    with(combo2_new, {
        wrong_codata  <- new_codata
        levels_wrong  <- levels(codata_combo2$group_id)
        levels_wrong[1:2] <- levels_wrong[2:1]
        wrong_codata  <- mutate(wrong_codata, group_id=factor(as.character(group_id), levels=levels_wrong))
    })
    expect_true(sum(levels(combo2_new$wrong_codata$group_id) != levels(codata_combo2$group_id)) == 2)

    expect_error(new_blrmfit_5 <- with(combo2_new, update(blrmfit, add_data=wrong_codata)),
                 "Mismatch in factor level defintion of grouping", fixed=TRUE)

})

test_that("update.blrmfit does regular updating", {
    single_agent_new  <- single_agent
    single_agent_new$only_cohort_SA <- data.frame(group_id="trial_A", num_patients=4, num_toxicities=2, drug_A=50, stringsAsFactors = TRUE)
    single_agent_new$only_blrmfit_1 <- with(single_agent_new, update(blrmfit, data=only_cohort_SA))
    expect_true(nrow(summary(single_agent_new$only_blrmfit_1)) == 1)
})

test_that("update.blrmfit combines data and add_data", {
    single_agent_new  <- single_agent
    single_agent_new$only_cohort_SA <- data.frame(group_id="trial_A", num_patients=4, num_toxicities=2, drug_A=50)
    single_agent_new$hist_SA_sub <- hist_SA[1:3,]
    single_agent_new$only_blrmfit_1 <- with(single_agent_new, update(blrmfit, data=hist_SA_sub, add_data=only_cohort_SA))
    expect_true(nrow(summary(single_agent_new$only_blrmfit_1)) == 4)
})
