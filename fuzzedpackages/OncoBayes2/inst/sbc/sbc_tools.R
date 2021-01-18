#'
#' Utilities for SBC validation
#'
library(devtools)
devtools::load_all("../..")
library(rstan)
library(mvtnorm)
library(checkmate) # only needed inside blrm_exnex
library(Formula) # only needed inside blrm_exnex
library(abind) # only needed inside blrm_exnex
library(dplyr)
library(tidyr)
library(assertthat)
source("lkj.R")

assert_that(Sys.getenv("CMDSTAN") != "", msg="CMDSTAN environment variable must be set.")

options(OncoBayes2.abbreviate.min = 0)

## Sample prior

## First sample EX hyperparameters
## prior_EX_mu_mean_comp + prior_EX_mu_sd_comp => EX_mu_comp (group means EX)
## prior_EX_tau_mean_comp + prior_EX_tau_sd_comp => EX_tau_comp (group taus EX)
## prior_EX_corr_eta_comp => rho (only 0 for now)

## prior_EX_mu_mean_inter + prior_EX_mu_sd_inter => EX_mu_inter (group means EX)
## prior_EX_tau_mean_inter + prior_EX_tau_sd_inter => EX_tau_comp (group taus EX)
## prior_EX_corr_eta_inter => rho (only 0 for now)
## => EX_eta (group specific means for inter (mvn normal))

## Then sample group parameters, EX
## => log_beta (group specific means for comp (mvn normal))
## => eta (group specific means for comp (mvn normal))

## and NEX parameters are in the same data structure with index starting at num_groups+1

## prior_NEX_mu_mean_comp + prior_NEX_mu_sd_comp => NEX_beta = NEX_comp (group means NEX iid groupwise)
## prior_NEX_mu_mean_inter + prior_NEX_mu_sd_inter => NEX_eta = NEX_inter (group means NEX iid groupwise)

## prior_is_EXNEX_comp + prior_EX_prob_comp => pick which one
## prior_is_EXNEX_inter + prior_EX_prob_inter => pick which one

sample_prior <- function(model) {

    num_strata <- model$num_strata
    num_groups <- model$num_groups
    num_comp <- model$num_comp
    num_inter <- model$num_inter
    blrm_args  <- model$blrm_args
    standata <- model$base_fit$standata
    group_stratum <- standata$group_stratum_cid

    ## group specific parameters: EX, then NEX
    log_beta <- array(NA, dim=c(2*num_groups, num_comp, 2))
    eta <- array(NA, dim=c(2*num_groups, num_inter))

    ## sample EX hyperparameters

    ## mu
    EX_mu_comp <- array(NA, dim=c(num_comp, 2))
    EX_mu_inter <- array(NA, dim=c(num_inter))

    for (j in 1:num_comp) {
        EX_mu_comp[j,1] <- rnorm(1, standata$prior_EX_mu_mean_comp[j,1], standata$prior_EX_mu_sd_comp[j,1])
        EX_mu_comp[j,2] <- rnorm(1, standata$prior_EX_mu_mean_comp[j,2], standata$prior_EX_mu_sd_comp[j,2])
    }

    if (num_inter > 0) {
        for (j in 1:num_inter) {
            EX_mu_inter[j] <- rnorm(1, standata$prior_EX_mu_mean_inter[j], standata$prior_EX_mu_sd_inter[j])
        }
    }

    ## tau
    EX_tau_comp <- array(NA, dim=c(num_strata, num_comp, 2))
    EX_tau_inter <- array(NA, dim=c(num_strata, num_inter))
    ## correlation matrix
    EX_corr_comp <- array(NA, dim=c(num_strata, num_comp, 2, 2))
    EX_corr_inter <- array(NA, dim=c(num_strata, num_inter, num_inter))
    EX_Sigma_comp <- array(NA, dim=c(num_strata, num_comp, 2, 2))
    EX_Sigma_inter <- array(NA, dim=c(num_strata, num_inter, num_inter))

    sample_tau_prior <- function(dist, a, b) {
        if(dist == 0)
            return(a)
        if(dist == 1)
            return(rlnorm(1, a, b))
        if(dist == 2)
            return(abs(rnorm(1, a, b)))
        stop("Unsupported tau prior density.")
    }

    for (s in 1:num_strata) {
        for (j in 1:num_comp) {
            EX_tau_comp[s,j,1] <- sample_tau_prior(standata$prior_tau_dist,
                                                   standata$prior_EX_tau_mean_comp[s,j,1],
                                                   standata$prior_EX_tau_sd_comp[s,j,1] )
            EX_tau_comp[s,j,2] <- sample_tau_prior(standata$prior_tau_dist,
                                                   standata$prior_EX_tau_mean_comp[s,j,2],
                                                   standata$prior_EX_tau_sd_comp[s,j,2] )
            EX_corr_comp[s,j,,] <- rcorvine(2, standata$prior_EX_corr_eta_comp[j], FALSE)
            EX_Sigma_comp[s,j,,] <- diag(as.vector(EX_tau_comp[s,j,]), 2, 2) %*% matrix(EX_corr_comp[s,j,,,drop=FALSE], 2, 2) %*% diag(as.vector(EX_tau_comp[s,j,]), 2, 2)
        }
        if (num_inter > 0) {
            for (j in 1:num_inter) {
                EX_tau_inter[s,j] <- sample_tau_prior(standata$prior_tau_dist,
                                                      standata$prior_EX_tau_mean_inter[s,j],
                                                      standata$prior_EX_tau_sd_inter[s,j] )
            }
            EX_corr_inter[s,,] <- diag(num_inter)
            if (num_inter > 1) {
                EX_corr_inter[s,,] <- rcorvine(num_inter, standata$prior_EX_corr_eta_inter, FALSE)
            }
            EX_Sigma_inter[s,,] <- diag(as.vector(EX_tau_inter[s,,drop=FALSE]), num_inter, num_inter) %*% matrix(EX_corr_inter[s,,,drop=FALSE], num_inter, num_inter) %*% diag(as.vector(EX_tau_inter[s,,drop=FALSE]), num_inter, num_inter)
        }
    }

  ## EX - group-specific parameters
  for (g in 1:num_groups) {
    s <- group_stratum[g]
    for (j in 1:num_comp) {
        log_beta[g,j,1:2] <- rmvnorm(1, EX_mu_comp[j,], EX_Sigma_comp[s,j,,])
        ##log_beta[g,j,1] <- rnorm(1, EX_mu_comp[j,1], EX_tau_comp[s,j,1] )
        ##log_beta[g,j,2] <- rnorm(1, EX_mu_comp[j,2], EX_tau_comp[s,j,2] )
        ##assert_that(standata$prior_EX_corr_eta_comp[j] == 1, msg="LKJ correlation == 1 is only supported.")
    }
    if (num_inter > 0) {
        if(num_inter > 1) {
            eta[g,] <- rmvnorm(1, EX_mu_inter, EX_Sigma_inter[s,,])
        } else {
            eta[g,1] <- rnorm(1, EX_mu_inter, EX_tau_inter[s,1])
        }
      ##for (j in 1:num_inter) {
      ##  eta[g,j] <- rnorm(1, EX_mu_inter[j], EX_tau_inter[s,j] )
      ##}
      ##assert_that(standata$prior_EX_corr_eta_inter == 1, msg="LKJ correlation == 1 is only supported.")
    }
  }
  ## NEX - group-specific parameters
  for (g in 1:num_groups) {
    for (j in 1:num_comp) {
      log_beta[num_groups+g,j,1] <- rnorm(1, standata$prior_NEX_mu_mean_comp[j,1], standata$prior_NEX_mu_sd_comp[j,1])
      log_beta[num_groups+g,j,2] <- rnorm(1, standata$prior_NEX_mu_mean_comp[j,2], standata$prior_NEX_mu_sd_comp[j,2])
    }
    if (num_inter > 0) {
      for (j in 1:num_inter) {
        eta[num_groups+g,j] <- rnorm(1, standata$prior_NEX_mu_mean_inter[j], standata$prior_NEX_mu_sd_inter[j])
      }
    }
  }

  ## convert slope to natural scale (enforced positivity)
  beta <- log_beta
  for (g in 1:(2*num_groups)) {
    for (j in 1:num_comp) {
      beta[g,j,2] <- exp(beta[g,j,2])
    }
  }

  ## sample EX / NEX membership
  is_EX_comp <- array(NA, dim=c(num_groups,num_comp))
  is_EX_inter <- array(NA, dim=c(num_groups,num_inter))
  draw_beta  <- array(NA, dim=c(1, num_groups, num_comp, 2))
  draw_eta  <- array(NA, dim=c(1, num_groups, num_inter))
  for(g in 1:num_groups) {
    for (j in 1:num_comp) {
      if(standata$prior_is_EXNEX_comp[j] == 1) {
        is_EX_comp[g,j] <- rbinom(1, 1, standata$prior_EX_prob_comp[g,j])
      } else {
        is_EX_comp[g,j] <- 1
      }
      gidx <- ifelse(is_EX_comp[g,j] == 1, g, num_groups + g)
      draw_beta[1,g,j,1]  <- beta[gidx,j,1]
      draw_beta[1,g,j,2]  <- beta[gidx,j,2]
    }
    if (num_inter > 0) {
      for (j in 1:num_inter) {
        if(standata$prior_is_EXNEX_inter[j] == 1) {
          is_EX_inter[g,j] <- rbinom(1, 1, standata$prior_EX_prob_inter[g,j])
        } else {
          is_EX_inter[g,j] <- 1
        }
        gidx <- ifelse(is_EX_inter[g,j] == 1, g, num_groups + g)
        draw_eta[1,g,j] <- eta[gidx,j]
      }
    }
  }

    ## name the array indices accordingly using the prior_summary
    ## structures
    ps  <- prior_summary(model$base_fit)
    dimnames(EX_mu_comp) <- dimnames(ps$EX_mu_log_beta)[c(2,1)]
    dimnames(EX_tau_comp) <- dimnames(ps$EX_tau_log_beta)[c(2,3,1)]

    dimnames(EX_mu_inter) <- dimnames(ps$EX_mu_eta)[c(1)]
    dimnames(EX_tau_inter) <- dimnames(ps$EX_tau_eta)[c(2,1)]

    dimnames(is_EX_comp) <- dimnames(ps$EX_prob_comp)
    dimnames(draw_beta) <- c(list(NULL), dimnames(ps$EX_prob_comp), list(coefficient=c("intercept", "log_slope")))

    dimnames(is_EX_inter)  <- dimnames(ps$EX_prob_inter)
    dimnames(draw_eta)  <- c(list(NULL), dimnames(ps$EX_prob_inter))

    list(draw_beta=draw_beta,
         draw_eta=draw_eta,
         EX_mu_comp=EX_mu_comp,
         EX_mu_inter=EX_mu_inter,
         EX_tau_comp=EX_tau_comp,
         EX_tau_inter=EX_tau_inter,
         EX_corr_comp=EX_corr_comp,
         EX_corr_inter=EX_corr_inter,
         log_beta=log_beta,
         eta=eta,
         is_EX_comp=is_EX_comp,
         is_EX_inter=is_EX_inter
         )
}

#'
#' Simulates a draw from the prior and fake data for it. This will be
#' the data generating step in the simulation. The function recieves
#' the problem data, job specifics and a blrmfit object which defines
#' the prior to sample and the design matrix.
#'
simulate_fake <- function(data, job, model, ...) {

  model  <- data$models[[model]]

  prior_draw  <- sample_prior(model)

  standata  <- model$base_fit$standata

  ## logit by data-row
  draw_mu <- with(standata, blrm_logit_grouped_vec(group, stratum, X_comp, X_inter, prior_draw$draw_beta, prior_draw$draw_eta))

  num_trials <- standata$r + standata$nr

  yrep <- rbinom(length(num_trials), num_trials, inv_logit(draw_mu))

  list(yrep = yrep, draw = prior_draw)

}

extract_draws <- function(sims, draw) lapply(sims, asub, idx=draw, dim=1, drop=FALSE)

extract_draw <- function(sims, draw) {
    assert_that(length(draw) == 1)
    lapply(lapply(sims, asub, idx=draw, dim=1, drop=FALSE), adrop, drop=1, one.d.array=TRUE)
}

restore_draw_dims <- function(standata, draw) {
    num_comp <- standata$num_comp
    num_inter <- standata$num_inter
    num_strata <- standata$num_strata
    num_groups <- standata$num_groups

    draw$mu_log_beta <- array(draw$mu_log_beta, c(num_comp,2))
    draw$tau_log_beta_raw <- array(draw$tau_log_beta_raw, c(num_strata,num_comp,2))
    draw$L_corr_log_beta <- array(draw$L_corr_log_beta, c(num_comp,2,2))
    draw$log_beta_raw <- array(draw$log_beta_raw, c(2*num_groups, num_comp, 2))

    if(num_inter != 0) {
        draw$eta_raw <- array(draw$eta_raw, c(2*num_groups, num_inter))
        draw$mu_eta <- array(draw$mu_eta, c(num_inter))
        draw$tau_eta_raw <- array(draw$tau_eta_raw, c(num_strata,num_inter))
        draw$L_corr_eta <- matrix(draw$L_corr_eta, num_inter, num_inter)
    } else {
        draw$eta_raw <- array(0, c(2*num_groups, num_inter))
        draw$mu_eta <- array(0, c(num_inter))
        draw$tau_eta_raw <- array(0, c(num_strata,num_inter))
        draw$L_corr_eta <- matrix(1, num_inter, num_inter)
    }

    draw
}

#' extracts from a given fit the mass matrix, stepsize and a draw from
#' the typical set. The warmup info from multiple chains is being
#' averaged together to obtain less noisy estimates.
learn_warmup_info <- function(standata, stanfit) {
    gmean <- function(x) exp(mean(log(x)))
    draw  <- extract_draw(rstan::extract(stanfit)[1:8], 1)
    warmup_info  <- extract_warmup_info(stanfit)
    warmup_info$stepsize  <- gmean(warmup_info$stepsize)
    warmup_info$inv_metric  <- apply(warmup_info$inv_metric, 1, gmean)
    c(warmup_info, list(draw=restore_draw_dims(standata, draw)))
}

#'
#'
#' Procedure to fit each fake data set using our fitting
#' procedure. This method obtains the problem data, job details and an
#' **instance** of the scenario as generated by `simulate_fake`.
#'

fit_exnex <- function(data, job, instance, ..., save_fit=FALSE) {

    yrep <- instance$yrep
    draw <- instance$draw
    group_draws <- list()
    group_draws$draw_beta <- array(draw$draw_beta, dim=dim(draw$draw_beta)[-1], dimnames=dimnames(draw$draw_beta)[-1])
    group_draws$draw_eta <- array(draw$draw_eta, dim=dim(draw$draw_eta)[-1], dimnames=dimnames(draw$draw_eta)[-1])

    pars <- job$pars$prob.pars
    model  <- data$models[[pars$model]]

    dref <- model$dref
    sim_data <- model$base_fit$data
    sim_data$num_toxicities <- yrep

    blrm_args <- model$blrm_args

    have_warmup_info  <- c("warmup_info") %in% names(model)

    if(have_warmup_info) {
        ## use a randomly selected warmup info from the ones provided
        fit_warmup_info <- sample(model$warmup_info, 1)[[1]]
        blrm_args <- model$blrm_args_with_warmup_info
        blrm_args$init  <- rep(list(fit_warmup_info$draw), blrm_args$chains)
        blrm_args$control <- modifyList(blrm_args$control,
                                        list(##adapt_inv_metric=fit_warmup_info$inv_metric,
                                             stepsize=fit_warmup_info$stepsize)
                                        )
    }

    fit <- update(model$base_fit,
                  data = sim_data,
                  init = blrm_args$init,
                  iter = blrm_args$iter,
                  warmup = blrm_args$warmup,
                  chains = blrm_args$chains,
                  control = blrm_args$control,
                  verbose=TRUE
                  )

    sampler_params <- rstan::get_sampler_params(fit$stanfit, inc_warmup=FALSE)
    n_divergent <- sum(sapply(sampler_params, function(x) sum(x[,'divergent__'])) )

    ##params <- c("mu_log_beta", "tau_log_beta", "beta_group", "mu_eta", "tau_eta", "eta_group")
    params <- c("mu_log_beta", "tau_log_beta", "beta_group")
    if(fit$has_inter)
        params <- c(params, "mu_eta", "tau_eta", "eta_group")
    fit_sum <- rstan::summary(fit$stanfit)$summary
    samp_diags <- fit_sum[apply(sapply(params, grepl, x = rownames(fit_sum)), 1, any), c("n_eff", "Rhat")]
    min_Neff <- ceiling(min(samp_diags[, "n_eff"], na.rm=TRUE))
    max_Rhat <- max(samp_diags[, "Rhat"], na.rm=TRUE)

    post <- rstan::extract(fit$stanfit, pars = params, inc_warmup = FALSE)

    lp_ess  <- as.numeric(rstan::monitor(as.array(fit$stanfit, pars="lp__"), print=FALSE)[1, c("Bulk_ESS", "Tail_ESS")])

    post_thin <- lapply(post, function(A) {
        assert_that(dim(A)[1] > 1023)
        asub(A, idx = seq(1, dim(A)[1], length = 1024-1), dims = 1, drop = FALSE)
    })

    calc_rank  <- function(sample, draw) {
        sdims  <- dim(sample)
        assert_that(all(sdims[-1] == dim(draw)))
        draw_margins  <- 2:length(sdims)
        res <- array(apply(sweep(sample, draw_margins, draw) < 0, draw_margins, sum), dim=sdims[-1])
        dimnames(res)  <- dimnames(draw)
        res
    }

    rank1 <- mapply(calc_rank,
                    post_thin[params[1:3]],
                    c(draw[c("EX_mu_comp", "EX_tau_comp")], list(beta_group=group_draws$draw_beta)),
                    SIMPLIFY=FALSE)

    if(fit$has_inter) {
        rank1 <- c(rank1,
                   mapply(calc_rank,
                          post_thin[params[4:6]],
                          c(draw[c("EX_mu_inter", "EX_tau_inter")], list(eta_group=group_draws$draw_eta)),
                          SIMPLIFY=FALSE)
                   )
    }

    flatten_array  <- function(data, var) {
        if(missing(var))
            var  <- deparse(substitute(data))
        idx <- expand.grid(lapply(dim(data), seq))
        num_dim <- length(dim(data))
        size  <- nrow(idx)
        values  <- sapply(1:size, function(r) { asub(data, idx[r,], drop=FALSE) } )
        idx[1:num_dim] <- lapply(1:num_dim, function(col) dimnames(data)[[col]][idx[,col]])
        names(values) <- paste0(var, "[", do.call("paste", c(idx[1:num_dim], list(sep=","))), "]")
        res <- matrix(values, nrow=1, ncol=size)
        colnames(res) <- names(values)
        as.data.frame(res)
    }

    rank_wide <- bind_cols(mapply(flatten_array, rank1, names(rank1), SIMPLIFY=FALSE))

    res <- list(rank = rank_wide,
                min_Neff = min_Neff,
                n_divergent = n_divergent,
                max_Rhat=max_Rhat,
                lp_ess_bulk = lp_ess[1],
                lp_ess_tail = lp_ess[2]
                )

    if(save_fit)
        res$fit  <- fit

    if(!have_warmup_info) {
        res <- c(res, learn_warmup_info(fit$standata, fit$stanfit))
    }

    return(res)
}


extract_warmup_info <- function(fit) {
    info  <- sapply(get_adaptation_info(fit), strsplit, "\n")
    ex_stepsize <- function(chain_info) {
        stepsize_line <- which(grepl("Step size", chain_info))
        as.numeric(strsplit(chain_info[stepsize_line], " = ")[[1]][2])
    }
    ex_mass <- function(chain_info) {
        metric_line <- which(grepl("inverse mass matrix", chain_info)) + 1
        as.numeric(strsplit(sub("^#", "", chain_info[metric_line]), ", ")[[1]])
    }
    stepsize <- sapply(info, ex_stepsize)
    inv_metric <- do.call(cbind, lapply(info, ex_mass))
    colnames(inv_metric) <- names(stepsize) <- paste0("chain_", 1:length(info))
    list(stepsize=stepsize, inv_metric=inv_metric)
}



# AB: not currently using this function from rbest SBC...
# scale_ranks <- function(Nbins, scale=1) {
#   ## scale must evenly divide the total number of bins
#   assert_that(round(Nbins/scale) == Nbins/scale)
#   breaks <- (0:(Nbins/scale))
#   Nbreaks <- length(breaks)
#   function(scen) {
#     vars <- grep("^rank.", names(scen), value=TRUE)
#     res <- lapply(vars, function(v) hist(ceiling((scen[[v]]+1)/scale), breaks=breaks, plot=FALSE, include.lowest=FALSE)$counts)
#     names(res) <- gsub("^rank", "count", vars)
#     res$rank <- breaks[-Nbreaks]
#     res <- as.data.frame(do.call(cbind, res))
#     res
#   }
# }


## Submits to batchtools cluster with fault tolerance, i.e.
## resubmitting failed jobs max_num_tries times
auto_submit <- function(jobs, registry, resources=list(), max_num_tries = 10) {
  all_unfinished_jobs <- jobs

  num_unfinished_jobs <- nrow(all_unfinished_jobs)
  num_all_jobs <- num_unfinished_jobs
  remaining_tries <- max_num_tries
  all_jobs_finished <- FALSE
  while (remaining_tries > 0 && !all_jobs_finished) {
    remaining_tries <- remaining_tries - 1

    message("Submitting jobs at ", Sys.time())
    # Once things run fine let's submit this work to the cluster.
    submitJobs(all_unfinished_jobs, resources=resources)
    # Wait for results.
    waitForJobs()
    message("Finished waiting for jobs at ", Sys.time())

    # Check status:
    print(getStatus())

    # Ensure that all jobs are done
    if (nrow(findNotDone()) != 0) {
      not_done_jobs <- findNotDone()
      print(getErrorMessages(not_done_jobs))
      ##browser()
      ##invisible(readline(prompt="Press [enter] to continue"))

      message("Some jobs did not complete. Please check the batchtools registry ", registry$file.dir)
      all_unfinished_jobs <- inner_join(not_done_jobs, all_unfinished_jobs)

      if (num_unfinished_jobs == nrow(all_unfinished_jobs) &&  nrow(all_unfinished_jobs) > 0.25 * num_all_jobs)
      {
        # Unfinished job count did not change -> retrying will probably not help. Abort!
        warning("Error: unfinished job count is not decreasing. Aborting job retries.")
        remaining_tries <- 0
      }

      if (num_unfinished_jobs == nrow(jobs))
      {
        # All jobs errored -> retrying will probably not help. Abort!
        warning("Error: all jobs errored. Aborting job retries.")
        remaining_tries <- 0
      }

      num_unfinished_jobs <- nrow(all_unfinished_jobs)
      message("Trying to resubmit jobs. Remaining tries: ", remaining_tries, " / ", max_num_tries)
    } else {
      all_jobs_finished <- TRUE
    }
  }

  invisible(all_jobs_finished)
}
