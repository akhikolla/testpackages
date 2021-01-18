###################################
####    sample observations    ####
###################################

#----    sample_groups    ----

# Sample observations according to the required the sample size of the first
# group and second group (when needed), mean difference, and standard deviation
# ratio between the two groups.

sample_groups <- function(sample_n1, mean_diff, sample_n2=NULL, ratio_sd = 1){

  if(is.null(sample_n2)){
    res <- list(x = rnorm(sample_n1, mean = mean_diff, sd = 1),
                y = NULL)
  }else{
    res <- list(x = rnorm(sample_n1, mean = mean_diff, sd = ratio_sd),
                y = rnorm(sample_n2, mean = 0, sd = 1))
  }

  return(res)
}

#----    sample_obs_cor    ----

# Sample observations from a bivariate normal distribution according to the
# required correlation value (i.e., effect_target).

sample_obs_cor <- function(sample_n1, effect_target){

  obs <- mvrnorm(n=sample_n1,mu=c(0,0),
                 Sigma=matrix(c(1,effect_target,effect_target,1),ncol=2))

  return(list(x = obs[,1], y = obs[,2]))
}


#----    sample_effect    ----

# Sample effect size values from a function defined by the user and lower and
# upper truncation specification.

sample_effect <- function(FUN, B_effect, tl = -Inf, tu = Inf, tol = 1e4){
  if(!is.function(FUN) || length(formals(FUN))!=1L || !eval_rgn_function(FUN))
    stop(c("FUN has to be a function that allows to sample numeric values\n",
           "  The function has to be of the type 'function(n) my_function(n, ...)'\n",
           "  It requires only one single argument 'n' representing the number of sampled values\n",
           "  E.s. 'function(n) rnorm(n, mean = 0, sd = 1)'"))

  args <- list(x = B_effect)

  if(names(formals(FUN))!="x")
    names(args) <- names(formals(FUN))

  effect_function <- body(FUN)
  effect_samples <- do.call(FUN, args)

  # Truncate distribution
  if(is.finite(tl) || is.finite(tu)){
    message("Truncation could require long computational time")

    if(tl>tu) stop("Argument 'tl' has to be greater than argument 'tu'")

    # select out of bounds values
    sel_iter <- effect_samples < tl | effect_samples > tu
    i <- 1
    while(sum(sel_iter) != 0L && i < tol){
      args[[1]] <- sum(sel_iter)
      effect_samples[sel_iter] <- do.call(FUN, args)
      sel_iter <- effect_samples < tl | effect_samples > tu
      i <- i+1
    }

    if(i == tol) stop("Truncation requires too long computational time, consider possible misspecification.")
  }

  effect_summary <- summary(effect_samples)

  return(list(effect_function = effect_function,
              effect_summary = effect_summary,
              effect_samples = effect_samples))
}

#----

