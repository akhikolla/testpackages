
# from rstan
as.matrix.bayesGAMfit <- function(x, ...) {
  x <- as(x, "stanfit")
  if (x@mode != 0) return(numeric(0)) 
  e <- rstan::extract(x, permuted = FALSE, inc_warmup = FALSE, ...) 
  if (is.null(e)) return(e)
  enames <- dimnames(e)
  edim <- dim(e)
  dim(e) <- c(edim[1] * edim[2], edim[3])
  dimnames(e) <- enames[-2]
  e
}

# from rstan
as.array.bayesGAMfit <- function(x, ...) {
  x <- as(x, "stanfit")
  if (x@mode != 0) return(numeric(0)) 
  out <- rstan::extract(x, permuted = FALSE, inc_warmup = FALSE, ...)
  return(out)
}

# from rstan
as.data.frame.bayesGAMfit <- function(x, ...) {
  x <- as(x, "stanfit")
  as.data.frame(as.matrix(x, ...))
}

# from rstan
print.bayesGAMfit <- function(x, pars = x@sim$pars_oi, 
                              probs = c(0.025, 0.25, 0.5, 0.75, 0.975), 
                              digits_summary = 2, include = TRUE, ...) { 
  x <- as(x, "stanfit")
  if (x@mode == 1L) { 
    cat("Stan model '", x@model_name, "' is of mode 'test_grad';\n",
        "sampling is not conducted.\n", sep = '')
    return(invisible(NULL)) 
  } else if (x@mode == 2L) {
    cat("Stan model '", x@model_name, "' does not contain samples.\n", sep = '') 
    return(invisible(NULL)) 
  } 
  
  if(!include) pars <- setdiff(x@sim$pars_oi, pars)
  s <- summary(x, pars, probs, ...)  
  if (is.null(s)) return(invisible(NULL))
  n_kept <- x@sim$n_save - x@sim$warmup2
  cat("Inference for Stan model: ", x@model_name, '.\n', sep = '')
  cat(x@sim$chains, " chains, each with iter=", x@sim$iter, 
      "; warmup=", x@sim$warmup, "; thin=", x@sim$thin, "; \n", 
      "post-warmup draws per chain=", n_kept[1], ", ", 
      "total post-warmup draws=", sum(n_kept), ".\n\n", sep = '')
  
  # round n_eff to integers
  s$summary[, 'n_eff'] <- round(s$summary[, 'n_eff'], 0)
  
  print(round(s$summary, digits_summary), ...) 
  
  sampler <- attr(x@sim$samples[[1]], "args")$sampler_t
  
  if (!is.null(x@stan_args[[1]]$method) && 
      x@stan_args[[1]]$method == "variational") {
    if ("diagnostics" %in% names(x@sim)
        & "ir_idx" %in% names(x@sim$diagnostics)
        & !is.null(x@sim$diagnostics$ir_idx)) {
      cat("\nApproximate samples were drawn using VB(", x@stan_args[[1]]$algorithm, 
          ") + PSIS at ", x@date, ".\n", sep = '')
    } else {
      cat("\nApproximate samples were drawn using VB(", x@stan_args[[1]]$algorithm, 
          ") at ", x@date, ".\n", sep = '')
      message("We recommend genuine 'sampling' from the posterior distribution for final inferences!")
    }
    return(invisible(NULL))
  } else {
    cat("\nSamples were drawn using ", sampler, " at ", x@date, ".\n",
        "For each parameter, n_eff is a crude measure of effective sample size,\n", 
        "and Rhat is the potential scale reduction factor on split chains (at \n",
        "convergence, Rhat=1).\n", sep = '')
    return(invisible(NULL))
  }
}
