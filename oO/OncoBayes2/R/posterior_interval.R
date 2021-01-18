#' Posterior intervals
#'
#' Posterior intervals of all model parameters.
#'
#' @param object fitted model object
#' @template args-prob
#' @template args-dots-ignored
#'
#' @details
#'
#' Reports the quantiles of posterior parameters which correspond to
#' the central probability mass specified. The output includes the
#' posterior of the hyper-parameters and the posterior of each group
#' estimate.
#'
#' @return Matrix of two columns for the central probability interval
#'     \code{prob} for all parameters of the model.
#'
#' @template start-example
#' @examples
#' example_model("single_agent")
#'
#' posterior_interval(blrmfit)
#'
#' @template stop-example
#' @method posterior_interval blrmfit
#' @aliases posterior_interval
#' @export
posterior_interval.blrmfit <- function(object, prob=0.95, ...) {
    post <- as.matrix(object$stanfit, pars=c("mu_log_beta", "tau_log_beta", "rho_log_beta", "mu_eta", "tau_eta", "Sigma_corr_eta", "beta_group", "eta_group"))
    rstantools::posterior_interval(post, prob=prob)
}

#' @method posterior_interval blrm_trial
#' @export
posterior_interval.blrm_trial <- function(object, prob=0.95, ...) {
    .assert_is_blrm_trial_and_prior_is_set(object)
    return(posterior_interval.blrmfit(object$blrmfit, prob=prob, ...))
}
