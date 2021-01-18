#' Posterior of predictive
#'
#' Simulation of the predictive distribution.
#'
#' @template args-methods
#' @template args-posterior
#' @template args-dots-ignored
#'
#' @details
#'
#' Simulates the posterior predictive of the model \code{object} for
#' the specified data set.
#'
#' @template note-groups
#' @template return-samples
#'
#' @template start-example
#' @examples
#'
#' example_model("single_agent")
#'
#' post_pred  <- posterior_predict(blrmfit)
#' ## turn DLT counts into DLT rates
#' post_pred_rate <- sweep(post_pred, 2, hist_SA$num_patients, "/")
#'
#' library(bayesplot)
#' library(ggplot2)
#'
#' ## compare posterior predictive of the model for the response rates
#' ## with observed data
#' with(hist_SA,
#'     ppc_intervals(num_toxicities / num_patients, post_pred_rate, x=drug_A, prob_outer=0.95)) +
#'     xlab("Dose [mg]")
#'
#' @template stop-example
#'
#' @method posterior_predict blrmfit
#' @aliases posterior_predict
#' @export
posterior_predict.blrmfit <- function(object, newdata, draws, ...) {
    dat <- inv_logit(pp_data(object, newdata=newdata, draws=draws))
    num_sims <- nrow(dat)
    num_obs <- ncol(dat)
    num_trials <- pp_binomial_trials(object, newdata)

    ##pr <- t(apply(dat, 1, rbinom, n=num_obs, size=num_trials))
    pr <- matrix(rbinom(num_sims*num_obs, matrix(num_trials, nrow=num_sims, ncol=num_obs, byrow=TRUE), dat), nrow=num_sims, ncol=num_obs)
    colnames(pr) <- colnames(dat)
    pr
}

#' @method posterior_predict blrm_trial
#' @export
posterior_predict.blrm_trial <- function(object, newdata, draws, ...) {
    .assert_is_blrm_trial_and_prior_is_set(object)
    if(missing(newdata)) {
        return(posterior_predict.blrmfit(object$blrmfit, newdata=object$data, draws=draws, ...))
    } else {
        return(posterior_predict.blrmfit(object$blrmfit, newdata=newdata, draws=draws, ...))
    }
}
