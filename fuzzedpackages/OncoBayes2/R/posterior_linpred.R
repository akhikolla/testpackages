#' Posterior of linear predictor
#'
#' Calculates the posterior of the linear predictor.
#'
#' @template args-methods
#' @template args-posterior
#' @template args-transform
#' @template args-dots-ignored
#'
#' @details
#'
#' Simulates the posterior of the linear predictor of the model
#' \code{object} for the specified data set.
#'
#' @template note-groups
#' @template return-samples
#'
#' @template start-example
#' @examples
#'
#' ## run single-agent analysis which defines blrmfit model object
#' example_model("single_agent")
#'
#' ## obtain posterior of linear prediction on 0-1 scale
#' post_prob_dlt <- posterior_linpred(blrmfit, TRUE, newdata=hist_SA)
#' ## name columns to obtain nice bayesplot labels
#' colnames(post_prob_dlt) <- hist_SA$drug_A
#'
#' library(bayesplot)
#' library(ggplot2)
#' mcmc_intervals(post_prob_dlt, prob=0.5, prob_outer=0.95) +
#'     coord_flip() +
#'     vline_at(c(0.16, 0.33), linetype=2) +
#'     ylab("Dose [mg]") +
#'     ggtitle("Posterior Probability of a DLT") +
#'     scale_x_continuous(breaks=c(0.1,0.16,0.33, 0.5, 0.75))
#'
#' @template stop-example
#'
#' @method posterior_linpred blrmfit
#' @aliases posterior_linpred
#' @export
posterior_linpred.blrmfit <- function(object, transform=FALSE, newdata, draws, ...) {
    dat <- pp_data(object, newdata=newdata, draws=draws)
    if (transform)
        dat <- inv_logit(dat)
    return(dat)
}

#' @method posterior_linpred blrm_trial
#' @export
posterior_linpred.blrm_trial <- function(object, transform=FALSE, newdata, draws, ...) {
    .assert_is_blrm_trial_and_prior_is_set(object)
    if(missing(newdata)) {
        return(posterior_linpred.blrmfit(object$blrmfit, transform=transform, newdata=object$data, draws=draws, ...))
    } else {
        return(posterior_linpred.blrmfit(object$blrmfit, transform=transform, newdata=newdata, draws=draws, ...))
    }
}
