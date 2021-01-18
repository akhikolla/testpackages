#' Return the number of posterior samples
#'
#' @param object fitted model object
#' @template args-dots-ignored
#'
#'
#' @template start-example
#' @examples
#'
#' ## run single-agent analysis which defines blrmfit model object
#' example_model("single_agent")
#'
#' nsamples(blrmfit)
#'
#' @template stop-example
#'
#' @method nsamples blrmfit
#' @aliases nsamples
#' @export
nsamples.blrmfit <- function(object, ...) {
    return(object$stanfit@sim$chains * (object$stanfit@sim$iter - object$stanfit@sim$warmup))
}

#' @method nsamples blrm_trial
#' @export
nsamples.blrm_trial <- function(object, ...) {
    .assert_is_blrm_trial_and_prior_is_set(object)
    return(nsamples.blrmfit(object$blrmfit))
}
