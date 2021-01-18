#' Summarise model results
#'
#' Provides model summaries for \code{\link{blrm_exnex}} and
#' \code{\link{blrm_trial}} analyses.
#'
#' @template args-methods
#' @template args-prob
#' @param interval_prob optional vector of sorted quantiles for which
#'     the interval probabilities are calculated
#' @template args-transform
#' @param ... not used in this function
#'
#' @details
#' The calculated posterior summaries are returned as a
#' \code{data.frame} and contain optional interval probabilites for
#' the specified vector of sorted quantiles. These summaries are
#' calculated on the response scale by default and can be obtained on
#' the link scale when setting \code{transform=FALSE}.
#'
#' @return Returns a \code{data.frame} of the key summaries of the
#'     posterior mean, standard deviation, central probability
#'     interval, median and optional interval probabilities. Each row
#'     of the \code{data.frame} corresponds to the respective input
#'     data which is by default the same data set as used for the
#'     \code{\link{blrm_exnex}} analysis or the data specified in the
#'     \code{newdata} argument.
#'
#' @template start-example
#' @examples
#' example_model("single_agent")
#'
#' ## obtain underdosing (0-0.16), target dosing (0.16-0.33) and
#' ## overdosing (0.33-1) probabilities
#' summary(blrmfit, interval_prob=c(0,0.16,0.33,1))
#'
#' @template stop-example
#'
#' @method summary blrmfit
#' @export
summary.blrmfit <- function(object, newdata, transform=TRUE, prob=0.95, interval_prob, ...) {
    assert_numeric(prob, lower=0, upper=1, finite=TRUE, any.missing=FALSE, min.len=1)

    data <- object$data
    if(!missing(newdata))
        data <- newdata

    rn <- rownames(data)
    post_logit <- posterior_linpred(object, transform=transform, newdata=data)
    colnames(post_logit) <- 1:nrow(data)

    post_logit <- as_tibble(post_logit) %>%
        gather(id, value, convert=TRUE)

    ## dummy definitions to silence R CMD check
    . <- value <- Var1 <- Freq <- NULL

    probs <- sort(unique(c(0.5, 0.5- prob/2, 0.5+ prob/2)))
    qsummaries <- lapply(probs, function(p) function(x) quantile(x, probs=p)) %>%
      setNames(paste0(100 * probs, "%"))

    out_sum <- post_logit %>%
        group_by(id) %>%
        summarise_all(c(list(mean=mean, sd=sd), qsummaries)) %>%
        arrange(id) %>%
        ungroup()

    if(!missing(interval_prob)) {
        assert_numeric(interval_prob, any.missing=FALSE, sorted=TRUE)
        if(transform)
            assert_numeric(interval_prob, lower=0, upper=1, finite=TRUE)
        out_interval <- post_logit %>%
            group_by(id) %>%
            mutate(value=cut(value, breaks=c(-Inf, interval_prob, Inf))) %>%
            do(as.data.frame(prop.table(table(.$value)))) %>%
            spread(Var1, Freq)
        out_interval <- out_interval[-c(2,ncol(out_interval))]
        out_sum <- out_interval %>%
            right_join(y=out_sum, by="id") %>%
            select(unique(c(names(out_sum), names(out_interval))))
    }

    out_sum$id <- NULL
    out_sum <- as.data.frame(out_sum)
    rownames(out_sum) <- rn
    out_sum
}

#' Summarise trial
#'
#' Provides model summaries for \code{\link{blrm_trial}} analyses.
#' @param object \code{\link{blrm_trial}} object
#' @param summarize one of the following options:
#' \itemize{
#'   \item{}{\code{blrmfit}: summary of the underlying blrmfit object with further arguments ...}
#'   \item{}{\code{blrm_exnex_call}: blrm_exnex call used to create the \code{blrmfit} object}
#'   \item{}{\code{dose_info}: dose_info that were defined}
#'   \item{}{\code{dose_prediction} prediction for the defined dose_info}
#'   \item{}{\code{data}: data that were observed}
#'   \item{}{\code{data_prediction}: prediction for the observed data}
#'   \item{}{\code{dimensionality}: numeric vector with entries "num_components", "num_interaction_terms", "num_groups", "num_strata" }
#' }
#' @param ... further arguments for summary.blrmfit
#'
#' @template start-example
#' @examples
#' # construct initial blrm_trial object from built-in example datasets
#' combo2_trial_setup <- blrm_trial(
#'   data = hist_combo2,
#'   dose_info = dose_info_combo2,
#'   drug_info = drug_info_combo2,
#'   simplified_prior = TRUE
#' )
#'
#' # extract blrm_call to see setup of the prior as passed to blrm_exnex
#' summary(combo2_trial_setup, "blrm_exnex_call")
#'
#' @template stop-example
#'
#' @method summary blrm_trial
#' @export
summary.blrm_trial <- function (object, summarize=c("blrmfit", "blrm_exnex_call", "dose_info", "dose_prediction", "data", "data_prediction", "dimensionality"), ...) {
  summarize <- match.arg(summarize)

  if (summarize %in% c("blrmfit", "blrm_exnex_call", "dose_prediction", "data_prediction")) {
    .assert_is_blrm_trial_and_prior_is_set(object)
  } else {
    .assert_is_blrm_trial(object)
  }

  switch(summarize,
         blrmfit = summary(object$blrmfit, ...),
         blrm_exnex_call = object$blrmfit$call,
         dose_info = object$dose_info,
         dose_prediction = .blrm_trial_predict(object, object$dose_info, ...),
         data = object$data,
         data_prediction = .blrm_trial_predict(object, object$data, ...),
         dimensionality = object[c("num_components", "num_interaction_terms", "num_groups", "num_strata")])
}

