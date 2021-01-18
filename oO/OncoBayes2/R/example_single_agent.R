#' @title Single Agent Example
#'
#' @name example-single-agent
#'
#' @description
#' Example using a single experimental drug.
#'
#' @details
#' The single agent example is described in the reference
#' Neuenschwander, B. et al (2008). The data are described
#' in the help page for \code{hist_SA}. In this case, the data
#' come from only one study, with the treatment being only single
#' agent. Hence the model specified does not involve a hierarchical
#' prior for the intercept and log-slope parameters. The model
#' described in Neuenschwander, et al (2008) is adapted as follows:
#' \deqn{\mbox{logit}\, \pi(d) = \log\, \alpha + \beta \, \log\, \Bigl(\frac{d}{d^*}\Bigr),}
#' where \eqn{d^* = 250}, and the prior for
#' \eqn{\boldsymbol\theta = (\log\, \alpha, \log\, \beta)} is
#' \deqn{\boldsymbol\theta \sim \mbox{N}(\boldsymbol m, \boldsymbol S),}
#' and \eqn{\boldsymbol m = (\mbox{logit}\, 0.5, \log\, 1)} and
#' \eqn{\boldsymbol S = \mbox{diag}(2^2, 1^2)} are constants.
#'
#' In the \code{blrm_exnex} framework, in which
#' the prior must be specified as a hierarchical model
#' \eqn{\boldsymbol\theta \sim \mbox{N}(\boldsymbol \mu, \boldsymbol \Sigma)}
#' with additional priors on \eqn{\boldsymbol\mu} and \eqn{\boldsymbol\Sigma},
#' the simple prior distribution above is accomplished
#' by fixing the diagonal elements \eqn{\tau^2_\alpha} and \eqn{\tau^2_\beta}
#' of \eqn{\boldsymbol\Sigma} to zero, and taking
#' \deqn{\boldsymbol\mu \sim \mbox{N}(\boldsymbol m, \boldsymbol S).}
#'
#' The arguments \code{prior_tau_dist} and \code{prior_EX_tau_mean_comp}
#' as specified below ensure that the \eqn{\tau}'s are fixed at zero.
#'
#' @template ref-critical_aspects
#'
#' @template start-example
#' @template example-single_agent
#' @template stop-example
#'
NULL
