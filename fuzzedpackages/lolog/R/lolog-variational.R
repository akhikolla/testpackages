


#' Fits a latent ordered network model using Monte Carlo variational inference
#'
#'
#' @param formula A lolog formula. See \code{link{lolog}}
#' @param nReplicates An integer controlling how many dyad ordering to perform.
#' @param dyadInclusionRate Controls what proportion of dyads in each ordering should be dropped.
#' @param targetFrameSize Sets dyadInclusionRate so that the model frame for the logistic regression will have on average this amount of observations.
#'
#'
#' @details
#' This function approximates the maximum likelihood solution via a variational inference on the
#' graph (y) over the latent edge variable inclusion order (s). Specifically, it replaces
#' the conditional probability p(s | y) by p(s). If the LOLOG model contains only dyad independent
#' terms, then these two probabilities are identical, and thus variational inference is
#' exactly maximum likelihood inference. The objective function is
#'
#' \deqn{E_{p(s)}\bigg(\log p(y| S, \theta) \bigg)}
#'
#' This can be approximated by drawing samples from p(s) to approximate the expectation. The
#' number of samples is controlled by the nReplicates parameter. The memory required is on the
#' order of nReplicates * (# of dyads). For large networks this can be impractical, so
#' adjusting dyadInclusionRate allows one to down sample the # of dyads in each replicate.
#'
#' If the model is dyad independent, replicates are redundant, and so nReplicates is set to
#' 1 with a note.
#'
#' The functional form of the objective function is equivalent to logistic regression, and so
#' the \code{\link{glm}} function is used to maximize it. The asymptotic covariance of the parameter
#' estimates is calculated using the methods of Westling (2015).
#'
#'
#' @return An object of class c('lologVariationalFit','lolog','list') consisting of the following
#' items:
#' \item{formula}{ The model formula}
#' \item{method}{"variational"}
#' \item{theta}{The fit parameter values}
#' \item{vcov}{The asymptotic covariance matrix for the parameter values.}
#' \item{nReplicates}{The number of replicates}
#' \item{dyadInclusionRate}{The rate at which dyads are included}
#' \item{allDyadIndependent}{Logical indicating model dyad independence}
#' \item{likelihoodModel}{An object of class *LatentOrderLikelihood at the fit parameters}
#' \item{outcome}{The outcome vector for the logistic regression}
#' \item{predictors}{The change statistic predictor matrix for the logistic regression}
#'
#'
#' @examples
#' library(network)
#' data(ukFaculty)
#' 
#' # Delete vertices missing group
#' delete.vertices(ukFaculty, which(is.na(ukFaculty %v% "Group")))
#' 
#' fit <- lologVariational(ukFaculty ~ edges() + nodeMatch("GroupC"),
#'                        nReplicates=1L, dyadInclusionRate=1)
#' summary(fit)
#'
#'
#' @references
#' Westling, T., & McCormick, T. H. (2015). Beyond prediction: A framework for inference with variational approximations in mixture models. arXiv preprint arXiv:1510.08151.
lologVariational <- function(formula,
                             nReplicates = 5L,
                             dyadInclusionRate = NULL,
                             targetFrameSize = 500000) {
  lolik <- createLatentOrderLikelihood(formula)
  nReplicates <- as.integer(nReplicates)
  
  dyadIndependent <- lolik$getModel()$isIndependent(TRUE, TRUE)
  dyadIndependentOffsets <-
    lolik$getModel()$isIndependent(TRUE, FALSE)
  allDyadIndependent <-
    all(dyadIndependent) & all(dyadIndependentOffsets)
  if (allDyadIndependent & nReplicates != 1L) {
    cat(
      "\n Model is dyad independent. Replications are redundant. Setting nReplicates <- 1L.\n"
    )
    nReplicates <- 1L
  }
  network <- lolik$getModel()$getNetwork()
  n <- network$size()
  ndyads <- n * (n - 1)
  if (!network$isDirected())
    ndyads <- ndyads / 2
  if (is.null(dyadInclusionRate)) {
    dyadInclusionRate <- min(1, targetFrameSize / ndyads)
  }
  samples <-
    lolik$variationalModelFrame(nReplicates, dyadInclusionRate)
  predictors <- lapply(samples, function(x)
    as.data.frame(x[[2]],
                  col.names = 1:length(x[[2]])))
  predictors <- do.call(rbind, predictors)
  outcome <- do.call(c, lapply(samples, function (x)
    x[[1]]))
  
  logFit <-
    glm(outcome ~ as.matrix(predictors) - 1, family = binomial())
  theta <- logFit$coefficients
  lolik$setThetas(theta)
  names(theta) <- names(lolik$getModel()$statistics())
  result <- list(
    method = "variational",
    formula = formula,
    theta = theta,
    vcov = vcov(logFit) * nReplicates / dyadInclusionRate,
    nReplicates = nReplicates,
    dyadInclusionRate = dyadInclusionRate,
    allDyadIndependent = allDyadIndependent,
    likelihoodModel = lolik,
    outcome = outcome,
    predictors = predictors
  )
  class(result) <- c("lologVariationalFit", "lolog", "list")
  result
}


#' Print of a lologVariationalFit object
#' @param x the object
#' @param ... additional parameters (unused)
#' @method print lologVariationalFit
print.lologVariationalFit <- function(x, ...) {
  if (x$allDyadIndependent)
    cat("MLE Coefficients:\n")
  else
    cat("Variational Inference Coefficients:\n")
  print(x$theta)
  if (x$dyadInclusionRate != 1) {
    cat("Inclusion rate:", x$dyadInclusionRate, "\n")
  }
  if (!x$allDyadIndependent)
    cat("# of replicates:", x$nReplicates, "\n")
}
