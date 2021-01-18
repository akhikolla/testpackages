#' emMixHMM implemens the EM (Baum-Welch) algorithm to fit a mixture of HMM
#' models.
#'
#' emMixHMM implements the maximum-likelihood parameter estimation of a mixture
#' of HMM models by the Expectation-Maximization (EM) algorithm, known as
#' Baum-Welch algorithm in the context of mixHMM.
#'
#' @details emMixHMM function implements the EM algorithm. This function starts
#'   with an initialization of the parameters done by the method `initParam` of
#'   the class [ParamMixHMM][ParamMixHMM], then it alternates between the E-Step
#'   (method of the class [StatMixHMM][StatMixHMM]) and the M-Step (method of
#'   the class [ParamMixHMM][ParamMixHMM]) until convergence (until the relative
#'   variation of log-likelihood between two steps of the EM algorithm is less
#'   than the `threshold` parameter).
#'
#' @param Y Matrix of size \eqn{(n, m)} representing the observed
#'   responses/outputs. `Y` consists of \emph{n} functions of `X` observed at
#'   points \eqn{1,\dots,m}.
#' @param K The number of clusters (Number of HMM models).
#' @param R The number of regimes (HMM components) for each cluster.
#' @param variance_type Optional character indicating if the model is
#'   "homoskedastic" or "heteroskedastic". By default the model is
#'   "heteroskedastic".
#' @param order_constraint Optional. A logical indicating whether or not a mask
#'   of order one should be applied to the transition matrix of the Markov chain
#'   to provide ordered states. For the purpose of segmentation, it must be set
#'   to `TRUE` (which is the default value).
#' @param init_kmeans Optional. A logical indicating whether or not the curve
#'   partition should be initialized by the K-means algorithm. Otherwise the
#'   curve partition is initialized randomly.
#' @param n_tries Optional. Number of runs of the EM algorithm. The solution
#'   providing the highest log-likelihood will be returned.
#'
#'   If `n_tries` > 1, then for the first run, parameters are initialized by
#'   uniformly segmenting the data into K segments, and for the next runs,
#'   parameters are initialized by randomly segmenting the data into K
#'   contiguous segments.
#' @param max_iter Optional. The maximum number of iterations for the EM
#'   algorithm.
#' @param threshold Optional. A numeric value specifying the threshold for the
#'   relative difference of log-likelihood between two steps of the EM as
#'   stopping criteria.
#' @param verbose Optional. A logical value indicating whether or not values of
#'   the log-likelihood should be printed during EM iterations.
#' @return EM returns an object of class [ModelMixHMM][ModelMixHMM].
#' @seealso [ModelMixHMM], [ParamMixHMM], [StatMixHMM]
#' @export
#'
#' @examples
#' data(toydataset)
#' Y <- t(toydataset[,2:ncol(toydataset)])
#'
#' mixhmm <- emMixHMM(Y = Y, K = 3, R = 3, verbose = TRUE)
#'
#' mixhmm$summary()
#'
#' mixhmm$plot()
emMixHMM <- function(Y, K, R, variance_type = c("heteroskedastic", "homoskedastic"), order_constraint = TRUE, init_kmeans = TRUE, n_tries = 1, max_iter = 1000, threshold = 1e-6, verbose = FALSE) {

  fData <- FData(X = seq.int(from = 0, to = 1, length.out = ncol(Y)), Y = Y)

  try_EM <- 0
  best_loglik <- -Inf

  while (try_EM < n_tries) {
    try_EM <- try_EM + 1

    if (n_tries > 1 && verbose) {
      message("EM try number: ", try_EM, "\n")
    }

    # Initialization
    variance_type <- match.arg(variance_type)
    param <- ParamMixHMM(fData = fData, K = K, R = R, variance_type = variance_type, order_constraint = order_constraint)

    param$initParam(init_kmeans, try_EM)

    iter <- 0
    converged <- FALSE
    prev_loglik <- -Inf

    stat <- StatMixHMM(paramMixHMM = param)

    # EM
    while ((iter <= max_iter) & !converged) {

      # E-Step
      stat$EStep(param)

      # M-Step
      param$MStep(stat)

      iter <- iter + 1

      if (verbose) {
        message("EM - mixHMMs: Iteration: ", iter, " | log-likelihood: "  , stat$loglik)
      }

      if (prev_loglik - stat$loglik > 1e-4) {
        warning("EM log-likelihood is decreasing from ", prev_loglik, "to ", stat$loglik, " !")
      }

      converged <- (abs((stat$loglik - prev_loglik) / prev_loglik) <= threshold)
      if (is.na(converged)) {
        converged <- FALSE
      } # Basically for the first iteration when prev_loglik is Inf

      prev_loglik <- stat$loglik
      stat$stored_loglik <- c(stat$stored_loglik, stat$loglik)

    } # End of EM loop

    if (stat$loglik > best_loglik) {
      statSolution <- stat$copy()
      paramSolution <- param$copy()

      best_loglik <- stat$loglik
    }

    if (n_tries > 1 && verbose) {
      message("Max value of the log-likelihood: ", stat$loglik, "\n\n")
    }

  }

  if (n_tries > 1 && verbose) {
    message("Best value of the log-likelihood: ", statSolution$loglik, "\n")
  }

  # Finding the curve partition by using the MAP rule
  statSolution$MAP()

  # Finish computation of statSolution
  statSolution$computeStats(paramSolution)

  return(ModelMixHMM(param = paramSolution, stat = statSolution))

}
