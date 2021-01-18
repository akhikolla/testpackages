#' emMixHMMR implements the EM algorithm to fit a mixture if HMMR models.
#'
#' emMixHMMR implements the maximum-likelihood parameter estimation of a mixture
#' of HMMR models by the Expectation-Maximization (EM) algorithm.
#'
#' @details emMixHMMR function implements the EM algorithm. This function starts
#'   with an initialization of the parameters done by the method `initParam` of
#'   the class [ParamMixHMMR][ParamMixHMMR], then it alternates between the
#'   E-Step (method of the class [StatMixHMMR][StatMixHMMR]) and the M-Step
#'   (method of the class [ParamMixHMMR][ParamMixHMMR]) until convergence (until
#'   the relative variation of log-likelihood between two steps of the EM
#'   algorithm is less than the `threshold` parameter).
#'
#' @param X Numeric vector of length \emph{m} representing the covariates/inputs
#'   \eqn{x_{1},\dots,x_{m}}.
#' @param Y Matrix of size \eqn{(n, m)} representing the observed
#'   responses/outputs. `Y` consists of \emph{n} functions of `X` observed at
#'   points \eqn{1,\dots,m}.
#' @param K The number of clusters (Number of HMMR models).
#' @param R The number of regimes (HMMR components) for each cluster.
#' @param p Optional. The order of the polynomial regression. By default, `p` is
#'   set at 3.
#' @param variance_type Optional. character indicating if the model is
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
#' @return EM returns an object of class [ModelMixHMMR][ModelMixHMMR].
#' @seealso [ModelMixHMMR], [ParamMixHMMR], [StatMixHMMR]
#' @export
#'
#' @examples
#' data(toydataset)
#' x <- toydataset$x
#' Y <- t(toydataset[,2:ncol(toydataset)])
#'
#' mixhmmr <- emMixHMMR(X = x, Y = Y, K = 3, R = 3, p = 1, verbose = TRUE)
#'
#' mixhmmr$summary()
#'
#' mixhmmr$plot()
emMixHMMR <- function(X, Y, K, R, p = 3, variance_type = c("heteroskedastic", "homoskedastic"), order_constraint = TRUE, init_kmeans = TRUE, n_tries = 1, max_iter = 1000, threshold = 1e-6, verbose = FALSE) {

  fData <- FData(X = X, Y = Y)

  try_EM <- 0
  best_loglik <- -Inf

  while (try_EM < n_tries) {
    try_EM <- try_EM + 1
    if (n_tries > 1 && verbose) {
      message("EM try number: ", try_EM, "\n")
    }

    # Initialization
    variance_type <- match.arg(variance_type)
    param <- ParamMixHMMR(fData = fData, K = K, R = R, p = p, variance_type = variance_type, order_constraint = order_constraint)

    param$initParam(init_kmeans, try_EM)

    iter <- 0
    converged <- FALSE
    prev_loglik <- -Inf

    stat <- StatMixHMMR(paramMixHMMR = param)

    # EM
    while ((iter <= max_iter) & !converged) {

      # E-Step
      stat$EStep(param)

      # M-Step
      param$MStep(stat)

      iter <- iter + 1

      if (verbose) {
        message("EM - mixHMMR: Iteration: ", iter, " || log-likelihood: "  , stat$loglik)
      }

      if (prev_loglik - stat$loglik > 1e-4) {
        warning("EM log-likelihood is decreasing from ", prev_loglik, "to ", stat$loglik, " !")
      }

      converged <- (abs((stat$loglik - prev_loglik) / prev_loglik) < threshold)
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
  statSolution$computeStats(param)

  return(ModelMixHMMR(param = paramSolution, stat = statSolution))

}
