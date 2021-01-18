#' Posterior sampling from a hierarchical Unconstrained-Multinomial model
#'
#' MCMC sampling from a Dirichlet-Multinomial model using \code{\link[rstan]{stan}}.
#'
#' @param nsamples Number of posterior samples
#' @param X 4-column or 5-column matrix of observations in the correct format.  See \code{\link{UM.suff}}.
#' @param popId Optional vector of population identifiers.  See \code{\link{UM.suff}}.
#' @param rhoId Populations for which posterior samples of the genotype probability vector \code{rho} are desired.  Defaults to all populations.  Set \code{rhoId = NULL} not to output these for any populations.
#' @param full.stan.out Logical.  Whether or not to return the full \code{stan} output.  For monitoring convergence of the MCMC sampling.
#' @param ... Further arguments to be passed to the \code{\link[rstan]{sampling}} function in \pkg{rstan}.
#' @return A list with elements
#' \itemize{
#'   \item \code{A}: The unique allele names.
#'   \item \code{G}: The 4-column matrix Package libcurl was not found in the pkg-config search path.of unique genotype combinations.
#'   \item \code{rho}: A matrix with \code{ncol(rho) == nrow(G)}, where each row is a draw from the posterior distribution of inheritance probabilities.
#'   \item \code{sfit}: If \code{full.stan.out = TRUE}, the fitted \code{stan} object.
#' }
#' @details The hierarchical Dirichlet-Multinomial model is given by
#' \deqn{
#'   Y_k \mid \rho_k \sim_{\textrm{ind}} \textrm{Multinomial}(\rho_k, N_k),
#' }{
#'   Y_k | \rho_k ~ind Multinomial(\rho_k, N_k),
#' }
#' \deqn{
#'   \rho_k \sim_{\textrm{iid}} \textrm{Dirichlet}(\alpha).
#' }{
#'   \rho_k ~iid Dirichlet(\alpha).
#' }
#' where \eqn{\alpha_0 = \sum_{i=1}^C \alpha_i} and \eqn{\bar \alpha = \alpha/\alpha_0}{\alpha_bar = \alpha/\alpha_0}.  MCMC sampling is achieved with the \pkg{rstan} package, which is listed as a dependency for \pkg{MADPop} so as to expose \pkg{rstan}'s sophisticated tuning mechanism and convergence diagnostics.
#' @examples
#' # fit hierarchical model to fish215 data
#'
#' # only output posterior samples for lake Simcoe
#' rhoId <- "Simcoe"
#' nsamples <- 500
#' hUM.fit <- hUM.post(nsamples = nsamples, X = fish215,
#'                     rhoId = rhoId,
#'                     chains = 1) # number of MCMC chains
#'
#' # plot first 20 posterior probabilities in lake Simcoe
#' rho.post <- hUM.fit$rho[,1,]
#' boxplot(rho.post[,1:20], las = 2,
#'         xlab = "Genotype", ylab = "Posterior Probability",
#'         pch = ".", col = "grey")
#' @export
hUM.post <- function(nsamples, X, popId, rhoId,
                     full.stan.out = FALSE, ...) {
  suff <- UM.suff(X, popId)
  Xobs <- suff$tab
  if(missing(rhoId)) {
    rhoId <- rownames(Xobs)
  }
  if(is.null(rhoId)) {
    iLrho <- numeric(0)
  } else {
    iLrho <- sapply(rhoId, function(ri) {
      if(ri %in% rownames(Xobs)) {
        ans <- which(ri == rownames(Xobs))
      } else {
        ans <- 0
      }
      ans
    })
    if(any(iLrho == 0)) stop("rhoId must be a subset of popId.")
  }
  nLrho <- length(iLrho)
  sfit <- sampling(object = stanmodels$DirichletMultinomial,
                   data = list(nG = ncol(Xobs), nL = nrow(Xobs),
                     nLrho = nLrho, iLrho = array(iLrho, dim = nLrho),
                     X = Xobs),
                   iter = nsamples, ...)
  #alpha <- extract(sfit, permuted = TRUE)$alpha
  #colnames(alpha) <- colnames(Xobs)
  ans <- list(A = suff$A, G = suff$G)
  if(nLrho > 0) {
    rho <- extract(sfit, permuted = TRUE)$rho
    dimnames(rho)[[2]] <- rownames(Xobs)[iLrho]
    dimnames(rho)[[3]] <- colnames(Xobs)
    ans <- c(ans, list(rho = rho))
  }
  if(full.stan.out) ans <- c(ans, list(sfit = sfit))
  ans
}
