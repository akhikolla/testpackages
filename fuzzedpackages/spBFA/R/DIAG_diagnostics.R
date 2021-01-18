###Function to get model fit diagnostics given a spBFA object
#'
#' diagnostics
#'
#' Calculates diagnostic metrics using output from the \code{\link{spBFA}} model.
#'
#' @param object A \code{\link{spBFA}} model object for which diagnostics
#'  are desired from.
#'
#' @param diags A vector of character strings indicating the diagnostics to compute.
#'  Options include: Deviance Information Criterion ("dic"), d-infinity ("dinf") and
#'  Watanabe-Akaike information criterion ("waic"). At least one option must be included.
#'  Note: The probit model cannot compute the DIC or WAIC diagnostics due to computational
#'  issues with computing the multivariate normal CDF.
#'
#' @param keepDeviance A logical indicating whether the posterior deviance distribution
#'  is returned (default = FALSE).
#'
#' @param keepPPD A logical indicating whether the posterior predictive distribution
#'  at each observed location is returned (default = FALSE).
#'  
#' @param Verbose A boolean logical indicating whether progress should be output (default = TRUE).
#'
#' @param seed An integer value used to set the seed for the random number generator
#'  (default = 54).
#'  
#' @details To assess model fit, DIC, d-infinity and WAIC are used. DIC is based on the
#'  deviance statistic and penalizes for the complexity of a model with an effective
#'  number of parameters estimate pD (Spiegelhalter et al 2002). The d-infinity posterior
#'  predictive measure is an alternative diagnostic tool to DIC, where d-infinity=P+G.
#'  The G term decreases as goodness of fit increases, and P, the penalty term, inflates
#'  as the model becomes over-fit, so small values of both of these terms and, thus, small
#'  values of d-infinity are desirable (Gelfand and Ghosh 1998). WAIC is invariant to
#'  parametrization and is asymptotically equal to Bayesian cross-validation
#'  (Watanabe 2010). WAIC = -2 * (lppd - p_waic_2). Where lppd is the log pointwise
#'  predictive density and p_waic_2 is the estimated effective number of parameters
#'  based on the variance estimator from Vehtari et al. 2016. (p_waic_1 is the mean
#'  estimator).
#'
#' @return \code{diagnostics} returns a list containing the diagnostics requested and
#'  possibly the deviance and/or posterior predictive distribution objects.
#'
#' @examples
#' ###Load pre-computed regression results
#' data(reg.bfa_sp)
#' 
#' ###Compute and print diagnostics
#' diags <- diagnostics(reg.bfa_sp)
#' print(unlist(diags))
#' 
#' @author Samuel I. Berchuck
#'
#' @references Gelfand, A. E., & Ghosh, S. K. (1998). Model choice: a minimum posterior predictive loss approach. Biometrika, 1-11.
#' @references Spiegelhalter, D. J., Best, N. G., Carlin, B. P., & Van Der Linde, A. (2002). Bayesian measures of model complexity and fit. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 64(4), 583-639.
#' @references Vehtari, A., Gelman, A., & Gabry, J. (2016). Practical Bayesian model evaluation using leave-one-out cross-validation and WAIC. Statistics and Computing, 1-20.
#' @references Watanabe, S. (2010). Asymptotic equivalence of Bayes cross validation and widely applicable information criterion in singular learning theory. Journal of Machine Learning Research, 11(Dec), 3571-3594.
#'
#' @export
diagnostics <- function(object, diags = c("dic", "dinf", "waic"), keepDeviance = FALSE, keepPPD = FALSE, Verbose = TRUE, seed = 54) {

  ###Check Inputs
  if (missing(object)) stop('"object" is missing')
  if (!is.spBFA(object)) stop('"object" must be of class spBFA')
  if (sum((!diags %in% c("dic", "dinf", "waic"))) > 0) stop('"diags" must contain at least one of "dic", "dinf" or "waic"')
  if (!is.logical(keepDeviance)) stop('"keepDeviance" must be a logical')
  if (!is.logical(keepPPD)) stop('"keepPPD" must be a logical')
  if (!is.logical(Verbose)) stop('"Verbose" must be a logical')
  
  ###Unload spBFA objects
  DatObj <- object$datobj
  DatAug <- object$dataug

  ###Set seed for reproducibility
  set.seed(seed)

  ###Set mcmc object
  NKeep <- dim(object$psi)[1]
  
  ###Set data objects
  M <- DatObj$M
  O <- DatObj$O
  K <- DatObj$K
  EyeNu <- DatObj$EyeNu
  FamilyInd <- DatObj$FamilyInd
  Nu <- DatObj$Nu
  YObserved <- DatObj$YObserved
  X <- DatObj$X

  ###Construct parameter object
  Para <- list()
  Para$Lambda <- object$lambda
  Para$Eta <- object$eta
  if (DatObj$P > 0) Para$Beta <- object$beta
  if (DatObj$P == 0) Para$Beta <- matrix(0, ncol = 0, nrow = NKeep)
  if (is.null(object$sigma2)) Para$Sigma2 <- matrix(1)
  if (!is.null(object$sigma2)) Para$Sigma2 <- object$sigma2
  LambdaMean <- apply(object$lambda, 2, mean)
  EtaMean <- apply(object$eta, 2, mean)
  if (DatObj$P > 0) BetaMean <- apply(object$beta, 2, mean)
  if (DatObj$P == 0) BetaMean <- matrix(0, ncol = 1, nrow = 0)
  if (!is.null(object$sigma2)) Sigma2Mean <- apply(object$sigma2, 2, mean)
  if (is.null(object$sigma2)) Sigma2Mean <- matrix(1, nrow = 1, ncol = 1)
  Lambda <- matrix(LambdaMean, nrow = M * O, ncol = K, byrow = TRUE)
  Eta <- matrix(EtaMean, ncol = 1)
  MuMean <- array(kronecker(EyeNu, Lambda) %*% Eta + X %*% BetaMean, dim = c(M, O, Nu))
  Sigma2 <- t(matrix(Sigma2Mean, nrow = O, ncol = M))
  CovMean <- array(0, dim = c(M, O, Nu))
  count <- 1
  for (f in 1:O) {
    if (FamilyInd[f] %in% 0:2) {
      CovMean[ , f, ] <- do.call("cbind", rep(list(Sigma2[, count])), Nu)
      count <- count + 1
    }
  }
  Para$MuMean <- MuMean
  Para$CovMean <- CovMean

  ###Compute Log-likelihood using Rcpp function GetLogLik
  LogLik <- NULL
  if (("dic" %in% diags) | ("waic" %in% diags)) LogLik <- GetLogLik(DatObj, Para, NKeep, Verbose)

  ###Compute DIC diagnostics
  dic <- NULL
  if ("dic" %in% diags) {

    ###Compute mean log-likelihood
    LogLikMean <- GetLogLikMean(DatObj, Para)

    ###Calculate DIC objects
    DBar <- -2 * mean(LogLik)
    DHat <- -2 * LogLikMean
    pD <- DBar - DHat
    DIC <- DBar + pD
    dic <- list(dic = DIC, pd = pD)

  }

  ###Compute PPD diagnostics
  ppd <- PPD <- NULL
  if ("dinf" %in% diags) {

    ###Get PPD
    PPD <- SamplePPD(DatObj, Para, NKeep, Verbose)

    ###Compute PPD Diagnostics
    PPDMean <- apply(PPD, 1, mean)
    PPDVar <- apply(PPD, 1, var)
    P <- sum(PPDVar)
    G <- sum( (PPDMean - YObserved) ^ 2)
    DInf <- G + P
    ppd <- list(p = P, g = G, dinf = DInf)

  }

  ###Compute WAIC diagnostics
  waic <- NULL
  if ("waic" %in% diags) {

    ###Get WAIC
    # The calculation of Waic!  Returns lppd, p_waic_1, p_waic_2, and waic, which we define
    # as 2*(lppd - p_waic_2), as recommmended in BDA
    lppd <- log( apply(exp(LogLik), 2, mean) )
    if (!is.finite(lppd)) {
      M <- max(LogLik)
      lppd <- -log(dim(LogLik)[1]) + (M - log(sum(exp(LogLik - M))))
    }
    p_waic_1 <- 2 * (lppd - apply(LogLik, 2, mean) )
    p_waic_2 <- apply(LogLik, 2, var)
    waic <- -2 * lppd + 2 * p_waic_2
    waic <- list(waic = waic, p_waic = p_waic_2, lppd = lppd, p_waic_1 = p_waic_1)

  }

  ###Output diagnostics
  if (!keepDeviance & !keepPPD) diags <- list(dic = dic, dinf = ppd, waic = waic)
  if (!keepDeviance & keepPPD) diags <- list(dic = dic, dinf = ppd, waic = waic, PPD = t(PPD))
  if (keepDeviance & !keepPPD) diags <- list(dic = dic, dinf = ppd, waic = waic, deviance = -2 * LogLik)
  if (keepDeviance & keepPPD) diags <- list(dic = dic, dinf = ppd, waic = waic, deviance = -2 * LogLik, PPD = t(PPD))
  return(diags)

}
