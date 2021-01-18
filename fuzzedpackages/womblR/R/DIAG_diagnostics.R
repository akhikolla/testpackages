###Function to get model fit diagnostics given a STBDwDM object
#'
#' diagnostics
#'
#' Calculates diagnostic metrics using output from the \code{\link{STBDwDM}} model.
#'
#' @param obj A \code{\link{STBDwDM}} model object for which diagnostics
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
#' @author Samuel I. Berchuck
#'
#' @references Gelfand, A. E., & Ghosh, S. K. (1998). Model choice: a minimum posterior predictive loss approach. Biometrika, 1-11.
#' @references Spiegelhalter, D. J., Best, N. G., Carlin, B. P., & Van Der Linde, A. (2002). Bayesian measures of model complexity and fit. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 64(4), 583-639.
#' @references Vehtari, A., Gelman, A., & Gabry, J. (2016). Practical Bayesian model evaluation using leave-one-out cross-validation and WAIC. Statistics and Computing, 1-20.
#' @references Watanabe, S. (2010). Asymptotic equivalence of Bayes cross validation and widely applicable information criterion in singular learning theory. Journal of Machine Learning Research, 11(Dec), 3571-3594.
#'
#' @export
diagnostics <- function(obj, diags = c("dic", "dinf", "waic"), keepDeviance = FALSE, keepPPD = FALSE) {

  ###Check Inputs
  if (missing(obj)) stop('"obj" is missing')
  if (!is.STBDwDM(obj)) stop('"obj" must be of class STBDwDM')
  if (sum((!diags %in% c("dic", "dinf", "waic"))) > 0) stop('"diags" must contain at least one of "dic", "dinf" or "waic"')
  if (!is.logical(keepDeviance)) stop('"keepDeviance" must be a logical')
  if (!is.logical(keepPPD)) stop('"keepPPD" must be a logical')

  ###Unload STBDwDM objects
  DatObj <- obj$datobj
  DatAug <- obj$dataug

  ###Check Inputs Again
  if ((DatObj$FamilyInd == 1) & (sum(diags %in% c("dic", "waic")) > 0))
    stop ('"probit" model cannot be used with "dic" or "waic"') #probit model can't do likelihood diagnostics

  ###Set seed for reproducibility
  set.seed(54)

  ###Set data objects
  M <- DatObj$M
  Z <- DatObj$Z
  AdjacentEdgesBoolean <- DatObj$AdjacentEdgesBoolean
  W <- DatObj$W
  EyeM <- DatObj$EyeM
  Rho <- DatObj$Rho
  FamilyInd <- DatObj$FamilyInd
  Nu <- DatObj$Nu
  YObserved <- DatObj$YObserved
  WeightsInd <- DatObj$WeightsInd

  ###Construct parameter object
  Para <- list()
  Para$Mu <- obj$mu
  Para$Tau2 <- obj$tau2
  Para$Alpha <- obj$alpha
  MuMean <- apply(obj$mu, 2, mean)
  Tau2Mean <- apply(obj$tau2, 2, mean)
  AlphaMean <- apply(obj$alpha, 2, mean)
  CovMean <- JointCovarianceCube(WAlphaCube(AlphaMean, Z, W, M, Nu, WeightsInd), Tau2Mean, EyeM, Rho, M, Nu)
  Para$MuMean <- MuMean
  Para$CovMean <- CovMean

  ###Set mcmc object
  NKeep <- dim(obj$phi)[1]

  ###Compute Log-likelihood using Rcpp function GetLogLik
  LogLik <- NULL
  if (("dic" %in% diags) | ("waic" %in% diags)) {

    ###Compute log-likelihood
    requireNamespace("mvtnorm", quietly = TRUE) #Requred for pmvnorm in Rcpp function
    if (DatObj$FamilyInd == 0) {
      NBelowCount <- c(0,0)
      YStarNonZero <- list()
      for (i in 1:DatObj$Nu) YStarNonZero[[i]] <- i
      DatAug$NBelowCount <- NBelowCount
      DatAug$YStarNonZero <- YStarNonZero
    }
    LogLik <- GetLogLik(DatObj, Para, DatAug, NKeep)

  }

  ###Compute DIC diagnostics
  dic <- NULL
  if ("dic" %in% diags) {

    ###Compute mean log-likelihood
    if (FamilyInd == 0) LogLikMean <- GetLogLikMean(DatObj, Para, DatAug)
    if (FamilyInd == 1) LogLikMean <- GetLogLikMean(DatObj, Para, DatAug)
    if (FamilyInd == 2) LogLikMean <- GetLogLikMean(DatObj, Para, DatAug)

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
    PPD <- SamplePPD(DatObj, Para, NKeep)

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
