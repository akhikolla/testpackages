###Function to get model fit diagnostics given a spCP object
#'
#' diagnostics
#'
#' Calculates diagnostic metrics using output from the \code{\link{spCP}} model.
#'
#' @param obj A \code{\link{spCP}} model object for which diagnostics
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
  if (!is.spCP(obj)) stop('"obj" must be of class "spCP"')
  if (missing(obj)) stop('"obj" is missing')
  if (sum((!diags %in% c("dic", "dinf", "waic"))) > 0) stop('"diags" must contain at least one of "dic", "dinf" or "waic"')
  if (!is.logical(keepDeviance)) stop('"keepDeviance" must be a logical')
  if (!is.logical(keepPPD)) stop('"keepPPD" must be a logical')

  ###Unload STBDwDM objects
  DatObj <- obj$datobj
  DatAug <- obj$dataug

  ###Set seed for reproducibility
  set.seed(54)

  ###Set data objects
  OneM <- DatObj$OneM
  tNu <- DatObj$tNu
  t1 <- DatObj$t1
  XThetaInd <- DatObj$XThetaInd
  TimeVec <- DatObj$TimeVec
  OneNu <- DatObj$OneNu
  OneN <- DatObj$OneN
  N <- DatObj$N
  M <- DatObj$M
  YObserved <- DatObj$YObserved

  ###Construct parameter object
  Para <- list()
  Para$Beta0 <- obj$beta0
  Para$Beta1 <- obj$beta1
  Para$Lambda0 <- obj$lambda0
  Para$Lambda1 <- obj$lambda1
  Para$Eta <- obj$eta
  Beta0Mean <- apply(Para$Beta0, 2, median)
  Beta1Mean <- apply(Para$Beta1, 2, median)
  BetaMean <- matrix(rbind(Beta0Mean, Beta1Mean), ncol = 1)
  Lambda0Mean <- apply(Para$Lambda0, 2, median)
  Lambda1Mean <- apply(Para$Lambda1, 2, median)
  LambdaMean <- matrix(rbind(Lambda0Mean, Lambda1Mean), ncol = 1)
  EtaMean <- apply(Para$Eta, 2, median)
  ThetaMean <- pmax(pmin(tNu * OneM, EtaMean), t1 * OneM)
  XThetaMean = GetXTheta(ThetaMean, XThetaInd, TimeVec, OneNu, OneN, tNu, N, M)
  MuMean = XThetaMean %*% BetaMean;
  Sigma2Mean = exp(2 * (XThetaMean %*% LambdaMean));
  Para$MuMean <- MuMean
  Para$Sigma2Mean <- Sigma2Mean

  ###Set mcmc object
  NKeep <- dim(obj$delta)[1]

  ###Compute Log-likelihood using Rcpp function GetLogLik
  LogLik <- NULL
  if (("dic" %in% diags) | ("waic" %in% diags)) {

    ###Compute log-likelihood
    requireNamespace("mvtnorm", quietly = TRUE) #Requred for pmvnorm in Rcpp function
    LogLik <- GetLogLik(DatObj, Para, NKeep)

  }

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
