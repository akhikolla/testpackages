inference <- function(fidsamples, param, alpha = 0.05){
  fipred <- inherits(fidsamples, "gfilinreg.pred")
  out <- numeric(4L)
  names(out) <- c("mean", "median", "lwr", "upr")
  sample <-
    if(fipred) fidsamples[["FPD"]][[param]] else fidsamples[["Theta"]][[param]]
  weights <- fidsamples[["weight"]]
  out[1L] <- sum(sample * weights) # mean
  h <- cbind(sample, weights)
  hsort <- h[order(h[,1L]), ]
  hsum <- cumsum(hsort[, 2L])
  ci_u <- min(which(hsum >= 1-alpha/2))
  ci_l <- min(which(hsum >= alpha/2))
  ci_m <- min(which(hsum >= 0.5))
  out[3L] <- hsort[ci_l, 1L] # lower bound
  out[4L] <- hsort[ci_u, 1L] # upper bound
  out[2L] <- hsort[ci_m, 1L] # estimate (median)
  out
}

#' Summary of fiducial samples
#' @description Summary of the fiducial samples.
#'
#' @param fidsamples fiducial samples, the output of \code{\link{gfilinreg}} or 
#'   \code{\link{gfilinregPredictive}}
#' @param conf confidence level
#'
#' @return A matrix with summary statistics: means, medians, and confidence 
#'   intervals.
#' @export
#'
#' @examples set.seed(666L)
#' dat <- data.frame(
#'   group = gl(2, 15), 
#'   y = c(2*rlogis(15L), 10 + 2*rlogis(15L))
#' )
#' gfi <- gfilinreg(y ~ 0 + group, distr = "logistic", data = dat, L = 30L)
#' gfiSummary(gfi)
gfiSummary <- function(fidsamples, conf = 0.95){
  sims <- if(inherits(fidsamples, "gfilinreg.pred")){
    fidsamples[["FPD"]]
  }else fidsamples[["Theta"]]
  seq_ <- 1L:ncol(sims)
  names(seq_) <- names(sims)
  out <-
    t(vapply(seq_, function(x) inference(fidsamples, x, 1-conf), numeric(4L)))
  attr(out, "confidence level") <- conf
  out
}

#' Fiducial cumulative distribution function
#' @description Fiducial cumulative distribution function of a parameter of 
#'   interest.
#'
#' @param parameter a right-sided formula defining the parameter of interest 
#' @param fidsamples fiducial samples, the output of \code{\link{gfilinreg}} or 
#'   \code{\link{gfilinregPredictive}}
#'
#' @return The fiducial cumulative distribution function of the parameter.
#'
#' @importFrom lazyeval f_eval_rhs
#' @importFrom spatstat ewcdf
#' @export
#'
#' @examples set.seed(666L)
#' dat <- data.frame(
#'   group = gl(2, 15), 
#'   y = c(2*rlogis(15L), 10 + 2*rlogis(15L))
#' )
#' gfi <- gfilinreg(y ~ 0 + group, distr = "logistic", data = dat, L = 30L)
#' fcdf <- gfiCDF(~ group1 - group2, gfi)
#' fcdf(0)
#' plot(fcdf)
gfiCDF <- function(parameter, fidsamples){
  dataName <- ifelse(inherits(fidsamples, "gfilinreg.pred"), "FPD", "Theta")
  data <- fidsamples[[dataName]]
  fsims <- f_eval_rhs(parameter, data = data)
  ewcdf(fsims, weights = fidsamples[["weight"]])
}

#' Fiducial confidence interval
#' @description Fiducial confidence interval of a parameter of interest.
#'
#' @param parameter a right-sided formula defining the parameter of interest 
#' @param fidsamples fiducial samples, the output of \code{\link{gfilinreg}} or 
#'   \code{\link{gfilinregPredictive}}
#' @param conf confidence level
#'
#' @return The fiducial confidence interval of the parameter.
#'
#' @importFrom spatstat quantile.ewcdf
#' @export
#'
#' @examples set.seed(666L)
#' dat <- data.frame(
#'   group = gl(2, 15), 
#'   y = c(2*rlogis(15L), 10 + 2*rlogis(15L))
#' )
#' gfi <- gfilinreg(y ~ 0 + group, distr = "logistic", data = dat, L = 30L)
#' gfiConfInt(~ group1 - group2, gfi)
gfiConfInt <- function(parameter, fidsamples, conf = 0.95){
  fcdf <- gfiCDF(parameter, fidsamples)
  alpha <- 1 - conf
  quantile.ewcdf(fcdf, c(alpha/2, 1-alpha/2))
}

#' Fiducial quantiles
#' @description Quantiles of the fiducial distribution of a parameter of 
#'   interest.
#'
#' @param parameter a right-sided formula defining the parameter of interest 
#' @param fidsamples fiducial samples, the output of \code{\link{gfilinreg}} or 
#'   \code{\link{gfilinregPredictive}}
#' @param probs numeric vector of probabilities
#'
#' @return Numeric vector of quantiles, of the same length as \code{probs}.
#'
#' @importFrom spatstat quantile.ewcdf
#' @export
#'
#' @examples set.seed(666L)
#' dat <- data.frame(
#'   group = gl(2, 15), 
#'   y = c(2*rlogis(15L), 10 + 2*rlogis(15L))
#' )
#' gfi <- gfilinreg(y ~ 0 + group, distr = "logistic", data = dat, L = 30L)
#' gfiQuantile(~ group1 - group2, gfi, c(25, 50, 75)/100)
gfiQuantile <- function(parameter, fidsamples, probs){
  fcdf <- gfiCDF(parameter, fidsamples)
  quantile.ewcdf(fcdf, probs = probs)
}
