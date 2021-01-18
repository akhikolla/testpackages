#' PhenoFlex_GDHwrapper
#' 
#' PhenoFlex wrapper function for the `phenologyFitter` function using the GDH heat accumulation model 
#'
#' @param x data.frame with at least columns `Temp` and `JDay`
#' @param par numeric vector of length 11 with the `PhenoFlex` fit
#' parameters in the following order: 1. yc, 2. zc, 3. s1, 4. Tu, 5. E0,
#' 6. E1, 7, A0, 8. A1, 9. Tf, 10. Tc, 11. Tb and 12. slope. For details see
#' \link{PhenoFlex}
#'
#' @return
#' A single numeric value with the JDay prediction for the
#' temperaturs in `x$Temp` and the \link{PhenoFlex} parameters
#' in `par`.
#' 
#' @export PhenoFlex_GDHwrapper
PhenoFlex_GDHwrapper <- function(x, par) {
  ## for GDH Tb > Tu and Tu > Tc
  if(par[4] <= par[11]) return(NA)
  if(par[10] <= par[4]) return(NA)
  bloomindex <- PhenoFlex(temp=x$Temp,
                          times=seq_along(x$Temp),
                          yc=par[1],
                          zc=par[2],
                          s1=par[3],
                          Tu=par[4],
                          E0=par[5],
                          E1=par[6],
                          A0=par[7],
                          A1=par[8],
                          Tf=par[9],
                          Tc=par[10],
                          Tb=par[11],
                          slope=par[12],
                          Imodel=0L,
                          basic_output=TRUE)$bloomindex
  if(bloomindex == 0) return(NA)
  ## attempt to smooth the function a bit
  JDay <- x$JDay[bloomindex]
  JDaylist <- which(x$JDay == JDay)
  n <- length(JDaylist)
  if(n == 1) return(JDay)
  return(JDay + which(JDaylist == bloomindex)/n - 1./(n/ceiling(n/2)))
}

#' PhenoFlex_fixedDynModelwrapper
#' 
#' PhenoFlex wrapper function for the `phenologyFitter` function using
#' the GDH heat accumulation model and parameters of the dynamical model
#' fixed. The default values for the dynamic model parameters are from
#' the excel file with unknown origin.
#'
#' @param x data.frame with at least columns `Temp` and `JDay`
#' @param par numeric vector of length 11 with the `PhenoFlex` fit
#' parameters in the following order: 1. yc, 2. zc, 3. s1, 4. Tu, 5. E0,
#' 6. E1, 7, A0, 8. A1, 9. Tf, 10. Tc, 11. Tb and 12. slope. For details see
#' \link{PhenoFlex}
#' @param A0 numeric. Parameter \eqn{A_0}{A0} of the dynamic model
#' @param A1 numeric. Parameter \eqn{A_1}{A1} of the dynamic model
#' @param E0 numeric. Parameter \eqn{E_0}{E0} of the dynamic model
#' @param E1 numeric. Parameter \eqn{E_1}{E1} of the dynamic model
#' @param slope numeric. Slope parameter for sigmoidal function
#' @param Tf numeric. Transition temperature for the
#' sigmoidal function
#' @return
#' A single numeric value with the JDay prediction for the
#' temperaturs in `x$Temp` and the \link{PhenoFlex} parameters
#' in `par`.
#' 
#' @export
PhenoFlex_fixedDynModelwrapper <- function(x, par,
                                           A0=139500,
                                           A1=2567000000000000000,
                                           E0 = 4153.5,
                                           E1 = 12888.8,
                                           slope=1.6,
                                           Tf=4) {
  ## for GDH Tb > Tu and Tu > Tc
  if((par[4] <= par[6]) || (par[5] <= par[4])) return(NA)
  bloomindex <- PhenoFlex(temp=x$Temp,
                          times=seq_along(x$Temp),
                          yc=par[1],
                          zc=par[2],
                          s1=par[3],
                          Tu=par[4],
                          E0=E0,
                          E1=E1,
                          A0=A0,
                          A1=A1,
                          Tf=Tf,
                          Tc=par[5],
                          Tb=par[6],
                          slope=slope,
                          Imodel=0L,
                          basic_output=TRUE)$bloomindex
  if(bloomindex == 0) return(NA)
  ## attempt to smooth the function a bit
  JDay <- x$JDay[bloomindex]
  JDaylist <- which(x$JDay == JDay)
  n <- length(JDaylist)
  if(n == 1) return(JDay)
  return(JDay + which(JDaylist == bloomindex)/n - 1./(n/ceiling(n/2)))
}



#' PhenoFlex_GAUSSwrapper
#' 
#' PhenoFlex wrapper function for the `phenologyFitter` function using
#' the Gaussian heat accumulation model
#'
#' @param x data.frame with at least columns `Temp` and `JDay`
#' @param par numeric vector of length 11 with the `PhenoFlex` fit
#' parameters in the following order: 1. yc, 2. zc, 3. s1, 4. Tu, 5. E0,
#' 6. E1, 7, A0, 8. A1, 9. Tf, 10. Delta, 11. s. For details see
#' \link{PhenoFlex}
#'
#' @return
#' A single numeric value with the JDay prediction for the
#' temperatures in `x$Temp` and the \link{PhenoFlex} parameters
#' in `par`.
#' 
#' @export
PhenoFlex_GAUSSwrapper <- function(x, par) {
  bloomindex <- PhenoFlex(temp=x$Temp,
                          times=seq_along(x$Temp),
                          yc=par[1],
                          zc=par[2],
                          s1=par[3],
                          Tu=par[4],
                          E0=par[5],
                          E1=par[6],
                          A0=par[7],
                          A1=par[8],
                          Tf=par[9],
                          Delta=par[10],
                          slope=par[11],
                          Imodel=1L,
                          basic_output=TRUE)$bloomindex
  if(bloomindex == 0) return(NA)
  ## attempt to smooth the function a bit
  JDay <- x$JDay[bloomindex]
  JDaylist <- which(x$JDay == JDay)
  n <- length(JDaylist)
  if(n == 1) return(JDay)
  return(JDay + which(JDaylist == bloomindex)/n - 1./(n/ceiling(n/2)))
}

#' PhenoFlex_fixedDynModelGAUSSwrapper
#' 
#' PhenoFlex wrapper function for the `phenologyFitter` function using
#' the GAUSS heat accumulation model and parameters of the dynamical model
#' fixed.
#'
#' @param x data.frame with at least columns `Temp` and `JDay`
#' @param par numeric vector of length 11 with the `PhenoFlex` fit
#' parameters in the following order: 1. yc, 2. zc, 3. s1, 4. Tu, 5. E0,
#' 6. E1, 7, A0, 8. A1, 9. Tf, 10. Tc, 11. Tb and 12. slope. For details see
#' \link{PhenoFlex}
#' @param A0 numeric. Parameter \eqn{A_0}{A0} of the dynamic model
#' @param A1 numeric. Parameter \eqn{A_1}{A1} of the dynamic model
#' @param E0 numeric. Parameter \eqn{E_0}{E0} of the dynamic model
#' @param E1 numeric. Parameter \eqn{E_1}{E1} of the dynamic model
#' @param slope numeric. Slope parameter for sigmoidal function
#' @param Tf numeric. Transition temperature (in degree Kelvin) for the
#' sigmoidal function
#' @return
#' A single numeric value with the JDay prediction for the
#' temperaturs in `x$Temp` and the \link{PhenoFlex} parameters
#' in `par`.
#' 
#' @export
PhenoFlex_fixedDynModelGAUSSwrapper <- function(x, par,
                                                A0=139500,
                                                A1=2567000000000000000,
                                                E0 = 4153.5,
                                                E1 = 12888.8,
                                                slope=1.6,
                                                Tf=4) {
  bloomindex <- PhenoFlex(temp=x$Temp,
                          times=seq_along(x$Temp),
                          yc=par[1],
                          zc=par[2],
                          s1=par[3],
                          Tu=par[4],
                          E0=E0,
                          E1=E1,
                          A0=A0,
                          A1=A1,
                          Tf=Tf,
                          Delta=par[5],
                          slope=slope,
                          Imodel=1L,
                          basic_output=TRUE)$bloomindex
  if(bloomindex == 0) return(NA)
  ## attempt to smooth the function a bit
  JDay <- x$JDay[bloomindex]
  JDaylist <- which(x$JDay == JDay)
  n <- length(JDaylist)
  if(n == 1) return(JDay)
  return(JDay + which(JDaylist == bloomindex)/n - 1./(n/ceiling(n/2)))
}

predictBloomDays <- function(par, SeasonList, modelfn, ...) {
  paravail <- requireNamespace("parallel")
  sres <- 
    if(FALSE) { ## currently, parallel seems to be slower
      mc.cores <- as.numeric(Sys.getenv("OMP_NUM_THREADS"))
      if(is.na(mc.cores)) mc.cores <- getOption("mc.cores", default=1L)
      parallel::mclapply(X=SeasonList, FUN=modelfn, par=par, ...,
                         mc.cores = mc.cores)
    }
    else {
      lapply(X=SeasonList, FUN=modelfn, par=par, ...)
    }
  return(invisible(simplify2array(x=sres)))
}

#' chifull
#'
#' function to compute the RSS
#'
#' @description
#' RSS to minimise by `phenologyFitter`
#' @param par numeric. vector of fit parameters
#' @param modelfn function. model function
#' @param bloomJDays numeric. vector of bloom hours! per year
#' @param SeasonList list. list of index vectors per year.
#' @param na_penalty numeric. penalty value for the residual if the model returns `NA`.
#' @param ... further parameters to pass on to `modelfn`.
chifull <- function(par,
                    modelfn,
                    bloomJDays,
                    SeasonList,
                    na_penalty=365,
                    ...) {

  sres <- predictBloomDays(par=par, SeasonList=SeasonList,
                              modelfn=modelfn, ...)
  s <- (sres-bloomJDays)
  ## we penalise if no bloomday could be found
  nai <- which(is.na(sres))
  s[nai] <- na_penalty
  ## Remaining NAs can only come from the data, not from the predictions
  return(sum(s[!is.na(s)]^2))
}


#' phenologyFit
#'
#' Constructor for class `phenologyFit`
#' 
#' @return
#' an empty object of class `phenologyFit`. It contains the named elements `model_fit` with 
#' the returned object from GenSA, `par` the best fit parameters, `pbloomJDays` the
#' predicted bloom JDays and the inputs
#' `par.guess`, `modelfn`, `bloomJDays`, and `SeasonList`. They are
#' all set to `NULL` by this function.
#' 
#' @export
phenologyFit <- function() {
  res <- list()
  attr(res, "class") <- c("phenologyFit", "list")
  res$model_fit <- NULL
  res$par <- NULL
  res$pbloomJDays <- NULL
  res$bloomJDays <- NULL
  res$par.guess <- NULL
  res$SeasonList <- NULL
  return(res)
}

#' @title phenologyFitter
#'
#' @param par.guess numeric vector. Initial guesses for fit parameters. This can be set to
#' `NULL`, in which case `GenSA` choses initial parameters.
#' @param modelfn function. Model function which computes the index in `temperatures` at
#' which blooming occures. It must have as first argument a data frame with at least the
#' two columns `Temp` and `JDays` for one season, see `SeasonList`. It can have further
#' arguments which can be passed via `...`. The `modelfn` must return a single numeric value
#' for the predicted bloom JDay for that season. `NA` is an allowed return value if no blooming
#' occures in that season.
#' The default is the \link{PhenoFlex} with GDH as heat accumulation. Alternative is
#' \link{PhenoFlex_GAUSSwrapper} with GAUSSian heat accumulation. But this function can also be user
#' defined. 
#' @param bloomJDays integer vector. vector of observed bloom JDays per year
#' @param SeasonList list. Must be a list of data frames, each data frame for one season. Each data.frame
#' must at least have a column `Temp` with the temperature vector and `JDays` with the corresponding
#' JDay vector. Can be generated by e.g. \link{genSeasonList}.
#' `length(SeasonList)` must be equal to `length(bloomJDays)`.
#' @param control control parameters to `GenSA`, see `GenSA::GenSA`
#' @param lower Vector with length of ‘par.guess’. Lower bounds for components.
#' @param upper Vector with length of ‘par.guess’. Upper bounds for components.
#' @param seed integer seed for the random number generator used by `GenSA`.
#' @param ... further parameters to be passed on to `modelfn`.
#'
#' @importFrom GenSA GenSA
#' @author Carsten Urbach <urbach@hiskp.uni-bonn.de>
#' 
#' @return
#' an object of class `phenologyFit`. It contains the named elements `model_fit` with
#' the returned object from GenSA, `par` the best fit parameters, `pbloomJDays` the
#' predicted bloom JDays and the inputs
#' `par.guess`, `modelfn`, `bloomJDays`, `lower`, `upper`, `control`, `SeasonList` and
#' `...`.
#'
#' @examples
#' ## this example does not make sense as a fit, but demonstrates
#' ## how to use `phenologyFitter`
#' data(KA_weather)
#' data(KA_bloom)
#' hourtemps <- stack_hourly_temps(KA_weather, latitude=50.4)
#' SeasonList <- genSeasonList(hourtemps$hourtemps, years=c(2007,2008))
#' par <- c(40, 190, 0.5, 25, 3372.8, 9900.3, 6319.5, 5.939917e13, 4, 36, 4, 1.6)
#' upper <- c(41, 200, 1, 30, 4000, 10000, 7000, 6.e13, 10, 40, 10, 50)
#' lower <- c(38, 180, 0.1, 0, 3000, 9000, 6000, 5.e13, 0, 0, 0, 0.05)
#' X <- phenologyFitter(par.guess=par, bloomJDays=KA_bloom$pheno[c(24,25)], 
#'   SeasonList=SeasonList, lower=lower, upper=upper,
#'   control=list(smooth=FALSE, verbose=TRUE, maxit=10, nb.stop.improvement=5))
#' summary(X)
#' plot(X)
#' @export
phenologyFitter <- function(par.guess=NULL,
                            modelfn=PhenoFlex_GDHwrapper,
                            bloomJDays,
                            SeasonList,
                            control=list(smooth=FALSE, verbose=TRUE, maxit=1000,
                                         nb.stop.improvement=250),
                            lower, upper,
                            seed = 1235433,
                            ...) {

  control$seed <- seed

  ## some sanity checks
  stopifnot(is.list(SeasonList))
  stopifnot(length(SeasonList) == length(bloomJDays))

  ## generate empty `phenologyFit` and store input objects
  res <- phenologyFit()
  res$par.guess <- par.guess
  res$modelfn  <- modelfn
  res$bloomJDays <- bloomJDays
  res$SeasonList <- SeasonList
  res$lower <- lower
  res$upper <- upper
  res$control <- control
  ##res$... <- ...
  
  ## cure the 365 to 1 jump at new years in JDays and, possibly, bloomJDays
  for(i in c(1:length(SeasonList))) {
    minJDay <- SeasonList[[i]]$JDay[1]
    maxJDay <- SeasonList[[i]]$JDay[length(SeasonList[[i]]$JDay)]
    ## some more sanity checks
    if(maxJDay > minJDay) {
      stop(paste0("Season ", i, " is overlapping with the previous or following one. Aborting!"))
    }
    if(bloomJDays[i] > maxJDay && bloomJDays[i] < minJDay) {
      stop(paste0("In season ", i, " the bloomJDay is outside the provided JDay vector. Aborting!"))
    }
    ## determine the index of the jump and its value
    dx <- diff(SeasonList[[i]]$JDay)
    mx <- min(dx)
    kmx <- which(dx == mx)
    SeasonList[[i]]$JDay[1:kmx] <- SeasonList[[i]]$JDay[1:kmx] + mx - 1
    if(bloomJDays[i] > minJDay) bloomJDays[i]  <- bloomJDays[i] + mx - 1
    ## store this also in res for later usage
    res$SeasonList[[i]]$JDayunwrapped <- SeasonList[[i]]$JDay
    res$bloomJDaysunwrapped <- bloomJDays
  }

  ## now we fit
  res$model_fit <- GenSA::GenSA(par=par.guess, fn=chifull, bloomJDays=bloomJDays,
                                SeasonList=SeasonList, 
                                modelfn=modelfn,
                                control=control,
                                lower=lower, upper=upper)
  res$par <- res$model_fit$par
  res$pbloomJDays <- predictBloomDays(par=res$par,
                                      SeasonList=res$SeasonList,
                                      modelfn=res$modelfn)
  return(res)
}

#' @title summary phenologyFit
#' 
#' @param object class phenologyFit. object to summarise
#' @param ... additional parameters, ignored here
#' 
#' @return
#' No return value.
#'
#' @export
summary.phenologyFit <- function(object, ...) {
  ## compute R^2
  ybar <- mean(object$bloomJDays)
  sqe <- sum((object$pbloomJDays - ybar)^2)
  sqt <- sum((object$bloomJDays - ybar)^2)
  cat("\n")
  cat("Phenology Fitter\n\n")
  cat("Final RSS: ", object$model_fit$value, "\n")
  cat("RMSE     : ", sqrt(object$model_fit$value/length(object$pbloomJDays)), "\n")
  cat("R^2      : ", sqe/sqt, "\n\n")
  cat("data versus predicted:\n")
  print(data.frame(data=object$bloomJDays, predicted=object$pbloomJDays, delta=object$bloomJDays-object$pbloomJDays))
}


#' @title print phenologyFit
#' 
#' @param x class phenologyFit. object to print
#' @param ... additional parameters, ignored here
#' 
#' @return
#' No return value.
#'
#' @export
print.phenologyFit <- function(x, ...) {
  summary(x)
}

#' predict phenologyFit
#'
#' Generic function to predict a `phenologyFit` object.
#' 
#' @param object object of class `phenologyFit` to predict.
#' @param SeasonList List with data frames per season, see
#' \link{phenologyFit} for more details.
#' @param ... additional parameters, ignored here
#'
#' @return
#' A numeric vector is returned with a predicted bloom day per
#' Season in `SeasonList`. If `SeasonList` is missing, the
#' original `SeasonList` is used for prediction.
#' 
#' @export
predict.phenologyFit <- function(object,
                                 SeasonList,
                                 ...) {
  x <- 
  if(missing(SeasonList)) {
    object$pbloomJDays
  }
  else {
    predictBloomDays(par=object$par,
                     SeasonList=SeasonList,
                     modelfn=object$modelfn)
  }
  return(x)
}

#' plot phenologyFit
#'
#' Generic function to plot a `phenologyFit` object
#'
#' @param x object of class `phenologyFit` to plot.
#' @param ylim numeric vector of length 2 with the limit for the y-axis
#' @param ... additional graphical parameters to pass on.
#'
#' @importFrom graphics legend
#' @return
#' No return value.
#' 
#' @export
plot.phenologyFit <- function(x,
                              ylim=c(0.9*min(c(x$bloomJDays, x$pbloomJDays)), 1.1*max(c(x$bloomJDays, x$pbloomJDays))),
                              ...) {
  plot(x=seq_along(x$SeasonList), y=x$bloomJDays,
       xlab="Season", ylab="JDay",
       col="blue", pch=21, ylim=ylim,
       ...)
  points(x=seq_along(x$SeasonList), y=x$pbloomJDays,
         col="red", pch=22)
  legend("topleft",
         legend=c("data", "predicted"),
         bty="n", pch=c(21,22), col=c("blue", "red"))
}


bootstrap_phenologyFit <- function() {
  res <- list()
  attr(res, "class") <- c("bootstrap_phenologyFit", "list")
  return(res)
}

