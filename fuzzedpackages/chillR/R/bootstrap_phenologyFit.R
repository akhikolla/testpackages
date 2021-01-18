#' bootstrap.phenologyFit
#'
#' bootstrap an object of S3 class `phenologyFit`
#'
#' @description
#' This function bootstraps the residuals of a `phenologyFit`. It
#' internally calls `phenologyFitter` on each bootstrap
#' replicate.
#' 
#' @param object class `phenologyFit`, the object to bootstrap
#' @param boot.R integer. The number of bootstrap replicates
#' @param control control parameters to `GenSA`, see `GenSA::GenSA`
#' @param lower Vector with length of ‘par.guess’. Lower bounds for components.
#' @param upper Vector with length of ‘par.guess’. Upper bounds for components.
#' If missing, `upper` in `object` is used.
#' @param seed integer seed for the random number generator used by `GenSA`.
#' If missing, `lower` in `object` is used.
#'
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom stats predict
#' 
#' @author Carsten Urbach <urbach@hiskp.uni-bonn.de>
#' @return
#' Invisibly returns a list with elements `boot.R`, `object`, `seed`, `residuals`,
#' `lower`, `upper`, and `res`. The latter list `res` has `boot.R` elements, which
#' are lists again. Each of these lists contains named elements `par`, `value`,
#' `bloomJDays`, and `pbloomJDays`. `par` are the best fit parameters on the particular bootstrap 
#' replicate, `value` the corresponding RSS, `bloomJDays` the re-sampled data and `pbloomJDays`
#' the predicted bloom JDays for this sample.
#' 
#' @export
bootstrap.phenologyFit <- function(object,
                                   boot.R=99,
                                   control=list(smooth=FALSE, verbose=FALSE, maxit=1000,
                                                nb.stop.improvement=250),
                                   lower, upper,
                                   seed=1766588
                                   ) {
  stopifnot(!missing(object))
  stopifnot(boot.R > 1)
  boot.R <- as.integer(boot.R)
  if(missing(lower)) lower <- object$lower
  if(missing(upper)) upper <- object$upper
  ## compute residuals
  residuals <- object$bloomJDays - object$pbloomJDays
  if(all(residuals == 0)) {
    stop("All residuals equal to zero, no variation. Aborting")
  }
  ## container for the results
  bootres <- list(res=list(), boot.R=boot.R, object=object,
                  seed=seed, lower=lower, upper=upper)
  attr(bootres, "class") <- c("bootstrap_phenologyFit", "list")
  ## loop over replicates
  pb <- utils::txtProgressBar(min=0, max=boot.R, initial=0)
  for(r in c(1:boot.R)) {
    bloomJDays <- object$pbloomJDays + sample(residuals, length(residuals), replace=TRUE)
    tmp <- phenologyFitter(par.guess=object$par,
                           modelfn=object$modelfn,
                           bloomJDays=bloomJDays,
                           SeasonList=object$SeasonList,
                           control=control,
                           lower=lower,
                           upper=upper,
                           seed=seed)
    bootres$res[[r]] <- list(par=tmp$par, value=tmp$model_fit$value,
                             bloomJDays=bloomJDays, pbloomJDays=tmp$pbloomJDays)
    utils::setTxtProgressBar(pb, r)
  }
  close(pb)
  return(invisible(bootres))
}

#' summary.bootstrap_phenologyFit
#'
#' Summarise a `bootstrap_phenologyFit` object
#' 
#' @param object class `bootstrap_phenologyFit` to summarise
#' @param ... generic parameters, ignored here
#'
#' @return
#' No return value.
#' 
#' @export
summary.bootstrap_phenologyFit <- function(object, ...) {
  Y <- t(array(unlist(lapply(object$res, FUN=function(x) return(x$par))), dim=c(length(object$object$par), sum(object$boot.R))))
  Err <- apply(Y, FUN=sd, MARGIN=2)
  q16 <- apply(Y, FUN=quantile, MARGIN=2, probs=0.16)
  q84 <- apply(Y, FUN=quantile, MARGIN=2, probs=0.84)
  ans <- data.frame(par=object$object$par, Err=Err, q16=q16, q84=q84)
  ans
}

#' plot bootstrap_phenologyFit
#'
#' Generic function to plot a `bootstrap_phenologyFit` object
#'
#' @param x object of class `bootstrap_phenologyFit` to plot.
#' @param ylim numeric vector of length 2 with the limit for the y-axis
#' @param ... additional graphical parameters to pass on.
#'
#' @importFrom graphics legend arrows
#' @return
#' No return value.
#' 
#' @export
plot.bootstrap_phenologyFit <- function(x,
                                        ylim=c(0.9*min(c(x$object$bloomJDays, x$object$pbloomJDays)), 1.1*max(c(x$object$bloomJDays, x$object$pbloomJDays))),
                                        ...) {
  plot(x=seq_along(x$object$SeasonList), y=x$object$bloomJDays,
       xlab="Season", ylab="JDay",
       col="blue", pch=21, ylim=ylim,
       ...)
  points(x=seq_along(x$object$SeasonList), y=x$object$pbloomJDays,
         col="red", pch=22)
  Y <- t(array(unlist(lapply(x$res, FUN=function(x) return(x$pbloomJDays))), dim=c(length(x$object$pbloomJDays), sum(x$boot.R))))
  Err <- apply(Y, FUN=sd, MARGIN=2)
  arrows(x0=seq_along(x$object$SeasonList), y0=x$object$pbloomJDays-Err,
         x1=seq_along(x$object$SeasonList), y1=x$object$pbloomJDays+Err,
         length=0, code=1, angle=90, col="red")
  legend("topleft",
         legend=c("data", "predicted"),
         bty="n", pch=c(21,22), lty=c(NA, 1), col=c("blue", "red"))
}

#' predict bootstrap_phenologyFit
#'
#' Generic function to predict a `bootstrap_phenologyFit` object.
#' 
#' @param object object of class `phenologyFit` to predict.
#' @param SeasonList List with data frames per season, see
#' \link{phenologyFit} for more details.
#' @param ... additional parameters, ignored here
#'
#' @return
#' A data.frame with one column `pbloomJDays` and a second one `Err`.
#' 
#' @export
predict.bootstrap_phenologyFit <- function(object, SeasonList, ...) {
  Y <- t(array(unlist(lapply(object$res,
                             FUN=function(x) return(predictBloomDays(par=x$par,
                                                                     SeasonList=SeasonList,
                                                                     modelfn=object$object$modelfn)))),
               dim=c(length(SeasonList), sum(object$boot.R))))
  return(data.frame(pbloomJDays=predict(object$object, SeasonList=SeasonList), Err=apply(Y, FUN=sd, MARGIN=2)))
}


#' Concatenate bootstrap_phenologyfit objects
#'
#' @param ... Zero or multiple objects of type `bootstrap_phenologyfit`.
#'
#' @return
#' An object of class `bootstrap_phenologyFit`, the concatenation of the
#' list of input object.
#' 
#' @export
c.bootstrap_phenologyFit <- function (...) {
  rval <- Reduce(concat_bootstrap_phenologyFit, list(...), NULL)
  return (invisible(rval))
}


concat_bootstrap_phenologyFit <- function(left, right) {
  if(is.null(left)) return(right)
  if(is.null(right)) return(left)
  stopifnot(inherits(left, "bootstrap_phenologyFit"))
  stopifnot(inherits(right, "bootstrap_phenologyFit"))
  stopifnot(left$seed != right$seed)
  left$boot.R <- c(left$boot.R, right$boot.R)
  left$res <- c(left$res, right$res)
  left$seed <- c(left$seed, right$seed)
  return(left)
}
