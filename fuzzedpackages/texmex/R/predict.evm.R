# Author: Harry Southworth
# Date: 2011-11-25
## Purpose: Create a predict method for objects of class evmOpt, evmSim
##          and evmBoot that
##          returns parameters, return levels or (maybe) return periods,
##          depending on arguments given.
#
# predict.evmOpt
# predict.evmSim
# predict.evmBoot
# rl
# rl.evm
# rl.evmSim
# rl.evmBoot
# linearPredictors
# linearPredictors.evm
# linearPredictors.evmSim
# linearPredictors.evmBoot

################################################################################
## evm



#' Predict return levels from extreme value models, or obtain the linear
#' predictors.
#'
#' Predict return levels from extreme value models, or obtain the linear
#' predictors.
#'
#' By default, return levels predicted from the unique values of the linear
#' predictors are returned. For the \code{evmBoot} method, estimates of
#' confidence intervals are simply quantiles of the bootstrap sample. The
#' \code{evmBoot} method is just a wrapper for the \code{evmSim} method.
#'
#' @param object An object of class \code{evmOpt}, \code{evmSim} or
#' \code{evmBoot}.
#' @param newdata The new data that you want to make the prediction for.
#' Defaults in \code{newdata = NULL} in which case the data used in fitting the
#' model will be used. Column names must match those of the original data
#' matrix used for model fitting.
#' @param type For the predict methods, the type of prediction, either "return
#' level" (or "rl") or "link" (or "lp"). Defaults to \code{type = "return
#' level"}. When a return level is wanted, the user can specify the associated
#' return period via the \code{M} argument. If \code{type = "link"} the linear
#' predictor(s) for \code{phi} and \code{xi} (or whatever other parameters are
#' in your \code{texmexFamily} are returned.
#'
#' For the plot methods for simulation based estimation of underlying
#' distributions i.e. objects derived from "evmSim" and "evmBoot" classes,
#' whether to use the sample median \code{type="median"} or mean
#' \code{type="mean"} estimate of the parameter.
#' @param se.fit Whether or not to return the standard error of the predicted
#' value. Defaults to \code{se.fit = FALSE} and is not implemented for
#' \code{predict.evmSim} or \code{predict.evmBoot}.
#' @param ci.fit Whether or not to return a confidence interval for the
#' predicted value. Defaults to \code{ci.fit = FALSE}. For objects of class
#' \code{evmOpt}, if set to \code{TRUE} then the confidence interval is a
#' simple symmetric confidence interval based on the estimated approximate
#' standard error. For the \code{evmSim} and \code{evmBoot} methods, the
#' confidence interval represents quantiles of the simulated distribution of
#' the parameters.
#' @param M The return period: units are number of observations. Defaults to
#' \code{M = 1000}. If a vector is passed, a list is returned, with items
#' corresponding to the different values of the vector \code{M}.
#' @param alpha If \code{ci.fit = TRUE}, a 100(1 - alpha)\% confidence interval
#' is returned. Defaults to \code{alpha = 0.050}.
#' @param unique.  If \code{unique. = TRUE}, predictions for only the unique
#' values of the linear predictors are returned, rather than for every row of
#' \code{newdata}. Defaults to \code{unique. = TRUE}.
#' @param all For the \code{evmSim} and \code{evmBoot} methods, if \code{all =
#' TRUE}, the predictions are returned for every simulated parameter vector.
#' Otherwise, only a summary of the posterior/bootstrap distribution is
#' returned. Defaults to \code{all = FALSE}.
#' @param full.cov Should the full covariance matrix be returned as part of a
#' \code{list} object. This is used internally and not intended for direct use.
#' Defaults to \code{full.cov = FALSE}
#' @param sumfun For the \code{evmSim} and \code{evmBoot} methods, a summary
#' function can be passed in. If \code{sumfun = FALSE}, the default, the
#' summary function used returns the estimated mean and median, and quantiles
#' implied by \code{alpha}.
#' @param x An object of class \code{lp.evmOpt}, \code{lp.evmSim} or
#' \code{lp.evmBoot}, to be passed to methods for these classes.
#' @param main,pch,ptcol,cex,linecol,cicol,polycol,plot,plot. Further arguments to plot
#' methods.
#' @param digits Number of digits to show when printing objects.
#' @param ... Further arguments to methods.
#' @return A list with two entries: the first being the call and the
#' second being a further list with one entry for each value of
#' \code{M}.
#' @note At present, the confidence intervals returned for an object of class
#' \code{evmOpt} are simple confidence intervals based on assumptions of
#' normality that are likely to be far from the truth in many cases. A better
#' approach would be to use profile likelihood, and we intend to implement this
#' method at a future date.  Alternatively, the credible intervals returned by
#' using Bayesian estimation and the predict method for class "evmSim" will
#' tend to give a better representation of the asymmetry of the estimated
#' intervals around the parameter point estimates.
#' @author Harry Southworth and Janet E. Heffernan
#' @keywords methods
#' @export
predict.evmOpt <-
    # Get predictions for an evm object. These can either be the linear predictors
    # or return levels.
function(object, M=1000, newdata=NULL, type="return level", se.fit=FALSE,
         ci.fit=FALSE, alpha=.050, unique.=TRUE, ...){
    theCall <- match.call()

    res <- list(obj = switch(type,
                             "rl"=, "return level" = rl.evmOpt(object, M, newdata,
                                                               se.fit=se.fit, ci.fit=ci.fit,
                                                               alpha=alpha, unique.=unique.),
                             "lp" =,"link" = linearPredictors.evmOpt(object, newdata, se.fit,
                                                                     ci.fit, alpha, unique.=unique.)
                             ),
                call = theCall)
    oldClass(res) <- class(res$obj)
    res$obj <- unclass(res$obj)
    res
}

## Linear predictor functions for GPD

##' @rdname predict.evmOpt
##' @export
linearPredictors.evmOpt <- function(object, newdata=NULL, se.fit=FALSE, ci.fit=FALSE,
                             alpha=.050, unique.=TRUE, full.cov=FALSE, ...){

    D <- texmexMakeNewdataD(object, newdata)

    if (unique.){
        z <- do.call('cbind', D)
        u <- !duplicated(z)
        D <- lapply(D, function(x, u) {
                           if(is.matrix(x[u,]))  x[u, ]
                           else if(ncol(x) == 1)  cbind(x[u,])
                           else t(cbind(x[u,]))
                       }, u=u )
    }

    res <- texmexMakeParams(coef(object), D)
    colnames(res) <- names(D)

    # Get the covariance matrices - one for every unique observation
    if(ci.fit | se.fit | full.cov){
      cov.se <- texmexMakeCovariance(object$cov, D)
      # Get standard errors
      ses <- t(sapply(cov.se, function(x){ sqrt(diag(x)) }))
      colnames(ses) <- paste(colnames(res), '.se', sep = '')
    }

    if (ci.fit){
        ci <- texmexMakeCI(res, ses, alpha)
        res <- cbind(res, ci)
    } # Close if(ci.fit

    if (se.fit){
        res <- cbind(res, ses)
    } # Close if(se.fit

    for (i in 1:length(D)){
      res <- addCov(res, D[[i]])
    }

    res <- list(link=res,family=object$family)

    if (full.cov){
        res$cov <- cov.se
    }

    oldClass(res) <- "lp.evmOpt"
    res
}

## Return level functions for GPD

## Reversing arguments M and newdata for anyone who wants to call these functions
## directly

## Will want to get return levels when using GEV rather than GPD, so make
## rl generic


#' Return levels
#'
#' Computation of return levels and confidence intervals for extreme
#' value models.
#'
#' @param object An object of class \code{evmOpt}, \code{evmSim} or
#'     \code{evmBoot}.
#' @param M The M-observation return level is computed by the
#'     function. Defaults to \code{M = 1000}.
#' @param newdata Data from which to calculate the return level. If
#'     not provided, the original data used to fit the model is used.
#'     Column names must match those of original data matrix used for
#'     model fitting.
#' @param se.fit Whether or not to return the standard error of the
#'     predicted value. Defaults to \code{se.fit = FALSE}.
#' @param ci.fit Whether or not to return a confidence interval for
#'     the predicted value. Defaults to \code{ci.fit = FALSE}. For
#'     objects of class \code{evmOpt}, if set to \code{TRUE} then the
#'     confidence interval is a simple symmetric confidence interval
#'     based on the estimated approximate standard error. For the
#'     \code{evmSim} and \code{evmBoot} methods, the confidence
#'     interval represents quantiles of the simulated distribution of
#'     the parameters.
#' @param alpha If \code{ci.fit = TRUE}, a 100(1 - alpha)\%
#'     confidence interval is returned. Defaults to \code{alpha =
#'     0.050}.
#' @param unique. If \code{unique. = TRUE}, predictions for only the
#'     unique values of the linear predictors are returned, rather
#'     than for every row of the original dataframe or of
#'     \code{newdata} if this latter is specified. Defaults to
#'     \code{unique. = TRUE}.
#' @param all For the \code{evmSim} and \code{evmBoot} methods, if
#'     \code{all = TRUE}, the predictions are returned for every
#'     simulated parameter vector. Otherwise, only a summary of the
#'     posterior/bootstrap distribution is returned. Defaults to
#'     \code{all = FALSE}.
#' @param sumfun For the \code{evmSim} and \code{evmBoot} methods, a
#'     summary function can be passed in. If \code{sumfun = FALSE},
#'     the default, the summary function used returns the estimated
#'     mean and median, and quantiles implied by \code{alpha}.
#' @param type For calls to plot methods for objects of class
#'     \code{rl.evmSim} or \code{rl.evmBoot}, specifies whether to
#'     use the sample mean (\code{type="mean"}) or median
#'     (\code{type="median"}) estimate of the return levels.
#' @param x Object passed to plot and print methods.
#' @param
#'     xlab,ylab,main,pch,ptcol,cex,linecol,cicol,polycol,smooth,sameAxes,ylim Further arguments to plot methods.
#' @param digits Number of digits to show when printing output.
#' @param ... Further arguments to be passed to methods.
#' @details The M-observation return level is defined as the value
#'     that is expected to be exceeded only once every M
#'     observations. Thus, it is an estimate of a high quantile of
#'     the fitted distribution.
#'
#' In models fit by the \code{evm} family of functions with
#' \code{family=gpd}, only a fraction of the data is actually
#' included in the model; the fitted GPD is a conditional model,
#' conditioning on the threshold having been exceeded. This
#' consideration is taken into account by \code{rl} which calculates
#' unconditional return levels from the entire distribution of
#' observations above and below the GPD fitting threshold.
#' @examples
#' mod <- evm(rain, qu=.8) # daily rainfall observations
#' rl(mod, M=100*365) # 100-year return level
#' @name rl
#' @export
rl <- function(object, M = 1000, newdata = NULL, se.fit = FALSE, ci.fit = FALSE, alpha = 0.050, unique. = TRUE, ...){
    UseMethod("rl")
}

#' @rdname predict.evmOpt
#' @export
linearPredictors <- function(object, newdata = NULL, se.fit = FALSE, ci.fit = FALSE, alpha = 0.050, unique. = TRUE, ...){
    UseMethod("linearPredictors")
}

#' @rdname rl
#' @export
rl.evmOpt <- function(object, M=1000, newdata=NULL, se.fit=FALSE, ci.fit=FALSE,
                       alpha=.050, unique.=TRUE, ...){
    co <- linearPredictors.evmOpt(object, newdata=newdata, unique.=unique., full.cov=TRUE)
    covs <- co$cov # list of covariance matrices, one for each (unique) observation
    co <- co$link
    X <- co[,-(1:length(object$data$D)), drop=FALSE]
    if(is.null(dim(X))){
      X <- matrix(X)
      dimnames(X) <- list(dimnames(co)[[1]], dimnames(co)[[2]][-(1:length(object$data$D))])
    }

    delta <- object$family$delta
    rl <- object$family$rl

    if (any(M * object$rate < 1)){
      stop("M * rate must be > 1 for every M")
    }

    res <- lapply(M, rl, param=co, model=object)

    getse <- function(o, co, M, delta, covs){
        dxm <- lapply(split(co, 1:nrow(co)), delta, m=M, model=o)

        # Get (4.15) of Coles, page 82, adjusted for phi = log(sigma)
        se <- sapply(1:length(covs),
                     function(i, dxm, covs){
                        covs <- covs[[i]]; dxm <- c(dxm[[i]])
                        sqrt(mahalanobis(dxm, center=rep(0, ncol(covs)), cov=covs, inverted=TRUE))
                     }, dxm=dxm, covs=covs)
        se
    }

    if (ci.fit){
        ci.fun <- function(i, object, co, M, res, alpha, delta, covs){
            wh <- res[[i]];
            se <- getse(object, co, M[i], delta=delta, covs=covs)
            lo <- wh - qnorm(1 - alpha/2)*se
            hi <- wh + qnorm(1 - alpha/2)*se
            wh <- cbind(wh, lo=lo, hi=hi)

            colnames(wh) <- c("RL", paste(100*alpha/2, "%", sep = ""),
                              paste(100*(1 - alpha/2), "%", sep = ""))
            wh
        } # ci.fun
        res <- lapply(1:length(M), ci.fun, object=object, co=co,
                                           M=M, res=res, alpha=alpha,
                                           delta=delta, covs=covs)
    } # Close if (ci.fit

    if (se.fit){
        se.fun <- function(i, object, co, M, res, alpha, delta, covs){
            wh <- res[[i]]
            se <- getse(object, co, M[i], delta=delta, covs=covs)
            wh <- cbind(RL=wh, se.fit=se)
            wh
        } # ci.fun
        res <- lapply(1:length(M), se.fun, object=object, co=co,
                                           M=M, res=res, alpha=alpha,
                                           delta=delta, covs=covs)
    }

    cov.fun <- function(i,res){
      wh <- res[[i]]
      wh <- addCov(wh,X)
      wh
    }
    res <- lapply(1:length(M), cov.fun,res=res)

    names(res) <- paste("M.", M, sep = "")
    oldClass(res) <- "rl.evmOpt"
    res
}

################################################################################
## evmSim

#' @rdname predict.evmOpt
#' @export
predict.evmSim <- function(object, M=1000, newdata=NULL, type="return level",
                         se.fit=FALSE, ci.fit=FALSE, alpha=.050, unique.=TRUE,
                         all=FALSE, sumfun=NULL, ...){
    theCall <- match.call()

    res <- list(obj = switch(type,
                             "rl" = , "return level" = rl.evmSim(object, M=M, newdata=newdata,
                                                                 se.fit=se.fit, ci.fit=ci.fit,
                                                                 alpha=alpha, unique.=unique., all=all,
                                                                 sumfun=sumfun,...),
                             "lp" = , "link" = linearPredictors.evmSim(object, newdata=newdata,
                                                                       se.fit=se.fit, ci.fit=ci.fit,
                                                                       alpha=alpha, unique.=unique., all=all,
                                                                       sumfun=sumfun,...)
                             ),
                call = theCall)
    oldClass(res) <- class(res$obj)
    res$obj <- unclass(res$obj)
    res
}

#' @rdname predict.evmOpt
#' @export
linearPredictors.evmSim <- function(object, newdata=NULL, se.fit=FALSE, ci.fit=FALSE,
                                     alpha=.050, unique.=TRUE, all=FALSE, sumfun=NULL, ...){
    if (se.fit){ warning("se.fit not implemented - ignoring") }
    D <- texmexMakeNewdataD(object$map, newdata)

    X.all <- do.call("cbind", D)
    ModelHasCovs <- ncol(X.all) > length(D)

    if(ModelHasCovs){
      covCols <- apply(X.all, 2, function(x) !all(x==1))
      Xnames <- colnames(X.all)
      X.all <- X.all[, covCols, drop=FALSE]
      colnames(X.all) <- Xnames[covCols]
    }

    if (unique.){
        X.all <- unique(X.all)
        u <- !duplicated(do.call("cbind", D))
        D <- lapply(D, function(x, u){ x[u,, drop=FALSE] }, u=u)
    }

    # Get matrices of parameters (i.e. split full parameter matrix into phi, xi whatever)
    param <- texmexGetParam(D, object$param)

    # Get linear predictors
    res <- lapply(1:nrow(D[[1]]), # For each observation get matrix of parameters
              function(i, x, p){
                  wh <- lapply(1:length(D),
                               function(j, x, p, i){
                                   rowSums(t(t(p[[j]]) * c(x[[j]][i, ])))
                               }, x=x, p=p, i=i)
                  wh <- do.call("cbind", wh)
                  colnames(wh) <- names(x)
                  wh
                }, x=D, p=param)
    # res should be a list containing a matrix for each observation.
    # The matrix represents the simulated posterior, one column for each
    # major parameter (i.e. linear predictors)

    ############################################################################
    ## Hard part should be done now. Just need to summarize

    if (ci.fit){
        # Need to get names by pasting together CI names and parameter names
        wh <- texmexMakeCISim(res[[1]], alpha=alpha, object=object, sumfun=sumfun)
        wh <- colnames(wh)
        wh <- paste(rep(names(D), ea=length(wh)), wh, sep = ":")

        # Need to get order of output correct, so need to faff about with transposing
        res <- sapply(res, function(x){
                               t(texmexMakeCISim(x, alpha=alpha, object=object, sumfun=sumfun))
                           })
        res <- t(res)
        colnames(res) <- wh
    } else if (all){
      res <- res
    } else { # Just point estimates
        res <- t(sapply(res, function(x){ apply(x, 2, mean) }))
    }

    if(!all){
      if(ModelHasCovs){
        for (i in 1:length(D)){
            res <- addCov(res,D[[i]])
        }
      }
    }
    else {
        if (ModelHasCovs & nrow(X.all) != length(res)){
            stop("Number of unique combinations of covariates doesn't match the number of parameters")
        }
        for (i in 1:length(res)){
          if(ModelHasCovs){
            res[[i]] <- cbind(res[[i]], matrix(rep(X.all[i,], nrow(res[[i]])),
                                               nrow=nrow(res[[i]]), byrow=TRUE))
            colnames(res[[i]]) <- c(names(D), colnames(X.all))
          } else {
            colnames(res[[i]]) <- names(D)
          }
        }
    }

    res <- list(link=res,family=object$map$family)
    oldClass(res) <- "lp.evmSim"
    res
}


#' @rdname rl
#' @export
rl.evmSim <- function(object, M=1000, newdata=NULL, se.fit=FALSE, ci.fit=FALSE, alpha=.050, unique.=TRUE, all=FALSE, sumfun=NULL,...){
    if (se.fit){ warning("se.fit not implemented") }

    co <- linearPredictors.evmSim(object, newdata=newdata, unique.=unique., all=TRUE, sumfun=NULL)$link
    # XXX Next line seems silly! Why not compute it from the line above?
    Covs <- linearPredictors.evmSim(object, newdata=newdata, unique.=unique., sumfun=NULL)$link
    X <- Covs[,-(1:length(object$map$data$D))]
    if(is.null(dim(X))){
      X <- matrix(X)
      dimnames(X) <- list(dimnames(Covs)[[1]],dimnames(Covs)[[2]][-(1:length(object$map$data$D))])
    }

    sim.rl <- function(m, param, model){
        rl <- model$family$rl
        rl(m=m, param, model)
    }

    # co is a list with one element for each unique item in
    # new data. Need to loop over vector M and the elements of co

    getrl <- function(m, co, ci.fit, alpha, all, object){
        res <- sapply(co, sim.rl, m=m, model=object$map)
        if (ci.fit){
            res <- texmexMakeCISim(res, alpha, object$map, sumfun, M=m)
        } # Close if (ci.fit
        else if (!all){
            res <- apply(res, 2, mean)
        }
        res
    }

    res <- lapply(M, getrl, co=co, ci.fit=ci.fit, alpha=alpha, all=all, object=object)

    if(!all){
      cov.fun <- function(i,res){
        wh <- res[[i]]
        wh <- addCov(wh,X)
        wh
      }
      res <- lapply(1:length(M), cov.fun,res=res)
    }

    names(res) <- paste("M.", M, sep = "")
    oldClass(res) <- "rl.evmSim"
    res
}

################################################################################
# evmBoot
#' @rdname predict.evmOpt
#' @export
predict.evmBoot <- function(object, M=1000, newdata=NULL, type="return level",
                            se.fit=FALSE, ci.fit=FALSE, alpha=.050, unique.=TRUE,
                            all=FALSE, sumfun=NULL, ...){
    theCall <- match.call()

    res <- list(obj = switch(type,
                             "rl" = , "return level" = rl.evmBoot(object, newdata=newdata, M=M,
                                                                  se.fit=se.fit, ci.fit=ci.fit,
                                                                  alpha=alpha, unique.=unique.,
                                                                  all=all, sumfun=sumfun,...),
                             "lp" = , "link" = linearPredictors.evmBoot(object, newdata=newdata,
                                                                  se.fit=se.fit, ci.fit=ci.fit,
                                                                  alpha=alpha, unique.=unique.,
                                                                  all=all, sumfun=sumfun,...)
                             ),
                call = theCall)
    oldClass(res) <- class(res$obj)
    res$obj <- unclass(res$obj)
    res
}

namesBoot2sim <- function(bootobject){
    names(bootobject) <- c("call", "param", "map")
    bootobject
}

#' @rdname predict.evmOpt
#' @export
linearPredictors.evmBoot <- function(object, newdata=NULL, se.fit=FALSE, ci.fit=FALSE, alpha=.050,
                                 unique.=TRUE, all=FALSE, sumfun=NULL,...){
    # This should just be the same as for an evmSim object, but some
    # names and stuff are different.
  object <- namesBoot2sim(object)
  res <- linearPredictors.evmSim(object, newdata=newdata, se.fit=se.fit, ci.fit=ci.fit, all=all, unique.=unique., alpha=alpha, sumfun=sumfun,...)
  oldClass(res) <- "lp.evmBoot"
  res
}

#' @rdname rl
#' @export
rl.evmBoot <- function(object, M=1000, newdata=NULL, se.fit=FALSE, ci.fit=FALSE, alpha=0.050, unique.=TRUE, all=FALSE, sumfun=NULL,...){
    # This should just be the same as for an evmSim object, but some
    # names are different.
  object <- namesBoot2sim(object)
  res <- rl.evmSim(object, M=M, newdata=newdata, se.fit=se.fit, ci.fit=ci.fit,alpha=alpha, unique.=unique., all=all, sumfun=sumfun,...)
  oldClass(res) <- "rl.evmBoot"
  res
}

################################################################################
# Method functions

#' @export
#' @rdname rl
print.rl.evmOpt <- function(x, digits=3, ...){
	if(is.null(x$obj)) x <- list(obj=x)
    nms <- names(x$obj)
    newnms <- paste("M =", substring(nms, 3), "predicted return level: ")
    lapply(1:length(x$obj), function(i, o, title){
                                 cat(title[i])
    	                         temp <- o[[i]]
    	                         names(temp) <- NULL
    	                         cat(signif(temp,digits=digits))
                                 cat("\n")
                                 NULL}, o=x$obj, title=newnms)
    invisible(x)
}

#' @export
print.rl.evmSim    <- print.rl.evmOpt
#' @export
print.rl.evmBoot <- print.rl.evmOpt



#' @export
#' @rdname predict.evmOpt
print.lp.evmOpt <- function(x, digits=3, ...){
    cat("Linear predictors:\n")
    print(unclass(x$obj$link), digits=3,...)
    invisible(x)
}

#' @export
print.lp.evmSim <- print.lp.evmOpt
#' @export
print.lp.evmBoot <- print.lp.evmOpt

