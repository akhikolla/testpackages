#' @include tools.R
NULL

## 2016-11-07 (new) a helper function to print 'call' objects since the call
##
##      cat("\nCall:", deparse(x$call, width.cutoff = width), "", sep = "\n")
##
## doesn't necessarilly do the job (argument 'width.cutoff' sets a lower
## limit and after that seems to look for the next feasible break point.
##
## The solution below ignores argument 'width' (since print.default() doesn't have it)
.deparseCall <- function(call, width = getOption("width")){
    wrk <- deparse(call, width.cutoff = floor(width * 0.85))
    if(max(nchar(wrk)) <= width)
        wrk
    else # fall back on print(); it should fit in options("width") if possible.
        capture.output(print(call))
}

## =============================================================================
## ----------------------------- declare generic -------------------------------
## =============================================================================
#' Extract Standard Errors of Model Coefficients
#'
#' Extract standard errors of model coefficients from objects returned by
#' count-modeling functions.
#'
#' @param object object returned by one of the count-modeling functions
#' @param parm parameter's name or index
#' @param type type of standard error: asymtotic normal standard errors
#'     (\code{"asymptotic"}) or bootsrap (\code{"boot"}).
#' @param ... further arguments for methods.
#'
#' @return a named numeric vector
#' @export
se.coef <- function(object, parm, type, ...) {
    UseMethod("se.coef", object)
}

## ==============================================================================
## ------------------------------ methods definition ----------------------------
## ==============================================================================

#' Methods for renewal objects
#'
#' Methods for renewal objects.
#'
#' Objects from class \code{"renewal"} represent fitted count renewal models and
#' are created by calls to \code{"renewalCount()"}. There are methods for this class
#' for many of the familiar functions for interacting with fitted models.
#'
#' @param object an object from class \code{"renewal"}.
#' @param ... further arguments for methods
#' @param type,parm,level,bootType,x,digits see the corresponding generics and section
#'     Details.
#'
#' @examples
#' fn <- system.file("extdata", "McShane_Wei_results_boot.RDS", package = "Countr")
#' object <- readRDS(fn)
#' class(object) # "renewal"
#'
#' coef(object)
#' vcov(object)
#'
#' ## Pearson residuals: rescaled by sd
#' head(residuals(object, "pearson"))
#' ## response residuals: not rescaled
#' head(residuals(object, "response"))
#'
#' head(fitted(object))
#'
#' ## loglik, nobs, AIC, BIC
#' c(loglik = as.numeric(logLik(object)), nobs = nobs(object),
#'   AIC = AIC(object), BIC = BIC(object))
#'
#' asym <- se.coef(object, , "asymptotic")
#' boot <- se.coef(object, , "boot")
#' cbind(asym, boot)
#' @name renewal_methods
NULL

#' @rdname renewal_methods
#' @method coef renewal
#' @export
coef.renewal <- function(object, ...) {
    object$coefficients
}

#' @rdname renewal_methods
#' @method vcov renewal
## #' @param object object from class \code{renewal()}
## #' @param ... extra parameters to be passed. user can pass parameter
## #' \code{method} with option \code{asymptotic} (numerical hessian) or
## #' \code{boot} (bootsrap). The other parameters will be passed to the \code{boot}
## #' function called inside \code{addBootSampleObject()}.
#' @export
vcov.renewal <- function(object, ...) {
    v <- object$vcov
    pars <- list(...)
    method <- ifelse(is.null(pars$method), "asymptotic", pars$method)
    pars$method <- NULL
    theta_ml <- coef(object)
    nmPars <- names(theta_ml)

    if (nrow(v) == length(theta_ml))
        return(v)
    else {
        if (method == "asymptotic") {
            .obj <- function(theta) sum(object$score_fct(theta))

            hess <- numDeriv::hessian(.obj, theta_ml)
            varCovarcount <- try(-solve(hess))
            if ((inherits(varCovarcount, "try-error")) |
                (any(diag(varCovarcount) < 0)) ) {
                varCovarcount <- Matrix::nearPD(-ginv(hess))$mat
                warning(paste("variance-covariance matrix was computed",
                              "by smoothing the genralized inverse hessian !"))
            }

            dimnames(varCovarcount) <- list(nmPars, nmPars)
        } else if (method == "boot") { ## bootsrap
            if (is.null(object$boot)) {
                warning("boot sample not found. it will be created ...")
                pars$R <- ifelse(is.null(pars$R), 100, pars$R)
                object <- addBootSampleObject(object, pars)
            }

            ## compute the variance-covariance matrix
            b <- object$boot
            th_hat <- t(b$t)
            th_bar <- colMeans(b$t)
            th_bar_mat <- matrix(th_bar, nrow = length(th_bar), ncol = b$R)

            th <- th_hat - th_bar_mat
            varCovarcount <- (th %*% t(th)) / (b$R - 1)
            dimnames(varCovarcount) <- list(nmPars, nmPars)
        }
    }

    return(varCovarcount)
}

#' @rdname renewal_methods
#' @method residuals renewal
#' @export
residuals.renewal <- function (object, type = c("pearson", "response",
                                                "prob"), ...) {
    type <- match.arg(type)
    res <- object$residuals
    switch(type,
           "response" = {
               out <- res
           },
           "pearson" = {
               out <- res / sqrt(object$wi)
           }
           )

    out
}

#' Method to visualise the residuals
#'
#' A method to visualise the residuals
#'
#' @param object object returned by one of the count-modeling functions
#' @param type character type of residuals to be used.
#' @param ... further arguments for methods.
#'
#' @export
residuals_plot <- function(object,  type, ...) {
    UseMethod("residuals_plot", object)
}

#' @rdname renewal_methods
#' @method residuals_plot renewal
#' @export
residuals_plot.renewal <- function(object, type = c("pearson", "response",
                                                    "prob"), ...) {
    type <- match.arg(type)
    res <- residuals(object, type)

    ## Uses ideas from Cameron (2013) Regression Analysis ofcount data (CH 5)
    par(mfrow = c(2, 2))

    ## plot estimated density
    plot(density(res, ...), ylab = paste(type, "residuals density"),
         main = paste(type, "residuals density"))

    ## residuals qq plot
    qqnorm(res, ylim = range(res), main = "Residuals qq-plot")
    qqline(res, ylim = range(res))

    ## residuals plotted against predicted values of the dependent variable
    ## i.e residuals plotted against the predicted mean Yhat
    plot(x = object$fitted.values,  y = res,
         xlab = "predicted mean value", ylab = paste(type, "residuals"),
         main = paste(type, "residuals", "VS", "predicted mean")
         )

    ## residuals plotted against regressors score
    ## (included regressors to model the scale parameter)
    ## vecDistParsList$scale
    plot(x = object$vecDistParsList$scale,  y = res,
         xlab = "predicted scale", ylab = paste(type, "residuals"),
         main = paste(type, "residuals", "VS", "predicted scale")
         )
}


#' @rdname renewal_methods
#' @method fitted renewal
#' @export
fitted.renewal <- function (object, ...) {
    object$fitted.values
}

#' % Extract Standard Errors of Model Coefficients
#'
#' % Extract Standard Errors of Model Coefficients
#'
#' The method for class \code{"renewal"} extracts standard errors of model
#' coefficients from objects returned by \code{renewal}.  When bootsrap standard
#' error are requested, the function checks for the bootsrap sample in
#' \code{object}. If it is not found, the bootsrap sample is created and a
#' warning is issued. Users can choose between asymtotic normal standard errors
#' (\code{asymptotic}) or bootsrap (\code{boot}).
#'
#' % param object object returned by \code{renewal}.
#' % param parm parameters name or index
#' % param type type of standard error. User can choose between asymtotic normal
#' %    standard errors (\code{asymptotic}) or bootsrap (\code{boot}).
#' @inheritParams se.coef
#'
#' @examples
#' ## see examples for renewal_methods
#' @rdname se.coef
#' @export
se.coef.renewal <- function(object, parm, type = c("asymptotic", "boot"), ...) {
    type <- match.arg(type)
    switch(type,
           "asymptotic" = {
               se <- sqrt(diag(vcov(object)))
               names(se) <- names(object$coefficients)
           },
           "boot" = {
               if (is.null(object$boot)) {
                   cat("Precomputed boot sample not found, creating a new one.\n")
                   cat("\tThis may take some time...\n")

                   pars <- list(...)
                   pars$R <- ifelse(is.null(pars$R), 100, pars$R)
                   object <- addBootSampleObject(object, pars)
               }
               ## extract se from car::summary.boot
               se <- summary(object$boot)[, "bootSE"]
           }
           )
    names(se) <- names(coef(object))

    ind <- seq(along = se)
    if (!missing(parm))
        ind <- which(parm %in% names(se))

    return(se[parm])
}


#' @method confint renewal
#' @examples
#' ## CI for coefficients
#' asym <- confint(object, type = "asymptotic")
#' ## Commenting out for now, see the nite in the code of confint.renewal():
#' ## boot <- confint(object, type = "boot", bootType = "norm")
#' ## list(asym = asym, boot = boot)
#' @rdname renewal_methods
#' @export
confint.renewal <- function(object, parm, level = 0.95,
                            type = c("asymptotic", "boot"),
                            bootType = c("norm", "bca", "basic", "perc"), ...) {
    type <- match.arg(type)
    bootType <- match.arg(bootType)
    switch(type,
           "asymptotic" = {
               return(confint.default(object, parm, level, ...))
           },
           "boot" = {
               pars <- list(...)
               if (is.null(object$boot)) {
                   warning("boot sample not found. it will be created ...")
                   pars$R <- ifelse(is.null(pars$R), 100, pars$R)
                   object <- addBootSampleObject(object, pars)
               }

               ci <- confint(object$boot, parm, level, bootType)
               rownames(ci) <- names(coef(object))
               return(ci)
           }
           )
}

#' @method summary renewal
#' @examples
#' summary(object)
#' @rdname renewal_methods
#' @export
summary.renewal <- function(object, ...) {
    object$residuals <- residuals(object, type = "pearson")
    se <- sqrt(diag(vcov(object)))
    coef <- object$coefficients

    zstat <- coef/se
    pval <- 2 * pnorm(-abs(zstat))
    coef <- cbind(coef, se, zstat, pval)
    colnames(coef) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
    object$coefficients <- coef
    object$fitted.values <- object$model <- object$y <- object$x <- NULL
    object$start <- NULL
    class(object) <- "summary.renewal"
    object
}

#' @method print renewal
#' @examples
#' print(object)
#' @rdname renewal_methods
#' @export
#'
print.renewal <- function(x, digits = max(3, getOption("digits") - 3), ...) {
    ## see not at .deparseCall()
    ##    cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") *
    ##                                                         0.85)), "", sep = "\n")
    cat("\nCall:", sep = "\n")
    cat(.deparseCall(x$call), sep = "\n")

    ##  2017-10-25 Note by Georgi: Print the info also when optim. did not
    ##                             converge, with a note after that.
    ##        if (!x$converged) {
    ##            cat("optimisation did not converge\n")
    ##        } else {
        ##--------------- prepare links char
        textLink <- .summarizeLinkInformation(x$link)
        cat(paste0("\nCount model coefficients (inter-arrival ", x$dist,
                   " with ", textLink, "):\n"))
        print.default(format(x$coefficients, digits = digits),
                      print.gap = 2, quote = FALSE)
        cat("\nLog-likelihood:", sprintf(paste0("%.", digits, "f"), x$loglik),
            "on", x$n - x$df.residual, "Df\n")
    ##        }

    if (!x$converged) {
        cat("Note: optimisation did not converge.\n")
    }
    invisible(x)
}

#' @rdname renewal_methods
#' @method print summary.renewal
#' @param width numeric width length
#' @export
print.summary.renewal <- function(x, digits = max(3, getOption("digits") - 3),
                                  width = getOption("width"),
                                  ...) {
    ## 2016-11-07 Note by Georgi: argument 'width.cutoff' seems to be a lower
    ##            limit, it seems to look for the next feasible break point
    ##      cat("\nCall:", deparse(x$call, width.cutoff = width), "", sep = "\n")
    ##
    ## print(summary(gamModel), width = 300)

    cat("\nCall:", sep = "\n")
    cat(.deparseCall(x$call, width), sep = "\n")
    cat("\n")

    if (!x$converged) {
        cat("model did not converge\n")
    } else {
        ## residuals print
        cat("Pearson residuals:\n")
        print(structure(quantile(x$residuals, na.rm = TRUE), names = c("Min",
            "1Q", "Median", "3Q", "Max")), digits = digits, ...)

        ## dist & link functions
        ##--------------- prepare links char
        textLink <- .summarizeLinkInformation(x$link)

        ## cat(paste0("\nCount model coefficients (inter-arrival ", x$dist,
        ##            " with ", textLink, "):\n"))
        cat("Inter-arrival dist.:", x$dist, "\n")
        cat("              Links:", textLink, "\n")

        cat("\nCount model coefficients\n")

        ##---------------- print coefficients
        printCoefmat(x$coefficients, digits = digits, signif.legend = FALSE)
        if (getOption("show.signif.stars") &
            any(rbind(x$coefficients)[, 4] < 0.1, na.rm = TRUE)
            )
            cat("---\nSignif. codes: ",
                "0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1",
                "\n")

        ##------------------ optimization routine
        cat(paste("\nNumber of iterations in", x$method, "optimization:",
                  x$iterations, "\n"))

        ##-------------------- exec time
        cat(paste("\nExecution time", formatC(round(x$execTime), 3), "\n"))

        ##------------------- log-likelihood
        cat("Log-likelihood:", sprintf(paste0("%.", digits, "f"), x$loglik),
            "on", x$n - x$df.residual, "Df\n")
    }
    invisible(x)
}

#' @rdname renewal_methods
#' @method model.matrix renewal
#' @export
model.matrix.renewal <- function(object, ...) {
    if (!is.null(object$x))
        return(object$x)
    else if (!is.null(object$model))
        return(.getModelMatrix(formula = Formula(object$formula),
                               dist = object$dist,
                               mf = object$model, anc = object$anc,
                               raw_data = object$data,
                               standardise = object$standardise,
                               standardise_scale = object$standardise_scale, ...)
               )
    else
        stop("model should be saved to create model.matrix !")
}

## \code{logLik} method for class \code{renewal}
## @param object an object from class \code{renewal}
## @param ... not used
#' @method logLik renewal
#' @rdname renewal_methods
#' @export
logLik.renewal <- function (object, ...) {
    structure(object$loglik, df = object$n - object$df.residual,
              class = "logLik")
}

## \code{nobs} method for class \code{renewal}
## @param object an object from class \code{renewal}
## @param ... not used
#' @method nobs renewal
#' @examples
#' ## see renewal_methods
#' @rdname renewal_methods
#' @export
nobs.renewal <- function (object, ...) {
    object$n
}

#' % extractAIC method for class renewal
#'
#' % extractAIC method for class renewal
#'
#' % param fit an object from class \code{renewal}
#' % param scale TODO
#' % param k TODO
#' % param ... not used
#' @param fit,scale,k same as in the generic.
#'
#' @method extractAIC renewal
#' @examples
#' ## see renewal_methods
#' @rdname renewal_methods
#' @export
extractAIC.renewal <- function (fit, scale, k = 2, ...) {
     c(attr(logLik(fit), "df"), AIC(fit, k = k))
}

#' Predict method for renewal objects
#'
#' Compute predictions from renewal objects.
#'
#' @param type type of prediction.  If equal to \code{"response"}, give the mean
#'     probability associated with the individual covariates. If \code{"prob"},
#'     give the probability of the observed count.
#' @param time TODO
#' @inheritParams stats::predict.lm
#' @method predict renewal
#' @examples
#' fn <- system.file("extdata", "McShane_Wei_results_boot.RDS", package = "Countr")
#' object <- readRDS(fn)
#' data <- object$data
#' ## old data
#' predOld.response <- predict(object, type = "response", se.fit = TRUE)
#' predOld.prob <- predict(object, type = "prob", se.fit = TRUE)
#'
#' ## newData (extracted from old Data)
#' newData <- head(data)
#' predNew.response <- predict(object, newdata = newData,
#'                             type = "response", se.fit = TRUE)
#' predNew.prob <- predict(object, newdata = newData,
#'                         type = "prob", se.fit = TRUE)
#'
#' cbind(head(predOld.response$values),
#'            head(predOld.response$se$scale),
#'            head(predOld.response$se$shape),
#'            predNew.response$values,
#'            predNew.response$se$scale,
#'            predNew.response$se$shape)
#'
#' cbind(head(predOld.prob$values),
#'       head(predOld.prob$se$scale),
#'       head(predOld.prob$se$shape),
#'       predNew.prob$values,
#'       predNew.prob$se$scale,
#'       predNew.prob$se$shape)
#' @export
predict.renewal <- function(object, newdata = NULL, type = c("response", "prob"),
                            se.fit = FALSE, terms = NULL, na.action = na.pass,
                            time = 1.0, ...) {

    type <- match.arg(type)
    dist <- object$dist

    ## custom params
    customPars <- object$customPars

    ## link list
    linkList <- object$link

    ## check convolution parameters
    ## convPars <- renewal.convPars(list(...)$convPars, object$dist)
    convPars <- object$convPars

    ## extract standardise
    standardise <- ifelse(is.null(object$standardise), FALSE,
                          object$standardise)

    standardise_scale <- ifelse(is.null(object$standardise_scale), 1,
                                object$standardise_scale)


    ## prepare the modelMatrixList object as well as Y (response)
    if (is.null(newdata)) { ## no data provided
        ## build the modelMatrixList object
        if (!is.null(object$x))
            modelMatrixList <- object$x ## mf
        else if (!is.null(object$model)) {
            modelMatrixList <-
                .getModelMatrix(formula = Formula(object$formula),
                                dist = object$dist,
                                mf = object$model,
                                anc = object$anc,
                                raw_data = object$data,
                                standardise = standardise,
                                standardise_scale = standardise_scale,
                                customPars)
        } else
            stop(paste("predicted probabilities cannot be",
                       "computed with missing newdata")
                 )

        if (type == "prob") {
            ## check that the response is found in the object
            Y <- object$y
            if (is.null(Y))
                stop(paste("response should be saved in the fit object for",
                           "probability predictions !"))
        } else {
            out <- object$fitted.values
            se_out <- NA
            if (se.fit) {
                ## C <- .modelData(modelMatrixList, object$dist, customPars)
                ## se <- sqrt(diag(C %*% vcov(object) %*% t(C)))
                ## se_out <- .transformSE(se, object$link)
                se_out <- .getPredictionStd(modelMatrixList, vcov(object),
                                            object$dist, object$link, customPars)
            }

            return(list(values = as.numeric(out), se = se_out))
        }
    } else { ## data provided
        Fform <- Formula(object$formula)
        mf <- model.frame(Fform, data = newdata, na.action = na.action)
        Y <- model.response(mf)

        if (standardise)
            modelMatrixList <- predict(object$x, newdata)
        else
            modelMatrixList <-
                .getModelMatrix(formula = Fform, dist = object$dist,
                                mf = mf, anc = object$anc,
                                raw_data = newdata,
                                standardise = standardise,
                                standardise_scale = standardise_scale,
                                customPars)
    }

    key <- type == "response"

    out <-  .objectiveFunction(coef(object), dist, modelMatrixList,
                               linkList, time, convPars, Y, NULL, Ev = key,
                               seriesPars = object$seriesPars,
                               weiMethod = object$weiMethod,
                               summa = FALSE, customPars)

    if (type == "response") ## extract the mean response
        out <- sapply(out, .extractElem, ind = "ExpectedValue")

    se_out <- NA
    if (se.fit)
        se_out <- .getPredictionStd(modelMatrixList, vcov(object),
                                    object$dist, object$link, customPars)

    return(list(values = as.numeric(out), se = se_out))
}

#' Create a bootsrap sample for coefficient estimates
#'
#' Create a boostrap sample from coefficient estimates.
#'
#' The information in \code{object} is used to prepare the arguments and then
#' \code{boot} is called to generate the bootstrap sample.
#' The bootstrap sample is stored in \code{object} as component \code{"boot"}.
#' Arguments in \code{"..."} can be used customise the \code{boot()} call.
#'
#' @param object an object to add boot object to
#' @param ... extra parameters to be passed to the \code{boot::boot()} function
#'     other than \code{data} and \code{statistic}.
#' @return \code{object} with additional component \code{"boot"}
#'
#' @seealso \code{\link{renewal_methods}}
#' @examples
#' ## see renewal_methods
#' @export
#' @importFrom boot boot
addBootSampleObject <- function(object, ...) {
    UseMethod("addBootSampleObject", object)
}

#' @rdname renewal_methods
#' @method addBootSampleObject renewal
#' @export
addBootSampleObject.renewal <- function(object, ...) {
    formula <- object$formula
    data <- object$data
    weights <- object$weights
    dist <- object$dist
    anc <- object$anc
    convPars <- object$convPars
    link <- object$link
    time <- object$time
    control <- object$control
    control$trace <- FALSE
    customPars <- object$customPars
    standardise <- ifelse(is.null(object$standardise), FALSE, object$standardise)
    standardise_scale <- ifelse(is.null(object$standardise_scale), 1,
                                object$standardise_scale)

    bootIter <- 0

    bootFun <- function(data, indices) {
        bdata <- data[indices, , drop = FALSE]
        res <- renewalCount(formula = formula, data = bdata, weights = weights,
                            dist = dist, anc = anc, convPars = convPars,
                            link = link, time = time, computeHessian = FALSE,
                            control = control, customPars = customPars,
                            standardise = standardise,
                            standardise_scale = standardise_scale)
        coef(res)
    }

    bootList <- list(...)
    if (is.list(bootList) & length(bootList) == 1)
        bootList <- bootList[[1]]

    bootList$data <- data
    bootList$statistic <- bootFun

    if (!"R" %in% names(bootList))
        bootList$R <- 100

    object$boot <- do.call(boot, bootList)
    object
}

#' @rdname renewal_methods
#' @method df.residual renewal
#' @export
df.residual.renewal <- function(object, ...) {
    object$df.residual
}

## copied update.default and modified it.
update.renewal <- function(object, formula., anc, ..., evaluate = TRUE) {
    if (is.null(call <- getCall(object)))
        stop("need an object with call component")
    extras <- match.call(expand.dots = FALSE)$...
    if (!missing(formula.))
        call$formula <- update.formula(formula(object), formula.)
    if (length(extras)) {
        existing <- !is.na(match(names(extras), names(call)))
        for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
        if (any(!existing)) {
            call <- c(as.list(call), extras[!existing])
            call <- as.call(call)
        }
    }

    ## process argument `anc' - this is the chunk specific to 'renewal'
    if(!missing(anc)  &&  !is.null(anc)) {
        newAnc <- call$anc # TODO: check if need as.list or whatever
        if(is.null(newAnc)) {
            call$anc <- anc
        } else {
            existingAncNames <- names(newAnc)
            ## TODO: for speed use match()
            for(nam in names(anc)) {
                if(nam %in% existingAncNames) {
                    newAnc[[nam]] <- update.formula(newAnc[[nam]], anc[[nam]])
                } else {
                    newAnc[[nam]] <- anc[[nam]]
                }
            }
            call$anc <- newAnc
        }
    }

    if (evaluate)
        eval(call, parent.frame())
    else call
}

## =============================================================================
## ------------------------- goodness of fit diagnosis -------------------------
## =============================================================================
#' Pearson Chi-Square test
#'
#' Carry out Pearson Chi-Square test and compute the Pearson statistic.
#'
#' The computation is inspired from Cameron(2013) Chapter 5.3.4. Observed and
#' fitted frequencies are computed and the contribution of every observed cell
#' to the Pearson's chi-square test statistic is reported. The idea is to check if
#' the fitted model has a tendancy to over or under predict some ranges of data
#'
#' @param object an object from class \code{renewal}
#' @param ... currently not used
#' @return data.frame with 5 columns given the count values (\code{Counts}),
#'     observed frequencies (\code{Actual}), model's prediction
#'     (\code{Predicted}), the difference (\code{Diff}) and the contribution to
#'     the Pearson's statistic (\code{Pearson}).
#' @seealso \code{\link{chiSq_gof}}
#' @references \insertRef{cameron2013regression}{Countr}
#' @export
chiSq_pearson <- function(object, ...) {
    UseMethod("chiSq_pearson", object)
}

#' @rdname chiSq_pearson
#' @method chiSq_pearson renewal
#' @export
chiSq_pearson.renewal <- function (object, ...) {
    lhs <- as.character(object$formula)[2]

    ## compute empirical frequencies
    tab <- table(object$data[[lhs]])
    count <- as.numeric(names(tab))
    obs_freq <- as.numeric(tab) / object$n

    ## compute fitted frequncy
    data_tmp <- object$data
    .computeFreq <- function(count_value) {
        data_tmp[[lhs]] <- count_value
        mean(predict(object, newdata = data_tmp, type = "prob")$values,
             na.rm = TRUE)
    }

    fitted_freq <- sapply(count, .computeFreq)

    ## ## adjust last cell: predict 1 - sum of previous cells
    ## fitted_freq[length(fitted_freq)] <-
    ##     1 - sum(fitted_freq[1:(length(fitted_freq) - 1)])

    ## error computation
    diff <- abs(obs_freq - fitted_freq)
    pearson <- object$n * (fitted_freq - obs_freq)^2 / fitted_freq

    out <- data.frame(Counts = count,
                      Actual = obs_freq,
                      Predicted = fitted_freq,
                      Diff = diff,
                      Pearson = pearson
                      )
    out
}


#' @rdname chiSq_pearson
#' @method chiSq_pearson glm
#' @importFrom pscl predprob
#' @export
chiSq_pearson.glm <- function (object, ...) {
    ## compute empirical frequencies
    tab <- table(object$y)
    count <- as.numeric(names(tab))
    obs_freq <- as.numeric(tab) / length(object$y)

    ## compute fitted frequncy
    pbs <- pscl::predprob(object)[, as.character(count)]
    fitted_freq <- colMeans(pbs)

    ## error computation
    diff <- abs(obs_freq - fitted_freq)
    pearson <- object$n * (fitted_freq - obs_freq)^2 / fitted_freq

    out <- data.frame(Counts = count,
                      Actual = obs_freq,
                      Predicted = fitted_freq,
                      Diff = diff,
                      Pearson = pearson
                      )
    out
}


#' Formal Chi-square goodness-of-fit test
#'
#' Carry out the formal chi-square goodness-of-fit test described by Cameron
#' (2013).
#'
#' The test is a conditional moment test described in details in Cameron (2013,
#' Section 5.3.4). We compute the asymptotically equivalent outer product of the
#' gradient version which is justified for renewal models (fully parametric +
#' parameters based on MLE).
#' @param object an object from class \code{renewal}
#' @param breaks integer values at which the breaks shoudl happen. The function
#'     will compute the observed frequencies in the intervals \code{[breaks[i],
#'     breaks[i + 1])}.
#' @param ... currently not used
#' @return data.frame
#' @seealso \code{\link{chiSq_pearson}}
#' @references \insertRef{cameron2013regression}{Countr}
#' @export
chiSq_gof <- function(object, breaks, ...) {
    UseMethod("chiSq_gof", object)
}

#' @rdname chiSq_gof
#' @method chiSq_gof renewal
#' @export
chiSq_gof.renewal <- function (object, breaks, ...) {
    ## compute the score vector si
    si <- numDeriv::jacobian(object$score_fct, coef(object))
    colnames(si) <- paste0("si_", names(coef(object)))

    ## prepare the probabilities
    res <- chiSq_pearson(object)

    ## compute the d_{ij}(y_i) - p_{ij}(x_i, \theta), for j in 1, ..., J -1
    ## --- compute pij
    data_tmp <- object$data
    lhs <- as.character(object$formula)[2]

    .prob <- function(ind, counts) {
        dat <- data_tmp[ind, , drop = FALSE]

        .pred <- function(count) {
            dat[[lhs]] <- count
            predict(object, newdata = dat, type = "prob")$values
        }
        sapply(counts, .pred)
    }

    pij <- t(sapply(1:nrow(data_tmp), .prob, counts = res$Counts))
    colnames(pij) <- res$Counts

    if (!missing(breaks))
        pij <- .adjust_breaks(breaks, res, pij)
    ##---------------- end pij

    ## --- compute dij
    dij <- .get_dij_mat(data_tmp[[lhs]], colnames(pij))
    ## ------------ end dij

    dij_pij <- dij - pij
    colnames(dij_pij) <- paste0("di_", colnames(pij))

    ## build the regression matrix: 1_i ~ s_i + dij_pij, j =1,'' J-1
    out <- .run_chisq_reg(si, dij_pij, colnames(pij))

    attr(out, "pij") <- pij
    attr(out, "Pearson") <- res
    out
}

#' @rdname chiSq_gof
#' @method chiSq_gof negbin
#' @export
chiSq_gof.negbin <- function (object, breaks, ...) {
    ## compute the score vector si
    modelMat <- model.matrix(object)
    theta <- object$theta

    .obj <- function(th) {
        names(th) <- names(coef(object))
        mu <- object$family$linkinv(modelMat %*% th)
        dnbinom(object$y, mu = mu, size = theta, log = TRUE)
    }

    si <- numDeriv::jacobian(.obj, coef(object))
    colnames(si) <- paste0("si_", names(coef(object)))

    ## prepare the probabilities
    res <- chiSq_pearson(object)

    ## compute the d_{ij}(y_i) - p_{ij}(x_i, \theta), for j in 1, ..., J -1
    ## --- compute pij
    pij <- pscl::predprob(object)
    if (!missing(breaks))
        pij <- .adjust_breaks(breaks, res, pij)

     ## --- compute dij
    dij <- .get_dij_mat(object$y, colnames(pij))
    ## ------------ end dij

    dij_pij <- dij - pij
    colnames(dij_pij) <- paste0("di_", colnames(pij))

    ## build the regression matrix: 1_i ~ s_i + dij_pij, j =1,'' J-1
    out <- .run_chisq_reg(si, dij_pij, colnames(pij))

    attr(out, "pij") <- pij
    attr(out, "Pearson") <- res
    out
}

#' @rdname chiSq_gof
#' @method chiSq_gof glm
#' @export
chiSq_gof.glm <- function (object, breaks, ...) {
    ## compute the score vector si
    modelMat <- model.matrix(object)

    .obj <- function(th) {
        names(th) <- names(coef(object))
        lambda <- object$family$linkinv(modelMat %*% th)
        dpois(object$y, lambda = lambda, log = TRUE)
    }

    si <- numDeriv::jacobian(.obj, coef(object))
    colnames(si) <- paste0("si_", names(coef(object)))

    ## prepare the probabilities
    res <- chiSq_pearson(object)

    ## compute the d_{ij}(y_i) - p_{ij}(x_i, \theta), for j in 1, ..., J -1
    ## --- compute pij
    pij <- pscl::predprob(object)
    if (!missing(breaks))
        pij <- .adjust_breaks(breaks, res, pij)

     ## --- compute dij
    dij <- .get_dij_mat(object$y, colnames(pij))
    ## ------------ end dij

    dij_pij <- dij - pij
    colnames(dij_pij) <- paste0("di_", colnames(pij))

    ## build the regression matrix: 1_i ~ s_i + dij_pij, j =1,'' J-1
    out <- .run_chisq_reg(si, dij_pij, colnames(pij))

    attr(out, "pij") <- pij
    attr(out, "Pearson") <- res
    out
}


