#' Predicting based on a \code{rocTree} model.
#'
#' The function gives predicted values with a \code{rocTree} fit.
#'
#' @param object is an \code{rocTree} object.
#' @param newdata is an optional data frame in which to look for variables with which to predict.
#' If omitted, the fitted predictors are used.
#' If the covariate observation time is not supplied, covariates will be treated as at baseline.
#' @param type is an optional character string specifying whether to predict the survival probability or the cumulative hazard rate.
#' @param control a list of control parameters. See 'details' for important special
#' features of control parameters. See \code{\link{rocTree}} for more details.
#' @param ... for future developments.
#'
#' @return Returns a \code{data.frame} of the predicted survival probabilities or cumulative hazard. 
#'
#' @importFrom stats model.frame
#' @export
#' @example inst/examples/ex_predict_rocTree.R
predict.rocTree <- function(object, newdata, type = c("survival", "hazard"),
                            control = list(), ...) {
    type <- match.arg(type)
    if (!is.rocTree(object)) stop("Response must be a 'rocTree' object")
    if (missing(newdata)) stop("Argument 'newdata' is missing")
    if (!(object$rName %in% colnames(newdata)))
        stop(paste("Object '", object$rName, "' not found.", sep = ""))
    if (!all(object$vName %in% colnames(newdata))) {
        missingName <- which(!(object$vName %in% colnames(newdata)))
        if (length(missingName) == 1)
            stop(paste("Object '", object$vName[missingName], "' not found.", sep = ""))
        if (length(missingName) > 1)
            stop(paste("Objects '", paste(object$vName[missingName], collapse = ", "),
                       "' not found.", sep = ""))
    }    
    type <- match.arg(type)
    control0 <- object$control
    control0[names(control0) %in% names(control)] <- control[names(control) %in% names(control0)]
    control <- control0
    raw <- newdata[findInt(object$data$.Y0, unlist(newdata[object$rName])), object$vNames]
    rownames(raw) <- NULL
    cutoff <- (1:control$nc) / (control$nc + 1)
    if (type %in% "survival") {
        if (object$ensemble)
            pred <- predict_rocForest_C(t(raw), object$data$.Y0, object$data$.D0, object,
                                        object$data$.X, object$disc, cutoff)
        else
            pred <- predict_rocTree_C(t(raw), object$data$.Y0, object$data$.D0, object,
                                      object$data$.X, object$disc, cutoff)
        object$survFun <- stepfun(object$data$.Y0, c(1, pred))
        object$pred <- data.frame(Time = unlist(newdata[,object$rName]),
                                  Survival = object$survFun(unlist(newdata[,object$rName])))
    }
    if (type %in% "hazard") {
        ## t0 <- seq(min(object$data$.Y0), max(object$data$.Y0), length.out = control$K)
        ## t0 <- seq(quantile(object$data$.Y0, .05), quantile(object$data$.Y0, .95),
        ##           length.out = control$K)
        t0 <- unlist(newdata[,object$rName])
        t0 <- seq(quantile(t0, .05), quantile(t0, .95), length.out = control$K)
        knots <- findInt(t0, object$data$.Y0)
        ## .mat1f2 <- sapply(object$data$.Y0[knots], K2, vec = object$data$.Y0, h = control$h) /
        .mat1f2 <- sapply(t0, K2, vec = object$data$.Y0, h = control$h) / control$h
        .mat1f2 <- .mat1f2[object$data$.D0 == 1,]
        if (object$ensemble)
            pred <- predict_rocForestHZ_C(t(raw[knots,]), t0, 
                                          object$data$.Y0, object$data$.D0, .mat1f2,
                                          control$h, object, object$data$.X,
                                          object$disc, cutoff)
        else
            pred <- predict_rocTreeHZ_C(t(raw[knots,]), t0, 
                                          object$data$.Y0, object$data$.D0, .mat1f2,
                                          control$h, object, object$data$.X,
                                          object$disc, cutoff)
        object$hazFun <- stepfun(t0, c(1, pred))
        object$pred <- data.frame(Time = t0, hazard = object$hazFun(t0))
    }
    rownames(object$pred) <- NULL
    class(object) <- "predict.rocTree"
    return(object)
}

is.predict.rocTree <- function(x) inherits(x, "predict.rocTree")

#' findInterval with 0 replaced with 1
#' @keywords internal
#' @noRd
findInt <- function(x, y) {
    pmax(1, findInterval(x, sort(y)))
}

#' findInterval with 0 replaced with 1, works with NA's in y
#' @keywords internal
#' @noRd
findInt.X <- function(x, y) {
    order(c(0, y))[pmax(1, findInterval(x, sort(c(0, y))))]
    ## order(y)[pmax(1, findInterval(x, sort(y)))]
}
