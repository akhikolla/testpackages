#' Get named vector of coefficients for renewal objects
#'
#' Get named vector of coefficients for renewal objects.
#'
#' This is a convenience function for constructing named vector of coefficients
#' for renewal count models. Such vectors are needed, for example, for starting
#' values in the model fitting procedures. The simplest way to get a suitably
#' named vector is to take the coefficients of a fitted model but if the fitting
#' procedure requires initial values, this is seemingly a circular situation.
#'
#' The overall idea is to take coefficients specified by \code{object} and
#' transform them to coefficients suitable for a renewal count model as
#' specified by the arguments \code{"..."}. The provided methods eliminate the
#' need for tedius manual preparation of such vectors and in the most common
#' cases allow the user to do this in a single line.
#'
#' The default method extracts the coefficients of \code{object} using
#'
#' \code{co <- coef(object)} (an error is raised if this fails). It prepares a
#' named numeric vector with names requested by the arguments in \code{"..."}
#' and assigns \code{co} to the first \code{length(co)} elements of the prepared
#' vector. The net effect is that the coefficients of a model can be initialised
#' from the coefficients of a nested model. For example a Poisson regression
#' model can be used to initialise a Weibull count model. Of course the non-zero
#' shape parameter(s) of the Weibull model need to be set separately.
#'
#' If \code{object} is from class \code{glm}, the method is identical to the
#' default method.
#'
#' If object is from class \code{\link{renewalCoefList}}, its elements are
#' simply concatenated in one long vector.
#'
#' @param object an object, there are methods for several classes, see Details.
#' @param ... further arguments to be passed to \code{renewalNames}, usually
#'     something like \code{target = "weibull"}.
#' @seealso renewalNames
#'
#' @references
#' \insertRef{CountrJssArticle}{Countr}
#' 
#' @export
renewalCoef <- function(object, ...) {
    UseMethod("renewalCoef", object)
}

#' @export
renewalCoef.glm <- function(object, ...) {
    co <- coef(object)

    nams <- renewalNames(object, ...)

    res <- numeric(length(nams))
    res[1:length(co)] <- co

    names(res) <- nams
    res
}

#' @export
renewalCoef.default <- function(object, ...) {
    co <- coef(object)

    nams <- renewalNames(object, ...)

    res <- numeric(length(nams))
    res[1:length(co)] <- co

    names(res) <- nams
    res
}

#' @export
renewalCoef.renewalCoefList <- function(object, ...) {
    names(object) <- NULL
    unlist(object)
}


#' Split a vector using the prefixes of the names for grouping
#'
#' Split a vector using the prefixes of the names for grouping.
#'
#' The names of the coefficients of renewal regression models are prefixed with
#' the names of the parameters to which they refer. This function splits such
#' vectors into a list with one component for each parameter. For example, for a
#' Weibull renewal regression model this will create a list with components
#' \code{"scale"} and \code{"shape"}.
#'
#' This is a convenience function allowing users to manipulate the coefficients
#' related to a parameter more easily. \code{\link{renewalCoef}} can convert
#' this list back to a vector.
#'
#' @param coef a named vector
#' @export
renewalCoefList <- function(coef) {
    m <- regexpr("^[^_]+_", names(coef))
    nams <- unique(regmatches(names(coef), m))

    f <- function(nam) {
            ind <- grep(paste0("^", nam), names(coef))
            coef[ind]
         }

    res <- lapply(nams, f)
    names(res) <- nams
    class(res) <- "renewalCoefList"

    res
}

#' Get names of parameters of renewal regression models
#'
#' Get names of parameters of renewal regression models
#'
#' \code{renewalNames} gives the a character vector of names of parameters for
#' renewal regression models. There are two main use scenarios:
#' \code{renewalNames(object, target = "dist")} and
#' \code{renewalNames(object,...)}. In the first scenario \code{target} can be a
#' count distribution, such as "weibull" or a parameter name, such as shape.  In
#' this case \code{renewalNames} transforms coefficient names of \code{object}
#' to those specified by \code{target}. In the second cenario the argument list
#' is the same that would be used to call \code{renewalCount}. In this case
#' \code{renewalNames} returns the names that would be used by renewalCount for
#' the coefficients of the fitted model.
#'
#' @param object an object.
#' @param ...    further arguments.
#' @export
renewalNames <- function(object, ...) {
    UseMethod("renewalNames", object)
}

#' @export
renewalNames.default <- function(object, ...) {
    nams <- names(coef(object))
    renewalNames.character(nams, ...)
}


.old_prefix_replace <- function(object, oldpat, newpat) {
    if(length(newpat) == 1) {
        gsub(oldpat, newpat, object)
    } else {
        res <- lapply(newpat, function(x) gsub(oldpat, x, object))
        c(unlist(res))
    }
}

.put_prefix <- function(object, prefix) {
    if(length(prefix) == 1) {
        paste0(prefix, object)
    } else {
        res <- lapply(prefix, function(x) paste0(x, object))
        c(unlist(res))
    }
}



#' @export
renewalNames.character <- function(object, target = "scale", old.target = NULL, ...) {
    nams <- object

    if(length(target) == 1  &&  target %in% names(.distDescr))
        target <- .distDescr[[target]]$parNames

    newpat <- paste0(target, "_")

    if(is.character(old.target)) {
        if(length(old.target) > 1)
            stop("old.target has length > 1.")
        oldpat <- paste0("^", old.target, "_")

        nams <- .old_prefix_replace(nams, oldpat, newpat)
    } else {
        nams <- .put_prefix(nams, newpat) # paste0(target, "_", nams)
    }
    nams <- gsub('\\(Intercept\\)', "", nams)
    nams
}


#' @export
renewalNames.formula <- function(object, data, subset, na.action, weights, offset,
                         dist = c("weibull",
                                  "weibullgam", "custom", ##"burr"
                                  "gamma", "gengamma"),
                         anc = NULL, convPars = NULL, link = NULL, time = 1.0,
                         control = renewal.control(...), customPars = NULL,
                         seriesPars = NULL, weiMethod = NULL,
                         computeHessian = TRUE,
                         standardise = FALSE, standardise_scale = 1,
                         model = TRUE, y = TRUE, x = FALSE, ...) {
    formula <- object
    dist <- match.arg(dist)
    ## check convolution parameters
    convPars <- renewal.convPars(convPars, dist)

    if (dist == "custom")
        customPars <- .checkcustomPars(customPars, convPars$extrap)
    else if (dist == "weibull") {
        seriesPars <- renewal.seriesPars(seriesPars)
        weiMethod <- renewal.weiMethod(weiMethod)
    } else if (dist == "weibullgam") {
        warning(
            "weibullgam should be used with care! no guarantee of convergence !")
        anc <- NULL ## no regression allowed on aux pars
        seriesPars <- renewal.seriesPars(seriesPars, TRUE)
        weigamMethod <- weiMethod
        weigamMethod <- ifelse(is.null(weigamMethod), "series_acc",
                               weigamMethod)
        if (! weigamMethod %in% c("series_acc", "series_mat")) {
            warning(paste(weiMethod,
                          "is not an accepted method for weibullgam dist!",
                          "accelerated series will be used !"))
            weigamMethod <- "series_acc"
        }
        weiMethod <- weigamMethod
    }

    ## prepare the formula setting
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "na.action"), names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE

    f <- Formula(formula)
    mf[[1]] <- as.name("model.frame")
    mf$formula <- f
    mf <- eval(mf, parent.frame())
    ## Not sure about this: copied from hurdle: CHECK
    cl <- match.call()

    Y <- model.response(mf)
    n <- length(Y)
    ## stop if a formula with multiple response is passed
    if (is.null(Y))
        stop("muti-response formula not accepted !")
    ## convert negative reponse to zeros
    if (length(Y) < 1)
        stop("empty model")
    if (!isTRUE(all.equal(as.vector(Y), as.integer(round(Y + 0.001)))))
        warning(paste("invalid dependent variable,",
                      "non-integer values!",
                      "will be transformed")
                )
    Y <- as.integer(round(Y + 0.001))
    if (any(Y < 0))
        stop("invalid dependent variable, negative counts")

    ## extract weights and reshape them
    weights <- model.weights(mf)
    if (is.null(weights))
        weights <- 1
    if (length(weights) == 1)
        weights <- rep.int(weights, n)
    weights <- as.vector(weights)
    names(weights) <- rownames(mf)


    ## get model matrices
    modelMatrixList <- .getModelMatrix(formula = f, dist = dist, mf = mf,
                                       anc = anc, raw_data = data,
                                       standardise = standardise,
                                       standardise_scale = standardise_scale,
                                       customPars)



    ## check initilas values
    start <- control$start
    start <- .checkInitialValues(dist, start, modelMatrixList, weights, Y,
                                 customPars)
    nmPars <- gsub('\\(Intercept\\)', "", names(start))
    ## run optimization routine

    nmPars
}
