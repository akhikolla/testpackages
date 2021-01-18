#' @include renewal_tools.R
NULL

#' Deprecated old name of \code{renewalCount()}
#'
#' Deprecated old name of \code{renewalCount()}
#'
#' This function is now renamed to \code{renewalCount()}.
#'
#' @param ... same arguments as \code{renewalCount()}.
#'
#' @keywords internal
#' @export
renewal <- function(...) {
    .Deprecated(msg = "Please use 'renewalCount' (the new name of 'renewal')")
    renewalCount(...)
}

.process_formula <- function(formula, dist){
    Fo <- Formula(formula)
    Fo.len <- length(Fo)
    if(Fo.len[2] < 2) {
        if(Fo.len[1] == 1)
            return(NULL)
        formula <- formula(Fo, lhs = 1)
        if(Fo.len[2] == 0)
            return(list(formula = formula, anc = NULL))

        anc <- lapply(2:Fo.len[1], function(k) formula(Fo, lhs = 0))
        lhsnames <-  sapply(2:Fo.len[1],
                            function(k) as.character(formula(Fo, lhs = k, rhs = 0))[2] ) #  .[1] is '~'
        names(anc) <- lhsnames
    } else {
        formula <- formula(Fo, lhs = 1, rhs = 1)
        anc <- lapply(2:Fo.len[2], function(k) formula(Fo, lhs = 0, rhs = k))

        lhsnames <- as.character(formula(Fo, rhs = 0))[2] # .[1] is '~'
        lhsnames <- lhsnames[-1] # drop the first name
        if(length(lhsnames) < length(anc)){
            if(length(lhsnames) == 0) {
                lhsnames <- getParNames(dist)[-1]
            } else {
                distnames <- getParNames(dist) # todo: 'custom' needs special treatment
                distnames <- distnames[! distnames %in% lhsnames]
                lhsnames <- c(lhsnames, distnames[1:(length(lhsnames) - 1)])
            }
        }
        stopifnot(length(lhsnames) == length(anc))

        names(anc) <- lhsnames
    }

    list(formula = formula, anc = anc)
}

#' Fit renewal count processes regression models
#'
#' Fit renewal regression models for count data via maximum likelihood.
#'
#' \code{renewal} re-uses design and functionality of the basic R tools for
#' fitting regression model (\code{lm}, \code{glm}) and is highly inspired by
#' \code{hurdle()} and \code{zeroinfl()} from package \code{pscl}. Package
#' \code{Formula} is used to handle formulas.
#'
#' Argument \code{formula} is a \code{formula} object. In the simplest case its
#' left-hand side (lhs) designates the response variable and the right-hand side
#' the covariates for the first parameter of the distribution (as reported by
#' \code{\link{getParNames}}. In this case, covariates for the ancilliary
#' parameters are specified using argument \code{anc}.
#'
#' The ancilliary regressions, can also be specified in argument \code{formula}
#' by adding them to the righ-hand side, separated by the operator \sQuote{|}.
#' For example \code{Y | shape ~ x + y | z} can be used in place of the pair
#' \code{Y ~ x + y} and \code{anc = list(shape = ~z)}. In most cases, the name
#' of the second parameter can be omitted, which for this example gives the
#' equivalent \code{Y ~ x + y | z}. The actual rule is that if the parameter is
#' missing from the left-hand side, it is inferred from the default parameter
#' list of the distribution.
#'
#' As another convenience, if the parameters are to to have the same covariates,
#' it is not necessary to repeat the rhs. For example, \code{Y | shape ~ x + y}
#' is equivalent to \code{Y | shape ~ x + y | x + y}. Note that this is applied
#' only to parameters listed on the lhs, so \code{Y ~ x + y} specifies
#' covariates only for the response variable and not any other parameters.
#'
#'
#' Distributions for inter-arrival times supported internally by this package
#' can be chosen by setting argument \code{"dist"} to a suitable character
#' string.  Currently the built-in distributions are \code{"weibull"},
#' \code{"weibullgam"}, \code{"gamma"}, \code{"gengamma"} (generalized-gamma)
#' and \code{"burr"}.
#'
#' Users can also provide their own inter-arrival distribution.  This is done by
#' setting argument \code{"dist"} to \code{"custom"}, specifying the initial
#' values and giving argument \code{customPars} as a list with the following
#' components:
#'
#' \describe{
#'
#' \item{parNames}{character, the names of the parameters of the distribution.
#'     The location parameter should be the first one.}
#'
#' \item{survivalFct}{function object containing the survival function. It
#'     should have signature \code{function(t, distPars)} where \code{t} is the
#'     point where the survival function is evaluated and \code{distPars} is the
#'     list of the distribution parameters. It should return a double value.}
#'
#' \item{extrapolFct}{function object computing the extrapolation values
#'     (numeric of length 2) from the value of the distribution parameters (in
#'     \code{distPars}). It should have signature \code{function(distPars)} and
#'     return a numeric vector of length 2. Only required if the extrapolation
#'     is set to \code{TRUE} in \code{convPars}.}
#' }
#'
#' Some checks are done to validate \code{customPars} but it is user's
#' responsibility to make sure the the functions have the appropriate
#' signatures.
#'
#' \strong{Note:} The Weibull-gamma distribution is an experimental version and
#' should be used with care! It is very sensitive to initial values and there is no
#' guarantee of convergence. It has also been reparameterized in terms of
#' \eqn{(1/r, 1/\alpha, c)}{(1/r, 1/alpha, c)} instead of \eqn{(r, \alpha,
#' c)}{(r, alpha, c)}, where \eqn{r} and \eqn{\alpha}{alpha} are the shape
#' and scale of the gamma distribution and \eqn{c} is the shape of the Weibull
#' distribution.
#'
#' \strong{(2017-08-04(Georgi) experimental feature:} probability residuals in
#' component 'probResiduals'. I also added type 'prob' to the method for
#' residuals() to extract them.
#'
#' probResiduals[i] is currently 1 - Prob(Y[i] given the
#' covariates). "one minus", so that values close to zero are "good".  On its
#' own this is probably not very useful but when comparing two models, if one of
#' them has mostly smaller values than the other, there is some reason to claim
#' that the former is superior. For example (see below), gamModel < poisModel in
#' 3:1
#'
#' @param formula a formula object. If it is a standard formula object, the left
#'     hand side specifies the response variable and the right hand sides
#'     specifies the regression equation for the first parameter of the
#'     conditional distribution. \code{formula} can also be used to specify the
#'     ancilliary regressions, using the operator `|`, see Details.
#'
#' @param data,subset,na.action, arguments controlling formula processing via
#'     \code{model.frame}.
#' @param weights optional numeric vector of weights.
#' @param offset optional numeric vector with an a priori known component to be
#'     included in the linear predictor of the count model. Currently not used.
#' @param dist character, built-in distribution to be used as the inter-arrival
#'     time distribution or \code{"custom"} for a user defined distribution, see
#'     Details. Currently the built-in distributions are \code{"weibull"},
#'     \code{"weibullgam"}, \code{"gamma"}, \code{"gengamma"}
#'     (generalized-gamma) and \code{"burr"}.
#' @param anc a named list of formulas for ancillary regressions, if any,
#'     otherwise \code{NULL}. The formulas associated with the (exact) parameter
#'     names are used. The left-hand sides of the formulas in \code{anc} are
#'     ignored.
#' @param link named list of character strings specifying the name of the link
#'     functions to be used in the regression. If \code{NULL}, the canonical
#'     link function will be used, i.e, \code{log} if the parameter is supposed
#'     to be positive, identity otherwise.
#' @param time numeric, time at which the count is observed; default to unity
#'     (1).
#' @param convPars a list of convolution parameters arguments with slots
#'     \code{nsteps}, \code{extrap} and \code{convMethod}, see
#'     \code{dCount_conv_bi}. If NULL, default parameters will be applied.
#' @param control a list of control arguments specified via
#'     \code{renewal.control}.
#' @param customPars list, user inputs if \code{dist = "custom"}, see details.
#' @param seriesPars list, series expansion input parameters with slots
#'     \code{terms} (number of terms in the series expansion), \code{iter}
#'     (number of iteration in the accelerated series expansion algorithm) and
#'     \code{eps} (tolerance in the accelerated series expansion algorithm),
#'     Only used if \code{dist = "weibull"} and \code{weiMethod =
#'     c("series_mat", "series_acc")}.
#' @param weiMethod character, computation method to be used if \code{dist =
#'     "weibull"} or \code{"weibullgam"}, see \code{dWeibullCount} and
#'     \code{dWeibullgammaCount}.
#' @param computeHessian logical, should the hessian (and hence the covariance
#'     matrix) be computed numerically at the fitted values.
#' @param standardise logical should the covariates be standardised using
#'     \code{standardize::standardize()} function.
#' @param standardise_scale numeric the desired scale for the covariates;
#'     default to 1
#' @param model,y,x logicals. If \code{TRUE} the corresponding components of the
#'     fit (model frame, response, model matrix) are returned.
#' @param ... arguments passed to \code{renewal.control} in the default setup.
#'
#' @return An \code{S3} object of class "renewal", which is a list with
#'     components including:
#'     \describe{
#'
#'     \item{coefficients}{values of the fitted coefficients.}
#'
#'     \item{residuals}{vector of weighted residuals \eqn{\omega * (observed -
#'     fitted)}.}
#'
#'     \item{fitted.values}{vector of fitted means.}
#'
#'     \item{optim}{data.frame output of \code{optimx}.}
#'
#'     \item{method}{optimisation algorithm.}
#'
#'     \item{control}{the control arguments, passed to \code{optimx}.}
#'
#'     \item{start}{starting values, passed to \code{optimx}.}
#'
#'     \item{weights}{weights to apply, if any.}
#'
#'     \item{n}{number of observations (with weights > 0).}
#'
#'     \item{iterations}{number of iterations in the optimisation algorithm.}
#'
#'     \item{execTime}{duration of the optimisation.}
#'
#'     \item{loglik}{log-likelihood of the fitted model.}
#'
#'     \item{df.residual}{residuals' degrees of freedom for the fitted model.}
#'
#'     \item{vcoc}{convariance matrix of all coefficients, computed numerically
#'     from the hessian at the fitted coefficients (if \code{computeHessian} is
#'     \code{TRUE}).}
#'
#'     \item{dist}{name of the inter-arrival distribution.}
#'
#'     \item{link}{list, inverse link function corresponding to each parameter in
#'     the inter-arrival distribution.}
#'
#'     \item{converged}{logical, did the optimisation algorithm converge?}
#'
#'     \item{data}{data used to fit the model.}
#'
#'     \item{formula}{the original formula.}
#'
#'     \item{call}{the original function call.}
#'
#'     \item{anc}{named list of formulas to model regression on ancillary
#'     parameters.}
#'
#'     \item{score_fct}{function to compute the vector of scores defined in
#'     Cameron(2013) equation 2.94.}
#'
#'     \item{convPars}{convolution inputs used.}
#'
#'     \item{customPars}{named list, user passed distribution inputs, see
#'     Details.}
#'
#'     \item{time}{observed window used, default is 1.0 (see inputs).}
#'
#'     \item{model}{the full model frame (if \code{model = TRUE}).}
#'
#'     \item{y}{the response count vector (if \code{y = TRUE}).}
#'
#'     \item{x}{the model matrix (if \code{x = TRUE}).}
#' }
#'
#' @examples
#' \dontrun{
#' ## may take some time to run depending on your CPU
#' data(football)
#' wei = renewalCount(formula = homeTeamGoals ~ 1,
#'                     data = football, dist = "weibull", weiMethod = "series_acc",
#'                     computeHessian = FALSE, control = renewal.control(trace = 0, 
#'                     method = "nlminb"))
#' }
#' @importFrom numDeriv hessian jacobian
#' @importFrom MASS ginv
#' @import optimx
#' @references
#'
#' \insertRef{CountrJssArticle}{Countr}
#' 
#' \insertRef{cameron2013regression}{Countr}
#'
#' @export
renewalCount <- function(formula, data, subset, na.action, weights, offset,
                         dist = c("weibull",
                                  "weibullgam", "custom", ##"burr"
                                  "gamma", "gengamma"),
                         anc = NULL, convPars = NULL, link = NULL, time = 1.0,
                         control = renewal.control(...), customPars = NULL,
                         seriesPars = NULL, weiMethod = NULL,
                         computeHessian = TRUE,
                         standardise = FALSE, standardise_scale = 1,
                         model = TRUE, y = TRUE, x = FALSE, ...) {
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

    wrk <- .process_formula(formula, dist)
    if(!is.null(wrk)){ # ancilliary spec's are in formula not 'anc'
        if(!is.null(anc))
            stop("ancilliary parameters should be in 'formula' or 'anc' but not in both")

        ## !!! overwrite 'formula' and 'anc'
        formula <- wrk$formula
        anc <- wrk$anc
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
# browser()
    ## extract weights and reshape them
    weights <- model.weights(mf)
    if (is.null(weights))
        weights <- 1
    if (length(weights) == 1)
        weights <- rep.int(weights, n)
    weights <- as.vector(weights)
    names(weights) <- rownames(mf)

    ## get the inverse link function
    if (class(link) == "InverseLink")
        linkList <- link
    else
        linkList <- .getLinkList(dist, link, customPars)

    ## get model matrices
    modelMatrixList <- .getModelMatrix(formula = f, dist = dist, mf = mf,
                                       anc = anc, raw_data = data,
                                       standardise = standardise,
                                       standardise_scale = standardise_scale,
                                       customPars)


    ## create objective function
    ## ============== countDist to pass to optim ===============================
    countDist <- function(params) {
        .objectiveFunction(params, dist, modelMatrixList,
                           linkList, time, convPars, Y, weights,
                           Ev = FALSE, summa = TRUE, seriesPars, weiMethod,
                           customPars)
    }

    ## check initilas values
    start <- control$start
    start <- .checkInitialValues(dist, start, modelMatrixList, weights, Y,
                                 customPars)
    nmPars <- gsub('\\(Intercept\\)', "", names(start))
    ## run optimization routine

    method <- control$method
    hessian <- ifelse(is.null(control$hessian), FALSE, control$hessian)
    control$method <- control$start <- control$hessian <- NULL
    ## remove warning from optimx
    control$dowarn <- FALSE

    if (control$trace) {
        print("calling optimx() for parameter estimation by ML ...")
        fitCount_ <- optimx(par = start, fn = countDist, method = method,
                            hessian = hessian, control = control)
    } else { ## remove the annoying message: Maximizing -- use negfn and neggr
        sink(tempfile())
        fitCount_ <-
            optimx(par = start, fn = countDist, method = method,
                   hessian = hessian, control = control)
        sink()
    }

    ## arrange results by best (largest) value and then fastest
    ## 2019-06-13 patch for issue#1; was:
    ##     fitCount_ <- fitCount_[with(fitCount_, order(-value, xtimes)), ]
    ##
    ##     now sort on 'xtimes' only if it is present.
    ##
    ## TODO: Note that column 'xtimes' is not missing on Windows, it is just called 'xtime'.
    ##       (Need to talk to the maintainer of optimx.)
    ##
    ## TODO: it would be better to get rid of the `with()` here...
    wrk2 <- if("xtimes" %in% names(fitCount_))
                with(fitCount_, order(-value, xtimes))
            else if("xtime" %in% names(fitCount_))
                with(fitCount_, order(-value, xtime))
            else
                with(fitCount_, order(-value))

    fitCount_ <- fitCount_[wrk2, ]

    fitCount <- fitCount_[1, , drop = FALSE]
    class(fitCount) <- class(fitCount_)
    ## coefficients
    coefs <- as.numeric(coef(fitCount))
    names(coefs) <- nmPars


    ## variance-covariance matrix
    if (computeHessian) {
        hess <- try(attr(fitCount, "details")[method, "nhatend"][[1]],
                    silent = TRUE)
        if (!.checkHess(hess, length(coefs))) {
            if (control$trace)
                print("computing a numerical approximation to the Hessian ...")
            hess <- numDeriv::hessian(countDist, coefs)
        }

        varCovarcount <- try(-solve(hess))
        if ((inherits(varCovarcount, "try-error")) |
            (any(diag(varCovarcount) < 0)) ) {
            varCovarcount <- Matrix::nearPD(-ginv(hess))$mat
            warning(paste("variance-covariance matrix was computed",
                          "by smoothing the genralized inverse hessian !"))
        }

        dimnames(varCovarcount) <- list(nmPars, nmPars)
    } else
        varCovarcount <- matrix()

    .obj_jac <- function(theta) {
        .objectiveFunction(theta, dist, modelMatrixList,
                           linkList, time, convPars, Y, weights,
                           Ev = FALSE, summa = FALSE, seriesPars, weiMethod,
                           customPars, log = TRUE)
    }

    ## residuals (Pearson)
    resTemp <-  .objectiveFunction(coefs, dist, modelMatrixList,
                                   linkList, time, convPars, Y, weights,
                                   TRUE, FALSE, seriesPars, weiMethod,
                                   customPars)

    ## 2017-08-03(Georgi) new - probabilities
    obsProbs <- .objectiveFunction(coefs, dist, modelMatrixList,
                                   linkList, time, convPars, Y, weights,
                                   FALSE, FALSE, seriesPars, weiMethod,
                                   customPars)


    vecDistParsList <- attr(resTemp, "distPars")
    Yhat <- sapply(resTemp, .extractElem, ind = "ExpectedValue")
    wi <- sapply(resTemp, .extractElem, ind = "Variance")
    res <- sqrt(weights) * (Y - Yhat)

    ## number of observations
    nobs <- sum(weights > 0)

    ## collect best method and number of iterations from fitCount_
    method <- rownames(fitCount_)[1]
    niter <- fitCount_$niter[1]

    ## value to be returned
    rval <- list(
        coefficients = coefs, residuals = res, fitted.values = Yhat, wi = wi,
        optim = fitCount_, method = method, control = control, start = start,
        vecDistParsList = vecDistParsList,
        weights =
        if (identical(as.vector(weights), rep.int(1L, n))) NULL else weights,
        n = nobs, iterations = niter, execTime = sum(fitCount_$xtimes),
        loglik = fitCount$value, df.residual = nobs - length(coefs),
        vcov = varCovarcount, score_fct = .obj_jac, dist = dist,
        link = linkList, converged = fitCount$convcode[1] == 0, data = data,
        formula = formula(f), call = cl, anc = anc, convPars = convPars,
        customPars = customPars, time = time, seriesPars = seriesPars,
        weiMethod = weiMethod,
        ## standardize option
        standardise = standardise, standardise_scale = standardise_scale,
        ## 2017-08-03(Georgi) (provisional) new component:
        ##    TODO: needs documenting if it is here to stay.
        probResiduals = 1 - obsProbs # subtracting from one to make zero the reference line.
        )
    if (model)
        rval$model <- mf
    if (y)
        rval$y <- Y

    ## save model matrix if required or standardise option is on
    if (x | standardise)
        rval$x <- modelMatrixList


    class(rval) <- "renewal"
    return(rval)
}
