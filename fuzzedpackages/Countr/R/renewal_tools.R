#' Plot a frequency chart
#'
#' Plot a frequency chart to compare actual and predicted values.
#'
#' In order to compare actual and fitted values, a barchart plot is created.  It
#' is the user's responsibility to provide the count, observed and fitted
#' values.
#'
#' @param count_labels character, labels to be used.
#' @param actual numeric, the observed probabilities for the different count
#'     specified in \code{count_labels}.
#' @param pred data.frame of predicted values. Should have the same number of
#'     rows as actual and one column per model. The columns' names will be used
#'     as labels for the different models.
#' @param colours character vector of colour codes with \code{length}
#'     \code{ncol(pred)} + 2.
#' @export
frequency_plot <- function(count_labels, actual, pred, colours) {
    df <- data.frame(count_range = rep(count_labels, each = ncol(pred) + 1),
                     groups = rep(c("Actual", colnames(pred)),
                         length(count_labels)),
                     preds = as.numeric(t(as.matrix(cbind(actual, pred))))
                     )

    ## generate colours
    nc <- ncol(pred) + 2
    if (missing(colours))
        colours <- RColorBrewer::brewer.pal(nc, "Blues")
    else if (length(colours) < nc) {
        ns <- nc - length(colours)
        colours <- c(colours, RColorBrewer::brewer.pal(ns, "Blues"))
    } else
        colours <- colours[1:nc]

    settings <- list(
        superpose.polygon = list(col = colours[1:(nc -1)],
                                 border="transparent"),
        strip.background = list(col = colours[nc]),
        strip.border = list(col = "black")
    )

    lattice::barchart(preds ~ count_range, data = df, groups = groups,
                      ylab = "probability",
                      scales=list(alternating = 1),
                      auto.key = list(columns = ncol(pred) + 1,
                                      space = "top", points = FALSE,
                                      rectangles = TRUE),
                      par.settings = settings,
                      par.strip.text = list(col = "white", font = 2)
                      )
}

#' Compare renewals fit to glm models fit
#'
#' Compare renewals fit to glm models fit on the same data.
#'
#' This function computes a data.frame similar to Table 5.6 in Cameron(2013),
#' using the observed frequencies and predictions from different
#' models. Supported models accepted are Poisson and negative binomial (fitted
#' using \code{MASS::glm.nb()}) from the glm family and any model from the
#' renewal family (passed in \code{...}).
#'
#' @param poisson_model fitted Poisson glm model
#' @param nbinom_model fitted negative binomial (fitted using
#'     \code{MASS::glm.nb()}). This argument is optional.
#' @param breaks integer values at which the breaks should happen. The function
#'     will compute the observed frequencies in the intervals \code{[breaks[i],
#'     breaks[i + 1])}.
#' @param ... renewal models to be considered.
#' @return data.frame with columns \code{Counts}, \code{Actual} (observed
#'     probability) and then 2 columns per model passed (predicted probability
#'     and pearson statistic) for the associated count value.
#'
#' @references \insertRef{cameron2013regression}{Countr}
#' @export
compareToGLM <- function(poisson_model, breaks, nbinom_model, ...) {
    models <- list(...)
    stopifnot(length(models) > 0)
    model_names <- names(models)
    mms <- c("poisson", model_names)
    data <- models[[1]]$data
    n <- models[[1]]$n

    .comp <- function(obj, model_name) {
        stopifnot(class(obj)[1] == "renewal")
        tp <- chiSq_pearson(obj)[, c("Counts", "Actual", "Predicted", "Pearson")]
        names(tp)[3:4] <- c(paste0(model_name, "_predicted"),
                            paste0(model_name, "_pearson")
                            )
        tp
    }

    .pears <- function(p, obs) {
        n * (p - obs)^2 / p
    }

    ## prepare first model
    res <- .comp(models[[1]], model_names[1])

    ## add the Poisson predictions
    res$poisson_predicted <-
        sapply(res$Counts,
               function(x) {
                   mean(
                       suppressWarnings(
                           dpois(x, predict(poisson_model, newdata = data,
                                            type = "response")
                                 )
                       ),
                       na.rm = TRUE)
               })
    res$poisson_pearson <- .pears(res$poisson_predicted, res$Actual)

    if (!missing(nbinom_model)) {
        ## add the negative binomial model
        res$nbinom_predicted <-
            sapply(res$Counts,
                   function(x) {
                       mean(
                           suppressWarnings(
                               dnbinom(x,
                                       mu = predict(nbinom_model,
                                                    newdata = data,
                                                    type = "response"),
                                       size = nbinom_model$theta
                                       )
                           ),
                           na.rm = TRUE)
                   })
        res$nbinom_pearson <- .pears(res$nbinom_predicted, res$Actual)
    }

    ## add results from other models
    if (length(models) > 1) {
        for (i in 2:length(models)) {
            tmp <- .comp(models[[i]], model_names[i])
            res <- left_join(res, tmp, by = c("Counts", "Actual"))
        }
    }

    if (!missing(breaks)) {
        counts <- character()
        actual <- numeric()
        ## identify model columns
        ind_models <- which(grepl("_predicted", colnames(res)))
        models <- pears <- matrix(NA, nrow = length(breaks),
                                  ncol = length(ind_models))

        colnames(models) <- colnames(res)[ind_models]

        for (i in 1:(length(breaks) - 1)) {
            counts <- c(counts, .getNames_(breaks[i], breaks[i + 1]))
            actual <- c(actual, .getValues2_(breaks[i], breaks[i + 1],
                                             res$Counts, res$Actual)
                        )
            for (j in 1:length(ind_models))
                models[i, j] <- .getValues2_(breaks[i], breaks[i + 1],
                                             res$Counts, res[, ind_models[j]])
        }

        ## produce the final counts
        counts <- c(counts, paste(">=", breaks[length(breaks)]))
        counts <- factor(counts, levels = counts) # to keep ">= X" last

        actual <- c(actual, .getValues3_(breaks[length(breaks)], res$Counts,
                                         res$Actual)
                    )

        i <- length(breaks)
        for (j in 1:length(ind_models)) {
            ## models[i, j] <- .getValues3_(breaks[i],
            ##                              res$Counts, res[, ind_models[j]])
            models[i, j] <- 1 - sum(models[1:(i - 1), j])
            pears[, j] <- .pears(models[, j], actual)
        }

        colnames(pears) <- gsub("predicted", "pearson", colnames(models))

        res <- data.frame(Counts = counts, Actual = actual, models, pears)
    }

    dplyr::select(res, Counts, Actual,
                  contains("_predicted"), contains("_pearson")
                  )
}

.distDescr <- list(
    "weibull" = list(
        parNames         = c("scale", "shape"),
        linkDefaultNames = c(scale = "log", shape = "log")
        ),

    "weibullgam"  = list(
        parNames = c("scale", "shape", "shapeGam", "scaleGam"),
        linkDefaultNames = c(scale = "log", shape = "log",
            shapeGam = "identity", scaleGam = "identity")
        ),

    "gamma" = list(
        parNames         = c("rate", "shape"),
        linkDefaultNames = c(rate = "log", shape = "log")
        ),

    "gengamma" = list(
        parNames         = c("mu", "sigma", "Q"),
        linkDefaultNames = c(mu = "identity", sigma = "log", Q = "identity")
        ),

    "burr" = list(
        parNames         = c("scale", "shape1", "shape2"),
        linkDefaultNames = c(scale = "log", shape1 = "log", shape2 = "log")
        )
    )

#' Return the names of distribution parameters
#'
#' Return the names of the parameters of a count distribution.
#'
#' @param dist character, name of the distribution.
#' @param ... parameters to pass when dist == "custom".
#' @return character vector with the names of the distribution parameters.
#' @export
getParNames <- function(dist, ...) {
    switch(dist,
           "custom" = {
               customPars <- list(...)[[1]]
               customPars[["parNames"]]
           },
           ## default
           .distDescr[[dist]]$parNames
           )
}

.getDefaultLinkNames <- function(dist, parNames) {
    switch(dist,
           "custom" = {
               linkDefaultNames <- rep("identity", length(parNames))
               names(linkDefaultNames) <- parNames
               linkDefaultNames
           },
           ## default
           .distDescr[[dist]]$linkDefaultNames
           )
}

## inverse links - no need to compute them every time.
.inverseLinks <- list(
    "log"      = function(x) exp(x),
    "cauchit"  = function(x) cauchit(x, inverse = TRUE),
    "cloglog"  = function(x) cloglog(x, inverse = TRUE),
    "probit"   = function(x) probit(x, inverse = TRUE),
    "logit"    = function(x) logit(x, inverse = TRUE),
    "identity" = function(x) x
    )

## compute the link function inverse.
#' @importFrom VGAM cauchit
#' @importFrom VGAM cloglog
#' @importFrom VGAM probit
#' @importFrom VGAM logit
.computeInverseLink <- function(link = c("log", "cauchit", "cloglog",
                                    "probit", "logit", "identity")) {
    link <- match.arg(link)
    obj <- switch(link,
                  "log" = {function(x) exp(x)},
                  "cauchit" = {function(x) cauchit(x, inverse = TRUE)},
                  "cloglog" = {function(x) cloglog(x, inverse = TRUE)},
                  "probit" = {function(x) probit(x, inverse = TRUE)},
                  "logit" = {function(x) logit(x, inverse = TRUE)},
                  "identity" = {function(x) x}
                  )
    attr(obj, "functionName") <- link
    obj
}


## set attribute "functionName" for the functions in .inverseLinks
##  (note the non-local assignment)
lapply(names(.inverseLinks),
       function(link) attr(.inverseLinks[[link]], "functionName") <<- link
       )

## TODO: remove after testing and replace the calls to  .computeInverseLink() accordingly.
stopifnot(all(sapply(names(.inverseLinks),
    function(link)
        identical(body(.inverseLinks[[link]]), body(.computeInverseLink(link))) &&
        identical(args(.inverseLinks[[link]]), args(.computeInverseLink(link))) &&
        identical(attr(.inverseLinks[[link]], "functionName"),
                  attr(.computeInverseLink(link), "functionName"))
    )))

## create a string with link function information
.summarizeLinkInformation <- function(linkObj) {
    res <- character()
    for (i in seq(along = linkObj))
        res[i] <- paste0(names(linkObj)[i], ": link ",
                         attr(linkObj[[i]], "functionName"))
    paste0(res, collapse = ", ")
}

#' Creates the renewal control list
#'
#' Creates the renewal control list used by \code{renewal}.
#'
#' The function takes the user passed inputs, checks them (todo: actually, it
#' doesn't!) and returns an appropriate list that is used inside \code{renewal}
#' by the optimization routine, such as \code{optimx} among others.
#'
#' @param method character, one or more of the optimization methods accepted by
#'     \code{optimx}. User can experiment with different algorithms and
#' the results associated with the best performing one will be retained.
#' @param maxit numeric, the maximum number of iterations in the optimization
#'     routine.
#' @param trace Non-negative integer. Should tracing information be printed to
#'     the screen.
#' @param start (named) numeric, vector of starting values.
#' @param kkt locical should the Kuhn, Karush, Tucker optimality conditions be
#'     tested? Default is \code{FALSE} to avoid numerical hessian computation.
#' @param ... TODO
#'
#' @return a list with the control parameters used by \code{renewal}
#' @keywords internal
#' @export
renewal.control <- function(method = "nlminb", maxit = 1000, trace = 1,
                            start = NULL, kkt = FALSE, ...) {

    rval <- list(method = method, maxit = maxit, trace = trace, kkt = kkt,
                 start = start)
    rval <- c(rval, list(...))
    rval$maximize <- TRUE

    rval
}

#' Creates the convolution inputs setting
#'
#' Checks and creates the convolution inputs list.
#'
#' @param convPars a list of convolution parameters arguments with slots
#'     \code{nsteps}, \code{extrap} and \code{convMethod}, see
#'     \code{dCount_conv_bi}. If NULL, default parameters will be applied.
#' @param dist TODO
#'
#' @return list convolution inputs.
#' @keywords internal
#' @export
renewal.convPars <- function(convPars, dist) {
    extrap <- ifelse(!is.null(convPars$extrap), convPars$extrap,
                     ifelse(dist == "custom", FALSE, TRUE)
                     )
    nsteps <- ifelse(!is.null(convPars$nsteps),
                     convPars$nsteps,
                     ifelse(extrap, 50, 150))
    convMethod <- ifelse(!is.null(convPars$convMethod),
                         convPars$convMethod, "dePril")

    list(nsteps = nsteps, extrap = extrap, convMethod = convMethod)
}

## check custom-parameters
.checkcustomPars <- function(customPars, extrap) {
    ## check that par names exists
    if (is.null(customPars$parNames))
        stop("parNames should be provided in customPars !")
    else {
        if (length(customPars$parNames) < 1)
            stop("parNames must be at least of length 1 !")
        if (!is.character(customPars$parNames))
            stop("parNames must be a character string !")
    }

    ## check survivalFct
    if (is.null(customPars$survivalFct))
        stop("survivalFct should be provided in customPars !")
    else if (class(customPars$survivalFct) != "function")
        stop("survivalFct must be a function object !")

    if (extrap) {
        ## check extrapolFct
        if (is.null(customPars$extrapolFct))
            stop("extrapolFct should be provided in customPars !")
        else if (class(customPars$extrapolFct) != "function")
            stop("extrapolFct must be a function object !")
    }

    customPars
}


#' Creates the series expansion inputs setting
#'
#' Check and creates the series expansion inputs list.
#'
#' @param seriesPars list series expansion input parameters with slots
#'     \code{terms} (number of terms in the series expansion), \code{iter}
#'     (number of iteration in the accelerated series expansion algorithm) and
#'     \code{eps} (tolerance in the accelerated series expansion algorithm),
#'     Only used if \code{dist = "weibull"} and \code{weiMethod =
#'     c("series_mat", "series_acc")}.
#' @param long TODO
#'
#' @return list series expansion inputs.
#' @keywords internal
#' @export
renewal.seriesPars <- function(seriesPars, long = FALSE) {
    terms <- ifelse(!is.null(seriesPars$terms), seriesPars$terms,
                    ifelse(long, 100, 50))
    iter <- ifelse(!is.null(seriesPars$iter), seriesPars$iter, 300)
    eps <- ifelse(!is.null(seriesPars$eps), seriesPars$eps, 1e-10)

    list(terms = terms, iter = iter, eps = eps)
}

#' Check weibull computation algorithm
#'
#' Check weibull computation algorithm.
#'
#' @param weiMethod character, desired  weibull method.
#' @return a valid weibull computation method.
#' @keywords internal
#' @export
renewal.weiMethod <- function(weiMethod) {
    if (is.null(weiMethod))
        weiMethod <- "series_acc"
    else if (!weiMethod %in% c("series_acc", "series_mat",
                               "conv_direct", "conv_naive", "conv_dePril")) {
        weiMethod <- "series_acc"
        warning(paste(weiMethod, "is not an accepted method for weibull dist!",
                      "accelerated series will be used !"))
    }

    weiMethod
}

## create a string with link function information
.summarizeLinkInformation <- function(linkObj) {
    res <- character()
    for (i in seq(along = linkObj))
        res[i] <- paste0(names(linkObj)[i], ": link ",
                         attr(linkObj[[i]], "functionName"))
    paste0(res, collapse = ", ")
}

## apply (inverse) link function to the standard errors of the parameters
.transformSE <- function(se, linkList) {
    nmPars <- names(linkList)
    ## i=1 always transform
    par1 <- nmPars[1]
    key1 <- paste0(par1, "_")
    ll1 <- which(grepl(key1, names(se)))
    se[ll1] <- linkList[[1]](se[ll1])

    for (i in 2:length(nmPars)) {
        pari <- nmPars[i]
        keyi <- paste0(pari, "_")
        lli <- which(grepl(keyi, names(se)))
        if (length(lli) > 1)
            se[lli] <- linkList[[i]](se[lli])
        else if (length(lli) == 1) {
            if (nchar(names(se)[lli]) > nchar(keyi)) ##
                se[lli] <- linkList[[i]](se[lli])
        }
    }
    se
}

.getDistParNames <- function(dist = c("weibull", "weibullgam",
                                 "gamma", "gengamma", "burr")) {
    dist <- match.arg(dist)

    switch(dist,
           "weibull" = {
               parNames <- c("scale", "shape")
           },
           "gamma" = {
               parNames <- c("rate", "shape")
           },
           "gengamma" = {
               parNames <- c("mu", "sigma", "Q")
           },
           "burr" = {
               parNames <- c("scale", "shape1", "shape2")
           },
           "weibullgam" = {
               parNames <- c("scale", "shape", "shapeGam", "scaleGam")
           }
           )

    parNames
}

.getLinkDefaultNames <- function(dist = c("weibull", "weibullgam",
                                     "gamma", "gengamma", "burr")) {
    dist <- match.arg(dist)
    switch(dist,
           "weibull" = {
               linkDefaultNames <- c(scale = "log", shape = "log")
           },
           "gamma" = {
               linkDefaultNames <- c(rate = "log", shape = "log")
           },
           "gengamma" = {
               linkDefaultNames <- c(mu = "identity", sigma = "log",
                                     Q = "identity")
           },
           "burr" = {
               linkDefaultNames <- c(scale = "log", shape1 = "log",
                                     shape2 = "log")
           },
           "weibullgam" = {
               linkDefaultNames <- c(scale = "log",
                                     shape = "identity",
                                     shapeGam = "identity",
                                     scaleGam = "identity")
           }
           )
    linkDefaultNames
}

## get the inverse link function for each distribution
.getLinkList <- function(dist, link, ...) {

    switch(dist,
           "custom" = {
               customPars <- list(...)[[1]]
               parNames <- customPars[["parNames"]]
               linkDefaultNames <- rep("identity", length(parNames))
               names(linkDefaultNames) <- parNames
           }, {
               parNames <- .getDistParNames(dist)
               linkDefaultNames <- .getLinkDefaultNames(dist)
           }
           )

    linkRes <- list()
    for (i in seq(along = parNames)) {
        pari <- parNames[i]
        linki <- link[[pari]]
        if(is.null(linki))
            linki <- .computeInverseLink(linkDefaultNames[pari])
        else
            linki <- .computeInverseLink(linki)

        linkRes[[pari]] <- linki
    }

    class(linkRes) <- "InverseLink"
    linkRes
}

## remove duplicated columns from model.matrix
.clean_model_matrix <- function(modelMat) {
    ind <- which(grepl("^`", colnames(modelMat)) &
                 grepl("`$", colnames(modelMat))
                 )
    dups_names <- colnames(modelMat)[ind]
    dups_names <- gsub("^`", "", dups_names)
    dups_names <- gsub("`$", "", dups_names)

    if (length(dups_names) > 0) {
        ind_remove <- numeric(0)
        for (i in 1:length(dups_names)) {
            namei <- dups_names[i]
            ii <- which(colnames(modelMat) == namei)
            if (all(modelMat[, ii] == modelMat[, ind[i]]))
                ind_remove <- c(ind_remove, ind[i])
        }

        return(modelMat[, -ind_remove, drop = FALSE])
    } else
        return(modelMat)
}

.extract_rhs <- function(formula) {
    tmp <- as.character(formula)
    if (length(tmp) == 3) return(tmp[2])
    else return(NA)
}

## model matrix returned as a list
.getModelMatrix <- function(formula, dist, mf, anc, raw_data,
                            standardise = FALSE, standardise_scale = 1, ...) {
    if (standardise) {
        tmp_formula <- terms.formula(formula, data = raw_data)
        tmp_standardise <-
            standardize::standardize(formula = formula(tmp_formula),
                                     data = raw_data,
                                     scale = standardise_scale,
                                     family = poisson)
        X <- .clean_model_matrix(
            model.matrix(tmp_standardise$formula,
                         data = tmp_standardise$data, rhs = 1)
            )
        attr(X, "standardize") <- tmp_standardise
    } else
        X <- .clean_model_matrix(model.matrix(formula, data = mf, rhs = 1))

    switch(dist,
           "custom" = {
               customPars <- list(...)[[1]]
               parNames <- customPars[["parNames"]]
           },
           {
               parNames <- .getDistParNames(dist)
           }
           )

    modelMatRes <- list()
    modelMatRes[[parNames[1]]] <- X
    if (length(parNames) > 1) {
        for (i in 2:length(parNames)) {
            pari <- parNames[i]
            if (pari %in% names(anc)) {
                fi <- Formula(anc[[pari]])
                if (standardise) {
                    tmp_formula <- terms.formula(fi, data = raw_data)
                    tmp_standardise <-
                        standardize::standardize(formula = formula(tmp_formula),
                                                 data = raw_data,
                                                 scale = standardise_scale,
                                                 family = poisson)
                    Xi <- .clean_model_matrix(
                        model.matrix(tmp_standardise$formula,
                                     data = tmp_standardise$data, rhs = 1)
                        )
                    attr(Xi, "standardize") <- tmp_standardise
                } else {
                    mfi <- model.frame(formula = fi, data = raw_data)
                    Xi <- .clean_model_matrix(model.matrix(fi, data = mfi, rhs = 1))
                }
                ## store results
                modelMatRes[[pari]] <- Xi
            }
        }
    }

    class(modelMatRes) <- c("modelMat", class(modelMatRes))
    modelMatRes
}

predict.modelMat <- function(object, newdata, ...) {
    .fct <- function(x) {
        stan <- attr(x, "standardize")
        response <- .extract_rhs(stan$formula)
        if (response %in% names(newdata))
            mf <- predict(stan, newdata, response = TRUE)
        else {
            mf <- predict(stan, newdata)
            mf[[response]] <- NA
        }

        .clean_model_matrix(model.matrix(stan$formula, data = mf, rhs = 1))
    }

    res <- lapply(object, .fct)
    names(res) <- names(object)
    class(res) <- c("modelMat", class(res))
    res
}

## return value is controlled by arguments Ev and summa:
##      if Ev = TRUE,  return the expected value of the count,
##      if Ev = FALSE and summa = TRUE, return loglikelihood,
##      if Ev = FALSE and summa = FALSE, return probabilities.
##
## 2017-08-03 - (a) new argument 'log' to pass on to probability computing
##                      functions (used only if ev = FALSE and summa = FALSE).
##              (b) made xMax an explicit argument rather than look for it
##                      using list(...).
## TODO: there are calls to .objectiveFunction() in the package which use unnamed
##       argument, clean up this (then log may go safely before '...').
.objectiveFunction <- function(params, dist, modelMatrixList, linkList,
                               time, convPars, Y = NULL, weights = NULL,
                               Ev = FALSE, summa = TRUE,
                               seriesPars = NULL, weiMethod = NULL,
                               ..., log = FALSE, xMax = NULL) {
    ## check if Y is provided when Ev = TRUE
    if (!Ev) {
        if(is.null(Y))
            stop("Y should be provided when Ev is FALSE!")
        else {
            ## weights
            if (is.null(weights))
                weights <- rep(1, length(Y))
        }
    }

    if (is.null(xMax)) {
        if (!is.null(Y))
            xMax <- max(max(Y), 15) ## 2017-08-02 magic number, 15, here and below.
        else {
            xMax <- 15
            print("xMax not provided! 15 will be used by default !")
        }
    }
    
    switch(dist,
           "custom" = {
               customPars <- list(...)[[1]]
               parNames <- customPars[["parNames"]]
               distParValues <- list()
               locationPar <- parNames[1]
               linkLocation <- linkList[[locationPar]]
               X <- modelMatrixList[[locationPar]]

               ## we assume regression for the first parameter (location)
               if(is.null(names(params)))
                   location <-
                       as.vector(linkLocation(X %*% params[1:ncol(X)]))
               else
                   location <-
                       as.vector(linkLocation(X %*%
                                              params[grepl(paste0(locationPar,
                                                                  "_"),
                                                           names(params))]
                                              )
                                 )

               distParValues[[locationPar]] <- location
               ix <- ncol(X)
               for (i in 2:length(parNames)) {
                   pari <- parNames[i]
                   linki <- linkList[[pari]]
                   if (pari %in% names(modelMatrixList)) {
                       Xi <- modelMatrixList[[pari]]
                       if (is.null(names(params)))
                           distParValues[[pari]] <-
                               as.vector(
                                   linki(Xi %*%
                                         params[(ix + 1):
                                                (ix + ncol(Xi))
                                                ]
                                         )
                                   )
                       else
                           distParValues[[pari]] <-
                               as.vector(
                                   linki(Xi %*%
                                         params[grepl(paste0(pari, "_"),
                                                      names(params))]
                                         )
                                   )

                       ix <- ix + ncol(Xi)
                   } else {
                       distParValues[[pari]] <-
                           linki(as.numeric(rep(params[ix + 1], nrow(X))))
                       ix <- ix + 1
                   }
               }

               ## seq along Y
               seky <- seq(along = Y)
               ## extraction tools
               survFct <- customPars[["survivalFct"]]
               ## create distPars in the appropriate format
               .fct0 <- function(xlist, i = 1) xlist[[i]]
               .fcti <- function(i) lapply(distParValues, .fct0, i = i)
               distPars <- lapply(seky, .fcti)

               if (convPars$extrap) {
                   .getextrapolPars <- customPars[["extrapolFct"]]
                   .fctExtrap <- function(i) .getextrapolPars(distPars[[i]])
               } else
                   .fctExtrap <- function(i) c(2, 2) ## not used

               extrapParsList <- lapply(seky, .fctExtrap)

               if (Ev) { ## expected value and variance
                   .evfct <- function(i)
                       evCount_conv_user(xmax = xMax, survR = survFct,
                                         distPars = distPars[[i]],
                                         extrapolPars = .fctExtrap(i),
                                         method = convPars$convMethod,
                                         nsteps = convPars$nsteps,
                                         time = time,
                                         extrap = convPars$extrap)

                   out <- lapply(seky, .evfct)
                   attr(out, "distPars") <- distParValues
                   return(out)

               } else { ## loglikelihood or probability value
                   if (summa)
                       return(dCount_conv_loglik_user(x = Y,
                                                      survR = survFct,
                                                      distPars = distPars,
                                                      extrapolPars =
                                                      extrapParsList,
                                                      nsteps = convPars$nsteps,
                                                      extrap = convPars$extrap,
                                                      method =
                                                      convPars$convMethod,
                                                      time = time,
                                                      weights = weights)
                              )
                   else {
                       return(
                           dCount_conv_user(x = Y, survR = survFct,
                                            distPars = distPars,
                                            extrapolPars = extrapParsList,
                                            nsteps = convPars$nsteps,
                                            extrap = convPars$extrap,
                                            method = convPars$convMethod,
                                            time = time, log = log)
                           )
                  }
               }
           },
           "weibull" = { ## we do weibull separately on purpose
               if (is.null(seriesPars))
                   seriesPars <- renewal.seriesPars(seriesPars)
               if (is.null(weiMethod))
                   weiMethod <- renewal.weiMethod(weiMethod)
               ## ------------- scale parameter -----------------------------
               linkScale <- linkList[["scale"]]
               X <- modelMatrixList[["scale"]]
               ## check parameters name
               if(is.null(names(params)))
                   scale <-
                       as.vector(linkScale(X %*% params[1:ncol(X)]))
               else
                   scale <-
                       as.vector(linkScale(X %*%
                                           params[grepl("scale_",
                                                        names(params))]
                                           )
                                 )

               ## ------------- shape parameter -----------------------------
               linkShape <- linkList[["shape"]]
               if ("shape" %in% names(modelMatrixList)) {
                   linkShape <- linkList[["shape"]]
                   Xshape <- modelMatrixList[["shape"]]
                   if(is.null(names(params)))
                       shape <- as.vector(
                           linkShape(Xshape %*%
                                     params[(ncol(X) + 1):length(params)]
                                     )
                           )
                   else
                       shape <- as.vector(
                           linkShape(Xshape %*%
                                     params[grepl("shape_", names(params))]
                                     )
                           )
               } else
                   shape <- linkShape(rep(params[length(params)], nrow(X)))

               if (Ev) {
                   .fct <- function(i)
                       evWeibullCount(xmax = xMax, shape = shape[i],
                                      scale = scale[i],
                                      method = weiMethod, time = time,
                                      conv_steps = convPars$nsteps,
                                      conv_extrap = convPars$extrap,
                                      series_terms = seriesPars$terms,
                                      series_acc_niter = seriesPars$iter,
                                      series_acc_eps = seriesPars$eps
                                      )

                   out <- lapply(1:nrow(X), .fct)
                   ## distParValues values
                   distParValues <- list()
                   distParValues$scale <- scale
                   distParValues$shape <- shape
                   attr(out, "distPars") <- distParValues
                   return(out)
               } else {
                   if (summa) {
                       return(
                           dWeibullCount_loglik(
                               x = Y, shape = shape, scale = scale,
                               conv_steps = convPars$nsteps,
                               conv_extrap = convPars$extrap,
                               method = weiMethod,
                               time = time, weights = weights,
                               series_terms = seriesPars$terms,
                               series_acc_niter = seriesPars$iter,
                               series_acc_eps = seriesPars$eps)
                              )
                   } else
                       return(
                           dWeibullCount(
                               x = Y, shape = shape, scale = scale,
                               conv_steps = convPars$nsteps,
                               conv_extrap = convPars$extrap,
                               method = weiMethod, time = time,
                               series_terms = seriesPars$terms,
                               series_acc_niter = seriesPars$iter,
                               series_acc_eps = seriesPars$eps,
                               log = log)
                          )
               }
           },
           "weibullgam" = { ## we do weibull-gamma separately as well
               ## ------------- scale parameter -----------------------------
               X <- modelMatrixList[["scale"]]
               ## check parameters name
               if(is.null(names(params))) {
                   scale <- as.vector(params[1:ncol(X)])
                   shape <- params[ncol(X) + 1]
                   shapeGam <- params[ncol(X) + 2]
                   scaleGam <- params[ncol(X) + 3]
               } else {
                   scale <- params[grepl("scale_", names(params))]
                   shape <- params["shape_"]
                   shapeGam <- params["shapeGam_"]
                   scaleGam <- params["scaleGam_"]
               }

               ## Note that we change parametrization here in terms of
               ## 1/r and 1/alpha

               if (Ev) {
                   .fct <- function(i) {
                       evWeibullgammaCount(xmax = xMax, shape = shape,
                                           shapeGam = 1.0 / shapeGam,
                                           scaleGam = 1.0 / scaleGam,
                                           Xcovar = X[i, ,drop = FALSE],
                                           beta = scale,
                                           method = weiMethod, time = time,
                                           series_terms = seriesPars$terms,
                                           series_acc_niter = seriesPars$iter,
                                           series_acc_eps = seriesPars$eps
                                      )
                   }
                   out <- lapply(1:nrow(X), .fct)
                   ## distParValues values
                   distParValues <- list()
                   distParValues$scale <- scale
                   distParValues$shape <- shape
                   distParValues$scaleGam <- scale
                   distParValues$shapeGam <- shape
                   attr(out, "distPars") <- distParValues
                   return(out)

               } else {
                   if (summa) {
                       return(
                           dWeibullgammaCount_loglik(
                               x = Y, shape = shape,
                               shapeGam = 1.0 / shapeGam, scaleGam = 1.0 / scaleGam,
                               Xcovar = X, beta = scale,
                               method = weiMethod,
                               time = time, weights = weights,
                               series_terms = seriesPars$terms,
                               series_acc_niter = seriesPars$iter,
                               series_acc_eps = seriesPars$eps)
                              )
                   } else
                       return(
                           dWeibullgammaCount(
                               x = Y, shape = shape,
                               shapeGam = 1.0 / shapeGam,
                               scaleGam = 1.0 / scaleGam,
                               Xcovar = X, beta = scale,
                               method = weiMethod, time = time,
                               series_terms = seriesPars$terms,
                               series_acc_niter = seriesPars$iter,
                               series_acc_eps = seriesPars$eps,
                               log = log)
                              )
               }
           },
           { ## default
               parNames <- .getDistParNames(dist)
               distParValues <- list()
               locationPar <- parNames[1]
               linkLocation <- linkList[[locationPar]]
               X <- modelMatrixList[[locationPar]]

               ## we assume regression for the first parameter (location)
               if(is.null(names(params)))
                   location <-
                       as.vector(linkLocation(X %*% params[1:ncol(X)]))
               else
                   location <-
                       as.vector(linkLocation(X %*%
                                              params[grepl(paste0(locationPar,
                                                                  "_"),
                                                           names(params))]
                                              )
                                 )

               distParValues[[locationPar]] <- location
               ix <- ncol(X)
               for (i in 2:length(parNames)) {
                   pari <- parNames[i]
                   linki <- linkList[[pari]]
                   if (pari %in% names(modelMatrixList)) {
                       Xi <- modelMatrixList[[pari]]
                       if(is.null(names(params)))
                           distParValues[[pari]] <-
                               as.vector(
                                   linki(Xi %*%
                                         params[(ix + 1):
                                                (ix + ncol(Xi))
                                                ]
                                         )
                                   )
                       else
                           distParValues[[pari]] <-
                               as.vector(
                                   linki(Xi %*%
                                         params[grepl(paste0(pari, "_"),
                                                      names(params))]
                                         )
                                   )

                       ix <- ix + ncol(Xi)
                   } else {
                       distParValues[[pari]] <-
                           linki(as.numeric(rep(params[ix + 1], nrow(X))))
                       ix <- ix + 1
                   }
               }

               ## seq along Y
               seky <- seq(along = Y)
               ## create distPars in the appropriate format

               .fct0 <- function(xlist, i = 1) xlist[[i]]
               .fcti <- function(i) lapply(distParValues, .fct0, i = i)
               distPars <- lapply(seky, .fcti)

               if (Ev) { ## expected value and variance
                   .evfct <- function(i)
                       evCount_conv_bi(xmax = xMax,
                                       distPars = distPars[[i]],
                                       dist = dist,
                                       method = convPars$convMethod,
                                       nsteps = convPars$nsteps,
                                       time = time,
                                       extrap = convPars$extrap)

                   out <- lapply(seky, .evfct)
                   attr(out, "distPars") <- distParValues
                   return(out)

               } else { ## loglikelihood or probability value
                   if (summa)
                       return(dCount_conv_loglik_bi(x = Y,
                                                    distPars = distPars,
                                                    dist = dist,
                                                    nsteps = convPars$nsteps,
                                                    extrap = convPars$extrap,
                                                    method =
                                                    convPars$convMethod,
                                                    time = time,
                                                    weights = weights)
                              )
                   else {
                       return(
                           dCount_conv_bi(x = Y,
                                          distPars = distPars,
                                          dist = dist,
                                          nsteps = convPars$nsteps,
                                          extrap = convPars$extrap,
                                          method = convPars$convMethod,
                                          time = time,
                                          log = log)
                           )
                   }
               }
           })
}

## extract apprpriate modelData that can be used later to construc se
.modelData <- function(modelMatrixList, dist, ...) {
    switch(dist,
           "custom" = {
               customPars <- list(...)[[1]]
               parNames <- customPars[["parNames"]]
           }, {
               parNames <- .getDistParNames(dist)
           }
           )

    locationPar <- parNames[1]
    X <- modelMatrixList[[locationPar]]
    nm <- paste0(locationPar, "_", colnames(X))
    colnames(X) <- gsub('\\(Intercept\\)', "", nm)

    if (length(parNames) > 1) {
        for (i in 2:length(parNames)) {
            pari <- parNames[i]
            if (pari %in% names(modelMatrixList)) {
                Xi <- modelMatrixList[[pari]]
                nm <- paste0(pari, "_", colnames(Xi))
                colnames(Xi) <- gsub('\\(Intercept\\)', "", nm)
            } else {
                Xi <- matrix(1, ncol = 1, nrow = nrow(X))
                colnames(Xi) <- paste0(pari, "_")
            }
            X <- cbind(X, Xi)
        }
    }
    X
}

.getPredictionStd <- function(modelMatrixList, vcov, dist, link, ...) {
    switch(dist,
           "custom" = {
               customPars <- list(...)[[1]]
               parNames <- customPars[["parNames"]]
           }, {
               parNames <- .getDistParNames(dist)
           }
           )

    seList <- list()
    ## do location parameter separately to avoid a call to if
    locationPar <- parNames[1]
    X <- modelMatrixList[[locationPar]]
    nm <- paste0(locationPar, "_", colnames(X))
    colnames(X) <- gsub('\\(Intercept\\)', "", nm)

    ## create standard error for the location parameter
    indX <- which(grepl(locationPar, colnames(X)))
    Xloc <- X[, indX]
    indCov <- which(grepl(locationPar, colnames(X)))
    covloc <- vcov[indCov, indCov]
    seloc <- sqrt(diag(Xloc %*% covloc %*% t(Xloc)))
    seList[[locationPar]] <- link[[locationPar]](seloc)

    ## anc parameters
    if (length(parNames) > 1) {
        for (i in 2:length(parNames)) {
            pari <- parNames[i]
            if (pari %in% names(modelMatrixList)) {
                Xi <- modelMatrixList[[pari]]
                nm <- paste0(pari, "_", colnames(Xi))
                colnames(Xi) <- gsub('\\(Intercept\\)', "", nm)
            } else {
                Xi <- matrix(1, ncol = 1, nrow = nrow(X))
                colnames(Xi) <- paste0(pari, "_")
            }
            indXi <- which(grepl(pari, colnames(Xi)))
            Xi <- Xi[, indXi, drop = FALSE]
            indCov <- which(grepl(pari, colnames(Xi)))
            covi <- vcov[indCov, indCov, drop = FALSE]
            sei <- sqrt(diag(Xi %*% covi %*% t(Xi)))
            seList[[pari]] <- link[[pari]](sei)
        }
    }
    seList
}


## temp function
.extractElem <- function(xList, ind = "ExpectedValue")
    as.numeric(xList[[ind]])

.checkHess <- function(hess, nPars) {
    ## 2019-08-01 was:
    ##     ans <- TRUE
    ##     
    ##     if (is.null(hess))
    ##         return(FALSE)
    ##     
    ##     if (any(is.na(hess)))
    ##         return(FALSE)
    ##     
    ##     if (!is.matrix(hess))
    ##         return(FALSE)
    ##     
    ##     if (ncol(hess) != nPars)
    ##         return(FALSE)
    ##     
    ##     if (nrow(hess) != nPars)
    ##         return(FALSE)
    ##     
    ##     return(ans)
    
    is.matrix(hess)  &&  ncol(hess) == nPars  && nrow(hess) == nPars  &&
        all(!is.na(hess))
}

## replace NA/NAN that appears in log-likelihood by logMin due to very small
## probability.
#' @keywords internal
.logNaReplace <- function() {
    log(.Machine$double.xmin)
}
