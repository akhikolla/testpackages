#' Do a D*M analysis
#'
#' @param formula A formula object of the form:
#' \code{binary response ~ reaction time + condition1 * condition2 * ... conditionN}.
#' @param data A dataframe for looking up data specified in formula.
#' For backwards compatibility this can also be with: a column named \code{rt} containing response times in ms,
#' a column named \code{response} containing at most 2 response options, and an
#' optional column named \code{condition} containing a numeric index as to which conditions
#' observations belong.
#' @param tt A time grid on which the density function will be evaluated.
#' Should be larger than the highest observed reaction time.
#' @param restr A restriction matrix where each column depicts one condition.
#' The number of rows should match the number of parameters (and be equal to the length of lower).
#' The contents of \code{restr} should be numbers, identical numbers means that these parameters
#' (either within or between condition) will be constrained. Different numbers means parameters
#' will not be constrained.
#' @param fixed A matrix that allows for fixing parameters to certain values.
#' @param lower Should be a vector containing lower bounds for each parameter.
#' Has a default if \code{fun.density == Voss.density}.
#' @param upper Should be a vector containing upper bounds for each parameter.
#' Has a default if \code{fun.density == Voss.density}.
#' @param DstarM If TRUE a D*M analysis is done, otherwise the Chi square distance
#' between data and model is minimized.
#' @param SE positive value, how many standard error to add to the variance to relax
#' the variance restriction a bit.
#' @param oscPdf Logical, if TRUE check for oscillations in calculated densities and
#' remove densities with oscillations.
#' @param Optim a named list with identical arguments to \code{\link{DEoptim.control}}.
#' In addition, if \code{verbose} == TRUE \code{Optim$steptol} can be a vector, i.e.
#' \code{c(200, 50, 10)} means: Do 200 iterations then check for convergence, do 50
#' iterations then check for convergence, check every 10 iterations for convergence until
#' itermax is reached. Defaults to \code{Optim = list(reltol = 1e-6, itermax = 1e3,
#' steptol = 50, CR = .9, trace = 0, parallelType = 0)}.
#' @param splits Numeric vector determining which conditions have an equal nondecision density.
#' Identical values in two positions indicate that the conditions corresponding to the indices
#' of those values have an identical nondecision distribution.
#' @param forceRestriction if TRUE the variance restriction is enforced.
#' @param mg Supply a data density, useful if a uniform kernel approximation does not suffice.
#' Take care that densities of response categories within conditions are degenerate and therefore integrate to the proportion a category was observed (and not to 1).
#' @param h bandwidth of a uniform kernel used to generate data based densities.
#' @param pars Optional parameter vector to supply if one wishes to evaluate the objective
#' function in a given parameter vector. Only used if \code{itermax} equal zero.
#' @param fun.density Function used to calculate densities. See details.
#' @param args.density A names list containing additional arguments to be send to fun.density.
#' @param fun.dist Function used to calculate distances between densities.
#' Defaults to a chi-square distance.
#' @param args.dist A named list containing additional arguments to be send to fun.dist.
#' @param verbose Numeric, should intermediate output be printed? Defaults to 1, higher values result in more progress output.
#' Estimation will speed up if set to 0. If set to TRUE, \code{Optim$trace} will be
#' forced to 0, hereby disabling the build in printing of \code{DEoptim}. To enable the
#' printing of \code{DEoptim}, set \code{verbose} to 0 and specify \code{Optim$trace}.
#' \code{Optim}. If set to 1, ETA refers to the expected maximum time until completion (when the iterations limit is reached).
#' @param useRcpp Logical, setting this to true will make the objective function use an Rcpp implementation
#' of \code{Voss.density} with the distance function \code{chisq}. This gains speed at the cost of flexibility.
#'
#' @details Response options will be alphabetically sorted and the first response option will be
#' treated as the 'lower' option. This means that if the observed proportion of the first
#' response options is higher, the drift speed will most likely be negative.
#'
#' \code{fun.density} allows a user to specify a custom density function. This function must (at least) take the following arguments:
#' \code{t}: a vector specifying at which time points to calculate the density
#' \code{pars}: a parameter vector
#' \code{boundary}: character 'upper' or 'lower' specifying for which response option the density will be calculated.
#' \code{DstarM}: Logical, if TRUE the density should not describe the nondecision density, if FALSE it should describe the nondecision density.
#' Any additional arguments can be passed to \code{fun.density} via the argument \code{args.density}.
#' If one intends to use a custom density function it is recommended to test the function first with \code{\link{testFun}}.
#' When specifying a custom density function it is probably also necessary to change the lower and upper bounds of the parameter space.
#'
#' For purposes of speed, the function can be run in parallel by providing the argument \code{Optim = list(parallelType = 1)}.
#' See \code{\link{DEoptim.control}} for details. Also, for Ratcliff models the objective function has been rewritten in Rcpp.
#' This limits some functionality but does result in a faster estimation. Usage of Rcpp can be enabled via \code{useRcpp = TRUE}.
#'
#' When verbose is set to 1, the ETA is an estimated of the time it takes to execute ALL iterations.
#' Convergence can (and is usually) reached before then.
#'
#' @return Returns a list of class \code{DstarM.fitD} that contains:
#' \item{Bestvals}{Named numeric vector. Contains the best parameter estimates.}
#' \item{fixed}{Numeric vector. Contains the best parameter estimates.}
#' \item{GlobalOptimizer}{List. All output from the call to \code{\link{DEoptim}}}
#' \item{Debug}{List. contains the number of DEoptim iterations, the number of function evaluation of the objective function, and the maximum number of iterations.}
#' \item{note}{String. A possible note that is used for summary purposes}
#' \item{tt}{Numeric vector. Contains the time grid used.}
#' \item{g.hat}{Numeric matrix. Named columns represent the (possibly smoothed) densities of the data distribution of each condition-response pair.}
#' \item{modelDist}{Numeric matrix. Named columns represent the densities of the model evaluated at grid \code{tt} with parameters \code{Bestvals}.}
#' \item{ncondition}{Numeric scalar. The number of conditions}
#' \item{var.data}{Numeric vector. The variance of each condition-response pair. There are as many values as hypothesized nondecision densities.}
#' \item{var.m}{Numeric vector. The variance of the model distributions. There are as many values as hypothesized nondecision densities.}
#' \item{restr.mat}{Numeric matrix. Contains the restrictions used.}
#' \item{splits}{Numeric vector. Equal to the input argument with the same name.}
#' \item{n}{Numeric scalar. The total number of observations.}
#' \item{DstarM}{Logical. Equal to the input argument with the same name.}
#' \item{fun.density}{Function. Equal to the input argument with the same name.}
#' \item{fun.dist}{Function. Equal to the input argument with the same name.}
#' \item{h}{Scalar. Equal to the input argument with the same name.}
#' \item{args.density}{Named list. Equal to the input argument with the same name.}
#' \item{args.dist}{Named list. Equal to the input argument with the same name.}
#'
#' @examples
#'
#'# simulate data with three stimuli of different difficulty.
#'# this implies different drift rates across conditions.
#'# define a time grid. A more reasonable stepsize is .01; this is just for speed.
#'tt = seq(0, 5, .1)
#'pars = c(.8, 2, .5, .5, .5, # condition 1
#'         .8, 3, .5, .5, .5,  # condition 2
#'         .8, 4, .5, .5, .5)  # condition 3
#'pdfND = dbeta(tt, 10, 30)
#'# simulate data
#'data = simData(n = 3e3, pars = pars, tt = tt, pdfND = pdfND)
#'# define restriction matrix
#'restr = matrix(1:5, 5, 3)
#'restr[2, 2:3] = 6:7 # allow drift rates to differ
#'# fix variance parameters
#'fixed = matrix(c('sz1', .5, 'sv1', .5), 2, 2)
#'\dontrun{
#'# Run D*M analysis
#'res = estDstarM(data = data, tt = tt, restr = restr, fixed = fixed)
#'coef(res)
#'summary(res)
#'}

#' @useDynLib DstarM, .registration = TRUE

#' @export
estDstarM <- function(formula = NULL, data, tt, restr = NULL, fixed = list(),
    lower, upper, Optim = list(), DstarM = TRUE, SE = 0, oscPdf = TRUE, splits = rep(0L,
        (ncondition)), forceRestriction = TRUE, mg = NULL, h = 1, pars, fun.density = Voss.density,
    args.density = list(), fun.dist = chisq, args.dist = list(tt = tt), verbose = 1L,
    useRcpp = TRUE) {

    # Error handling
    verbose <- as.integer(verbose)
    if (!is.numeric(verbose))
        stop("verbose must be numeric")
    if (verbose < 0L) {
        verbose <- 0L
    } else if (verbose > 2L) {
        verbose <- 2L
    }

    Optim <- errCheckOptim(Optim)
    by <- unique(zapsmall(diff(tt)))

    data <- getData(formula, data, verbose = verbose)
    rtime <- data[["rtime"]]
    response <- data[["response"]]
    condition <- data[["condition"]]
    hasConditions <- data[["hasConditions"]]
    data <- data[["data"]]

    note <- errCheckData(data = data, tt = tt, h = h, by = by, rtime = rtime,
        response = response, condition = condition)
    ncondition <- length(unique(data[[condition]]))  # get number of conditions

    # mm is a helper matrix. matrix multiplication with this matrix sums
    # every 2 columns often pdfs needs to be summed over conditions, i.e.
    # when normalizing them or calculating their variance.
    mm <- matrix(0, ncondition * 2, ncondition)
    mm[1:dim(mm)[1L] + dim(mm)[1L] * rep(1:dim(mm)[2L] - 1, each = 2)] <- 1

    # get defaults
    if (missing(lower) | missing(upper)) {
        if (DstarM) {
            parnames <- c("a", "v", "z", "sz", "sv")
            lower <- c(0.01, -6, 0.05, 0.01, 0)
            upper <- c(2, 6, 0.95, 0.99, 10)
        } else {
            parnames <- c("a", "v", "t0", "z", "sz", "sv", "st0")
            lower <- c(0.01, -6, 0.1, 0.05, 0.01, 0, 0)
            upper <- c(2, 6, 0.8, 0.95, 0.99, 10, 1)
        }
    } else {
        stopifnot(length(lower) == length(upper), all(lower < upper))
        if (is.null(names(lower))) {
            parnames <- paste0("par", seq_along(lower), "_")
        } else {
            parnames <- names(lower)
        }
    }
    npars <- length(lower)

    if (is.null(restr)) {
        restr.mat <- matrix(1:(npars * ncondition), npars, ncondition)
    } else {
        # convert numbers to 1:length(unique)
        restr.mat <- apply(as.matrix(c(restr)), 1, function(x, eq) which(x ==
            eq), unique(c(restr)))
        dim(restr.mat) <- dim(as.matrix(restr))
        if (dim(restr.mat)[2L] != ncondition)
            stop(sprintf("Number of columns of restr (%s) must match number of conditions (%s).",
                dim(restr.mat)[2L], ncondition), call. = FALSE)
        if (dim(restr.mat)[1L] != npars)
            stop(sprintf("Number of rows of restr must match number of parameters",
                dim(restr.mat)[1L], npars), call. = FALSE)
        uniqR <- unique(c(restr.mat))
        for (i in seq_along(uniqR)) {
            idxR <- which(restr.mat == uniqR[i], arr.ind = TRUE)
            if (!all(duplicated(idxR[, 1])[-1])) {
                idx <- !duplicated(idxR[, 1])
                warning(sprintf("Different parameters have been restricted to the same value (see restr %s). This is unorthodox but the analysis will continue.",
                  paste0("[", idxR[idx, 1], ", ", idxR[idx, 2], "]", collapse = "; ")),
                  immediate. = TRUE, call. = FALSE)
            }
        }
    }
    # rcpp checks
    if (useRcpp && !(identical(fun.density, Voss.density) || identical(fun.dist,
        chisq))) {
        stop("useRcpp only implemented for fun.density = Voss.density and fun.dist = chisq.")
    }

    # Indexing based on restrictions on nondecision distribution and number
    # of conditions
    group <- groups(nconditie = ncondition, splits = splits, restriction = FALSE)
    # mm2 contains which pdfs should be summed for calculating variances
    mm2 <- matrix(0, 2 * ncondition, dim(group)[2L])
    for (i in 1:dim(group)[2L]) {
        mm2[group[, i], i] <- 1
    }
    rt <- split(data[[rtime]], list(data[[response]], data[[condition]]))
    n <- (ql <- lengths(rt)) %*% mm2

    if (is.null(mg)) {
        # obtain data density
        g <- getGhat(rt = rt, tt, ncondition, mm, by)
    } else {
        # if a model or data density is supplied
        errCheckDatamg(mg, tt, ncondition)
        g <- mg
        # n = Inf # irrelevent number
        SE <- 0
        # ql = rep.int(1, dim(g)[2L])
    }
    # calculate first 5 moments - done this way to maybe add kurtosis
    # restrictions in the future
    g2 <- g %*% mm2
    g2 <- g2 %*% (diag(dim(g2)[2L])/apply(g2, 2, simpson, x = tt))
    moments.data <- matrix(0, 5, dim(g2)[2L])
    for (i in 1:dim(g2)[2L]) {
        moments.data[, i] <- unlist(lapply(0:4, nth.cmomentS, x = tt, fx = g2[,
            i]))
    }
    # convolve data density with uniform kernel - only do this when no own
    # density is supplied perhaps make it possible to also smooth? although
    # users can also do this themselves since they are already supplying
    # their own densities.
    if (DstarM && is.null(mg)) {
        kernel <- rev(by * stats::dunif(tt, 0, h))
        for (i in 1:dim(g)[2L]) {
            g[, i] <- customConvolveO(as.double(g[, i]), kernel)[seq_along(tt)]
        }
        g[g < 0] <- 0
    }
    g <- g %*% (diag(dim(g)[2L])/rep(apply(g %*% mm, 2, simpson, x = tt),
        each = 2))
    colnames(g) <- names(rt)  # name columns of g after conditions for clarity
    mu.data <- moments.data[1, ]
    var.data <- moments.data[3, ]
    # apply SE - gives the variance restriction some lenience
    var.data <- var.data * (1 + SE * sqrt(2/(n - 1)))  # TODO: include ref

    # ii and jj are vectors that contain all combinations of p and p' to
    # compare in the objective function
    iijj <- getCombinations(group)
    ii <- iijj[["ii"]]
    jj <- iijj[["jj"]]

    # Impose parameter constraints restrictions between conditions are
    # imposed via a numerical matrix fixing of parameters is done by looking
    # up the names of the vectors.  if both lower and upper missing test if
    # function actually works
    if (!all(missing(lower) & missing(upper))) {
        args <- list(t = tt, pars = lower, boundary = "lower", DstarM = DstarM)  # , moreArgs
        testFun(fun.density = fun.density, args = args, lower = lower, upper = upper)
    }
    # replicate by number of conditions
    lower <- rep(lower, ncondition)
    upper <- rep(upper, ncondition)
    names(lower) <- paste0(parnames, rep(seq_len(ncondition), each = npars))
    names(upper) <- names(lower)

    # search condition restrictions and implement these in the lower/ upper
    # vector
    search <- unique(c(restr.mat))
    ind <- rep.int(0L, (length(search)))
    for (i in search) {
        ind[i] <- (which(restr.mat == i, arr.ind = TRUE)[1, ] + c(0, -1)) %*%
            c(1, npars)
    }
    lower <- lower[ind]
    upper <- upper[ind]

    # get restrictions
    fixed <- getFixed(fixed, names(lower), useRcpp)
    indFixed <- fixed$indFixed
    if (!is.null(indFixed)) {
        lower <- lower[-indFixed]
        upper <- upper[-indFixed]
    }

    # keep only arguments in args.density that appear in fun.density
    args.density <- args.density[names(args.density) %in% names(formals(fun.density))]

    # collect arguments so far in a list
    if (useRcpp) {
        argsList <- list(tt = tt, restr = restr.mat - 1, ii = ii - 1L, jj = jj -
            1L, mm = mm, mm2 = mm2, g = g, varData = var.data, ql = ql, DstarM = DstarM,
            oscPdf = oscPdf, forceRestriction = forceRestriction, precision = 3)
    } else {
        # make every columnn of restr.mat temporarily a list
        restrList <- unlist(apply(restr.mat, 2, list), recursive = FALSE,
            use.names = FALSE)
        argsList <- list(tt = tt, g = g, ql = ql, restr.mat = restrList, DstarM = DstarM,
            oscPdf = oscPdf, ii = ii, jj = jj, mm = mm, mm2 = mm2, fun.density = fun.density,
            args.density = args.density, fun.dist = fun.dist, args.dist = args.dist,
            var.data = var.data, forceRestriction = forceRestriction, by = by,
            parnames = names(lower))
    }
    if (missing(pars)) {
        # change to is is.null(pars)?

        argsList$fixed <- fixed

        # do differential evolution
        nrep <- ceiling((Optim$itermax - sum(Optim$steptol[1:(length(Optim$steptol) -
            1)]))/Optim$steptol[length(Optim$steptol)])
        if (nrep <= 1L) {
            nrep <- which.max(cumsum(Optim$steptol) * (cumsum(Optim$steptol) <=
                Optim$steptol))
        } else if (nrep > 1L) {
            last <- length(Optim$steptol)
            Optim$steptol2 <- rep.int(0, (length(Optim$steptol) + nrep))
            Optim$steptol2[1L:(last)] <- Optim$steptol
            Optim$steptol2[(last + 1L):length(Optim$steptol2)] <- Optim$steptol[last]
            Optim$steptol <- Optim$steptol2
            nrep <- length(Optim$steptol)
        }

        oldval <- Inf
        nfeval <- 0  # number of function evaluations
        # add arguments for DEoptim to the list
        if (useRcpp) {

            argsList$fn <- totalobjectiveC
            argsList$fixed <- fixed$mat
            argsList$anyFixed <- fixed$anyFixed

            # imposeFixations(argsList$fixed, lower, names(lower))
            # imposeFixationsC(pars = lower, fixed = qq) browser() if (length(fixed)
            # == 0) { argsList$anyFixed = FALSE argsList$fixed = matrix(1, 1, 1) #
            # c++ does not accept NULL } else { argsList$anyFixed = TRUE
            # argsList$fixed = fixed2Rcpp(argsList$fixed, lower) }

        } else {
            argsList$all <- FALSE
            argsList$fn <- total.objective

        }
        argsList$lower <- lower
        argsList$upper <- upper

        argsList$control <- DEoptim::DEoptim.control()
        sharedNames <- names(argsList$control)[(names(argsList$control) %in%
            names(Optim))]
        argsList$control[sharedNames] <- Optim[sharedNames]

        if (Optim$parallelType == 1) {
            argsList$control$packages <- c(list("rtdists", "DstarM"), Optim$packages)
            argsList$control$parVar <- Optim$parVar
        } else if (Optim$parallelType == 2) {
            argsList$control$foreachArgs <- c(list("getPdf", "getVar", "oscCheck",
                "simpson", "nth.cmomentS", "customDdiffusion", "customApprox",
                "customConvolveO"), Optim$foreachArgs)
        }
        if (verbose > 0) {
            cat(sprintf("Starting %s analysis...\n", ifelse(DstarM, "D*M",
                "Traditional")))
            argsList$control$itermax <- Optim$steptol[1]
            argsList$control$storepopfrom <- Optim$steptol[1] - 1
            argsList$control$trace <- 0
            tempOut <- list(Bestvals = rep(NA, length(unique(unlist(restr.mat)))),
                restr.mat = restr.mat)
            class(tempOut) <- "DstarM"
            prevSize <- 0  # to avoid deleting notes/ Starting...

            # do DEoptim for steptol iterations with custom printing and custom
            # convergence checks
            before <- Sys.time()
            for (i in 1:nrep) {

                out <- do.call(DEoptim::DEoptim, argsList)
                argsList$control$initialpop <- out$member$pop  # update population value
                argsList$control$itermax <- Optim$steptol[i + 1] + 1  # update population value
                argsList$control$storepopfrom <- Optim$steptol[i + 1]
                nfeval <- nfeval + out$optim$nfeval
                newval <- out$optim$bestval

                improvement <- oldval - newval

                # print progress
                if (verbose == 1) {

                  diff <- Sys.time() - before
                  currentIter <- sum(Optim$steptol[1:i])
                  totIter <- sum(Optim$steptol)
                  # based on dplyr::progress_estimated
                  avg <- diff/currentIter
                  eta <- avg * (totIter - currentIter)

                  cat(sprintf("Iteration: %d | Rel. improvement: %.3g | ETA: %.2f %s \n",
                    currentIter, (oldval - newval)/newval, eta, attr(eta,
                      "unit")))

                } else if (verbose == 2) {
                  tempOut$Bestvals <- imposeFixations(fixed = fixed, pars = out$optim$bestmem,
                    parnames = names(lower))
                  names(tempOut$Bestvals)[nchar(names(tempOut$Bestvals)) ==
                    0] <- fixed$fixedMat[1, ]

                  if (i > 1)
                    cat(sprintf("\n%s\n", paste(rep("=", nCharRep), collapse = "")))
                  msg1 <- utils::capture.output(cat(sprintf("Total iterations done: %s \nImprovement over last %s iterations: %10g \nObjective function value: %10g \nCurrent parameter estimates:\n",
                    sum(Optim$steptol[1:i]), Optim$steptol[i], improvement,
                    newval)), print(tempOut))
                  cat(paste0(msg1, collapse = "\n"))
                  if (i == 1)
                    nCharRep <- max(nchar(msg1))
                  prevSize <- sum(nchar(msg1)) + length(msg1) - 1
                }
                if (oldval - newval < Optim$reltol * newval) {
                  # check for convergence
                  break
                } else {
                  oldval <- newval
                }
            }
            niter <- sum(Optim$steptol[1:i])

        } else {
            # do DEoptim with build in convergence

            argsList$control$itermax <- Optim$itermax
            argsList$control$steptol <- Optim$steptol[1]
            argsList$control$storepopfrom <- Optim$steptol[1] + 1
            argsList$control$reltol <- Optim$reltol

            out <- do.call(DEoptim::DEoptim, argsList)
            niter <- out$optim$iter
            nfeval <- out$optim$nfeval

        }

        # gather relevant output
        Bestvals <- imposeFixations(fixed = fixed, pars = out$optim$bestmem,
            parnames = names(lower))
        if (!is.null(fixed$fixedNames)) {
            names(Bestvals)[names(Bestvals) == ""] <- fixed$fixedNames
        }
        # calculate model densities at parameter estimates
        pars <- Bestvals[c(restr.mat)]  # extract all parameters



        out <- list(Bestvals = Bestvals, fixed = fixed, GlobalOptimizer = out,
            Debug = list(niter = niter, nfeval = nfeval, itermax = Optim$itermax),
            note = note)

        if (verbose > 0) {
            cat("\nAnalyses complete!\n")
        }

    } else {
        # calculate the objective function for a given set of parameters

        argsList$pars <- pars
        argsList$all <- TRUE
        out <- list(objVals = do.call(total.objective, argsList), pars = pars,
            Bestvals = pars)

    }
    # calculate model densities at parameter estimates
    dim(pars) <- dim(restr.mat)
    pars.list <- unlist(apply(pars, 2, list), recursive = FALSE)
    m <- getPdf(pars.list = pars.list, tt = tt, DstarM = DstarM, mm = mm,
        oscPdf = FALSE, fun.density = fun.density, args.density = args.density)
    if (!is.null(m)) {
        var.m <- getVar(m, tt, mm2)
        colnames(m) <- names(rt)
    } else {
        var.m <- NULL
        warning("Solution contained improper pdfs. No decision model distributions have been saved. Perhaps rerun the analysis with a more a narrow time grid.",
            immediate. = TRUE, call. = FALSE)
    }

    conditionNames <- sort(unique(data[[condition]]))

    out2 <- list(tt = tt, g.hat = g, modelDist = m, ncondition = ncondition,
        var.data = var.data, var.m = var.m, restr.mat = restr.mat, splits = splits,
        n = n, DstarM = DstarM, fun.density = fun.density, fun.dist = fun.dist,
        h = h, args.density = args.density, args.dist = args.dist, conditionNames = conditionNames,
        formula = formula)
    out <- c(out, out2)
    class(out) <- "DstarM.fitD"
    return(out)

}

# custom functions: do the same as their names but without error handling
# ~5.5 times faster than stats::approx but devtools::check() gives
# warnings :( customApprox = function(x, fx, n.pts) { xout =
# seq.int(x[1L], x[length(x)], length.out = n.pts) yout =
# .Call(stats:::C_Approx, as.double(x), as.double(fx), xout, 1L, NA, NA,
# 0L) return(list(x = xout, y = yout)) } customApprox = function(x, fx,
# n.pts) stats::approx(x, y = fx, n = n.pts)

# R version of approx that is still ~2.75x faster than approx(x, y = fx,
# n = n.pts) designed to only work for n.pts = 2L*n - 1L (for purposes of
# speed)
customApprox <- function(x, fx, n.pts) {
    xout <- seq.int(x[1L], x[length(x)], length.out = n.pts)
    yout <- rep.int(0L, n.pts)
    yout[seq.int(1L, n.pts, 2L)] <- fx
    yout[seq.int(2L, n.pts, 2L)] <- (fx[2L:length(fx)] - fx[1L:(length(fx) -
        1L)])/2L + fx[-length(fx)]
    return(list(x = xout, y = yout))
}


customConvolveO <- function(x, y) {
    n1 <- length(y) - 1L
    y <- c(y, rep.int(0L, (n1)))
    x <- stats::fft(stats::fft(c(rep.int(0L, (n1)), x)) * Conj(stats::fft(y)),
        inverse = TRUE)
    return(Re(x)/length(y))
}
## get functions: obtain what's in the name get data based densities;
## returns matrix
getGhat <- function(rt, tt, ncondition, mm, by) {
    tt2 <- c(0, tt + by/2)
    g <- matrix(0, length(tt), 2 * ncondition)
    condRespSize <- lengths(rt)
    for (i in which(condRespSize > 0)) {
        # get frequencies for nonzero observations
        g[, i] <- graphics::hist(rt[[i]], breaks = tt2, include.lowest = TRUE,
            plot = FALSE)$counts
    }
    if (any(condRespSize == 0)) {
        # return warnings for zero observations
        warning(paste0("Zero observations in condition response pair(s) ",
            paste0(which(condRespSize == 0), collapse = ", "), ". Analysis will continue but check for errors."),
            immediate. = TRUE, call. = FALSE)
    }

    g <- g %*% (diag(dim(g)[2L])/rep(apply(g %*% mm, 2, simpson, x = tt),
        each = 2))
    return(g)
}

# get kurtosis of nondecision density | DEPRECATED?
getKurt <- function(tt, kurtf, m, var.data, var.m, var.r) {
    c0 <- kurtf * (var.data)^2  # C(a+b)^2
    c1 <- var.m^2 * kurt(tt, m)  # a^2 A
    c2 <- 6 * var.m * var.r  # 6ab
    return(c0 - c1 - c2)  #/ var.r^2) # divide by var.r^2 is redundant since there is already a restriction preventing this form becoming negative
}

# get all model densities; returns matrix | via simpson
getPdf <- function(pars.list, tt, DstarM, mm, oscPdf = TRUE, fun.density = Voss.density,
    args.density = list()) {
    # get pdfs
    pdf <- Map(fun.density, pars = rep(pars.list, 1, NA, 2), boundary = rep(c("lower",
        "upper"), length(pars.list)), MoreArgs = c(args.density, list(t = tt,
        DstarM = DstarM)))
    pdf <- do.call(cbind, pdf)
    if (oscPdf) {
        if (!all(apply(pdf, 2, oscCheck)))
            {
                # if not all pdf dists pass the oscCheck
                return(NULL)
            }  # if they all pass oscCheck
    }

    cor <- 1/apply(pdf %*% mm, 2, simpson, x = tt)  # ensure pdfs integrate to 1.
    if (any(is.infinite(cor))) {
        # return NULL if pdfs are gibberish.
        return(NULL)
    }
    pdf <- pdf %*% (diag(dim(pdf)[2L]) * rep(cor, each = 2L))
    return(pdf)
}

# get all variances; returns vector | via simpson
getVar <- function(Pdf, tt, mm2) {
    Pdf <- Pdf %*% mm2
    Pdf <- Pdf %*% (diag(dim(Pdf)[2L])/apply(Pdf, 2, simpson, x = tt))
    return(apply(Pdf, 2, nth.cmomentS, x = tt, nth = 2))
}

# helper function
groups <- function(nconditie, splits, restriction = FALSE) {
    # error handling
    stopifnot(length(splits) == nconditie)
    # create matrix where everything is split
    if (restriction) {
        ncond <- matrix(1:nconditie, ncol = length(splits))
    } else {
        ncond <- matrix(1:(nconditie * 2), nrow = 2, ncol = length(splits))
    }
    # find unique values to split by
    uniq <- unique(splits)
    ncond.temp <- vector("list", length = length(uniq))
    for (i in seq_along(uniq)) {
        # concatenate values with equal split value and put them in a list.
        ncond.temp[[i]] <- matrix(c(ncond[, which(splits == uniq[i])]))
    }
    # find maximum length
    ind.temp <- which.max(lengths(ncond.temp))
    # add NA's to fill shorter columns
    for (i in (1:length(ncond.temp))[-ind.temp]) {
        ncond.temp[[i]] <- c(ncond.temp[[i]], rep(NA, lengths(ncond.temp)[ind.temp] -
            lengths(ncond.temp)[i]))
    }
    # bind list into a matrix
    ncond.final <- do.call(cbind, ncond.temp)
    return(ncond.final)
}

# kurt == nth.cmomentS(..., nth = 4) ? DEPRECATED?
kurt <- function(x, fx, n = length(x)) nth.cmomentS(x, fx, 4, n)/nth.cmomentS(x,
    fx, 2, n)^2

# ceil to base
mceil <- function(x, base) base * ceiling(x/base)

# round to base
mround <- function(x, base) base * round(x/base)

# Estimate n-th moment of a distribution via simpson integration
nth.momentS <- function(x, fx, nth = 1, n = length(x)) simpson(x, x^nth *
    fx, n)

# Estimate n-th central moment of a distribution via simpson integration
nth.cmomentS <- function(x, fx, nth = 1, n = length(x)) simpson(x, (x - simpson(x,
    x * fx, n))^nth * fx, n)

# checks if pdf is unimodal. Only works for pdfs, returns TRUE if
# unimodal, FALSE otherwise
oscCheck <- function(pdf) {
    check <- rle(sign(zapsmall(diff(pdf))))
    return(length(check$values[check$values != 0]) == 2)
}

# simpson integral adapted from Bolstad::sintegral(); returns scalar |
# customApprox
simpson <- function(x, fx, n = length(x)) {
    ap <- customApprox(x, fx, n.pts = 2L * n - 1L)  # was + 1L
    ind <- 1L:(n - 1L)
    return(sum(ap$y[2L * ind - 1L] + 4L * ap$y[2L * ind] + ap$y[2L * ind +
        1L]) * (ap$x[2L] - ap$x[1L])/3L)
}

# works with stats::approx simpson = function(x, fx, n = length(x)) { ap
# = customApprox(x, fx, n.pts = 2*n + 1) return(simpsonC(ap$x, ap$y)) }

# objective function for D*M and Chisq analyses separate these into 2
# objective functions?
total.objective <- function(pars, tt, g, ql, restr.mat, fixed = NULL, DstarM = TRUE,
    all = FALSE, forceRestriction, oscPdf, ii, jj, mm, mm2, fun.density,
    args.density, fun.dist, args.dist, var.data, parnames = NULL, by) {
    # impose parameter fixations browser()
    pars <- imposeFixations(fixed = fixed, pars = pars, parnames = parnames)
    # get all unique parameter configurations taking restrictions into
    # account
    pars.list <- lapply(restr.mat, function(ind, pars) pars[ind], pars)

    # get pdfs and check for oscillations
    m <- getPdf(pars.list = pars.list, tt = tt, DstarM = DstarM, mm = mm,
        oscPdf = oscPdf, fun.density = fun.density, args.density = args.density)
    # getPdf returns NULL if pdfs generated fail certain sanity checks
    if (is.null(m)) {
        return(1e+09)
    }
    if (!DstarM) {
        # traditional estimation
        out <- rep.int(0L, (dim(g)[2L]))
        for (i in 1:length(out)) {
            # calc chi square dist
            args.dist$a <- g[, i]
            args.dist$b <- m[, i]
            out[i] <- do.call(fun.dist, args.dist) * 100 * ql[i]/sum(ql)
        }
    } else {
        # DstarM estimation check for violations of restriction on variance
        if (forceRestriction) {
            var.m <- getVar(m, tt, mm2)
            var.r <- var.data - var.m
            if (any(var.r < 0L | abs(var.r) < .Machine$double.neg.eps)) {
                return(1e+09 - 20 * sum(var.r[var.r < 0]))
            }
        }
        out <- rep(0, (length(ii)))
        for (l in 1:length(ii)) {
            # args.dist$a = customConvolveO(g[, ii[l]], rev(m[,
            # jj[l]]))[seq_along(tt)] args.dist$b = customConvolveO(g[, jj[l]],
            # rev(m[, ii[l]]))[seq_along(tt)]
            args.dist$a <- customConvolveO(g[, ii[l]], by * rev(m[, jj[l]]))[seq_along(tt)]
            args.dist$b <- customConvolveO(g[, jj[l]], by * rev(m[, ii[l]]))[seq_along(tt)]
            out[l] <- do.call(fun.dist, args.dist) * 100 * (ql[ii[l]] + ql[jj[l]])/sum(ql)
        }
    }
    if (all) {
        # return unsummed distance vector for debugging purposes
        return(out)
    } else {
        # return summed vector for differential evolution algorithm
        return(sum(out))
    }
}

# imposes fixations from fixed onto pars
imposeFixations <- function(fixed, pars, parnames) {

    if (fixed$anyFixed) {
        for (i in 1:length(fixed$indFixed)) {
            # loop over everything to be looked up

            if (fixed$isNumeric[i]) {
                # numeric replacement

                insert <- as.numeric(fixed$fixedMat[3L, i])

            } else {
                # expression (i.e. z = a/2)

                replacement <- pars[which(parnames == fixed$fixedMat[3L,
                  i])][1L]
                insert <- eval(parse(text = gsub(pattern = fixed$fixedMat[3L,
                  i], replacement = replacement, x = fixed$fixedMat[2L, i])))

            }

            # execute operation to which parameter was fixed

            pars <- append(pars, insert, fixed$indFixed[i] - 1L)  # append fixations to parameter vector
            # names(pars)[fixed$indFixed[i]] = fixed$fixedMat[1L, i] # add names to
            # parameter vector
        }
    }

    return(pars)

}
