#' Bootstrap a conditional multivariate extreme values model
#'
#' Bootstrap a conditional multivariate extreme values model following the
#' method of Heffernan and Tawn, 2004.
#'
#' Details of the bootstrap method are given by Heffernan and Tawn (2004). The
#' procedure is semi-parametric.
#'
#' Firstly, values of all variables are simulated independently from the
#' parametric Gumbel or Laplace distributions (depending on the choice of
#' \code{margins} in the original call to \code{\link{mex}}). The sample size
#' and data dimension match that of the original data set.  Then an empirical
#' bootstrap sample is generated from the original data after its
#' transformation to the Gumbel/Laplace scale. Again, sample size and structure
#' match the original data set. The empirical bootstrap samples from each
#' margin are then sorted, and then replaced by their corresponding values from
#' the sorted Gumbel/Laplace samples. This procedure preserves the dependence
#' structure of the empirical bootstrap sample while ensuring the marginal
#' properties of the resulting semi-parametric bootstrap sample are those of
#' the parametric Gumbel/Laplace distribution.
#'
#' The simulated, ordered Laplace/Gumbel sample is then transformed to the
#' scale of the original data by using the Probability Integral Transform.
#' Values beneath the original thresholds for fitting of the GPD tail models
#' are transformed by using the empirical distribution functions and for values
#' above these thresholds, the fitted GPDs are used.  This completes the
#' semi-parametric bootstrap from the data.
#'
#' Parameter estimation is then carried out as follows: The parameters in the
#' generalized Pareto distributions are estimated by using the bootrap data,
#' these data are then transformed to the Laplace/Gumbel scale using the
#' orginal threshold, their empirical distribution function and these estimated
#' GPD parameters. The variables in the dependence structure of these variables
#' are then estimated.
#'
#' Note that maximum likelihood estimation will often fail for small samples
#' when the generalized Pareto distribution is being fit. Therefore it will
#' often be useful to use penalized likelihood estimation. The function
#' \code{bootmex} does whatever was done in the call to \code{migpd} or
#' \code{mex} that generated the object with which it is being called.
#'
#' Also note that sometimes (again, usually with small data sets) all of the
#' simulated Laplace/Gumbel random numbers will be beneath the threshold for
#' the conditioning variable. Such samples are abandoned by \code{bootmex} and
#' a new sample is generated. This probably introduces some bias into the
#' resulting bootstrap distributions.
#'
#' The \code{plot} method produces histograms of bootstrap gpd parameters (the
#' default) or scatterplots of dependence parameters with the point estimates
#' for the original data shown.
#'
#' By design, there is no \code{coef} method. The bootstrapping is done to
#' account for uncertainty. It is not obvious that adjusting the parameters for
#' the mean bias is the correct thing to do.
#'
#' @aliases bootmex print.bootmex plot.bootmex
#' @usage bootmex(x, R = 100, nPass=3, trace=10,referenceMargin=NULL)
#'
#' \method{plot}{bootmex}(x, plots = "gpd", main = "", ...)
#' \method{print}{bootmex}(x, ...)
#' @param x An object of class "mex" as returned by function \code{\link{mex}}.
#' @param R The number of bootstrap runs to perform. Defaults to \code{R}=100.
#' @param nPass An integer. Sometimes, particularly with small samples, the
#' estimation process fails with some bootstrap samples. The function checks
#' which runs fail and takes additional bootstrap samples in an attempt to get
#' parameter estimates. By default, it has nPass=3 attempts at this before
#' giving up.
#' @param trace How often to inform the user of progress. Defaults to
#' \code{trace=10}.
#' @param referenceMargin Optional set of reference marginal distributions to use
#'   for marginal transformation if the data's own marginal distribution is not
#'   appropriate (for instance if only data for which one variable is large is
#'   available, the marginal distributions of the other variables will not be
#'   represented by the available data).  This object can be created from a
#'   combination of datasets and fitted GPDs using the function
#'   \code{makeReferenceMarginalDistribution}.
#' @param plots What type of diagnostic plots to produce.  Defaults to "gpd" in
#' which case gpd parameter estimate plots are produced otherwise plots are
#' made for the dependence parameters.
#' @param main Title for plots.
#' @param ... Further arguments to be passed to methods.
#' @return An object of class 'bootmex'. Print and plot functions are
#' available.
#' @author Harry Southworth
#' @seealso \code{\link{migpd}}, \code{\link{mexDependence} },
#' \code{\link{bootmex}}, \code{\link{predict.mex}}.
#' @references J. E. Heffernan and J. A. Tawn, A conditional approach for
#' multivariate extreme values, Journal of the Royal Statistical society B, 66,
#' 497 -- 546, 2004
#' @keywords models multivariate
#' @examples
#'
#' \donttest{
#' mymex <- mex(winter , mqu = .7, dqu = .7, which = "NO")
#' myboot <- bootmex(mymex)
#' myboot
#' plot(myboot,plots="gpd")
#' plot(myboot,plots="dependence")
#' }
#' @export bootmex
bootmex <-
    # Bootstrap inference for a conditional multivaratiate extremes model.
function (x, R = 100, nPass = 3, trace = 10,referenceMargin=NULL) {
    theCall <- match.call()
    if (!inherits(x, "mex")){
      stop("object must be of type 'mex'")
    }

# Construct the object to be returned.
    ans <- list()
    ans$call <- theCall

    getTran <- function(i, x, data, mod, th, qu, margins) {
        param <- mod[[i]]$coefficients
        revTransform(margins$q2p(c(x[, i])), data = c(data[, i]), th = th[i],
                     qu = qu[i], sigma = exp(param[1]), xi = param[2])
    }

    mar <- x$margins
    dep <- x$dependence
    which <- dep$which
    constrain <- dep$constrain
    v <- dep$v
    dqu <- dep$dqu
    dth <- dep$dth
    margins <- dep$margins
    penalty <- mar$penalty
    priorParameters <- mar$priorParameters
    start <- 0.75* coef(x)$dependence[1:2,] # scale back towards zero in case point est on edge of original parameter space and falls off edge of constrained space for bootstrap sample

    n <- dim(mar$transformed)[[1]]
    d <- dim(mar$transformed)[[2]]
    dqu <- rep(dqu, d)
    dependent <- (1:d)[-which]

    ans$simpleDep <- dep$coefficients
    ans$dqu <- dqu
    ans$which <- which
    ans$R <- R
    ans$simpleMar <- mar
    ans$margins <- margins
    ans$constrain <- constrain

    innerFun <- function(i, x, which, dth, dqu, margins, penalty, priorParameters, constrain, v=v, start=start,
        pass = 1, trace = trace, n=n, d=d, getTran=getTran, dependent=dependent,referenceMargin=referenceMargin) {

        g <- sample(1:(dim(mar$transformed)[[1]]), size = n, replace = TRUE)
        g <- mar$transformed[g, ]
        ok <- FALSE

        while (!ok) {
          for (j in 1:(dim(g)[[2]])){
            u <- runif(nrow(g))
            g[order(g[, j]), j] <- sort(margins$p2q(u))
          }
          if (sum(g[, which] > dth) > 1  &   all(g[g[,which] > dth , which] > 0)){ ok <- TRUE }
        }

        g <- sapply(1:d, getTran, x = g, data = mar$data, margins=margins,
                    mod = mar$models, th = mar$mth, qu = mar$mqu)

        dimnames(g)[[2]] <- names(mar$models)

        ggpd <- migpd(g, mth = mar$mth,
                      penalty = penalty, priorParameters = priorParameters)

        gd <- mexDependence(ggpd, dqu = dqu, which = which, margins=margins[[1]], constrain=constrain, v=v, start=start,referenceMargin=referenceMargin)
        res <- list(GPD = coef(ggpd)[3:4, ],
                    dependence = gd$dependence$coefficients,
                    Z = gd$dependence$Z,
                    Y = g)

        if (pass == 1) {
            if (i%%trace == 0) {
                message(paste(i, "replicates done\n"))
            }
        }
        res
    } # Close innerFun

    res <- lapply(1:R, innerFun, x = x, which = which, dth = dth, margins=margins,
        dqu = dqu, penalty = penalty, priorParameters = priorParameters, constrain=constrain, v=v, start=start,
        pass = 1, trace = trace, getTran=getTran, n=n, d=d, dependent=dependent,referenceMargin=referenceMargin)

    # Sometimes samples contain no extreme values. Need to have another pass or two
    if (nPass > 1) {
        for (pass in 2:nPass) {
            rerun <- sapply(res, function(x) any(sapply(x, function(x) any(is.na(x)))))
            wh <- !unlist(lapply(res, function(x) dim(x$Z)[[1]] > 0))
            rerun <- apply(cbind(rerun, wh), 1, any)
            if (sum(rerun) > 0) {
                message("Pass ", pass, " : ", sum(rerun), " samples to rerun.\n")
                rerun <- (1:R)[rerun]
                res[rerun] <- lapply((1:R)[rerun], innerFun,
                  x = x, which = which, dth = dth, dqu = dqu, margins=margins,
                  penalty = penalty, priorParameters = priorParameters, constrain=constrain, v=v, start=start,
                  pass = pass, trace = trace, getTran=getTran, n=n, d=d, dependent=dependent,referenceMargin=referenceMargin)
            }
        }
    }

    ans$boot <- res
    oldClass(ans) <- c("bootmex", "mex")
    ans
}

