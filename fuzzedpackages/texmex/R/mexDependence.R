#' Estimate the dependence parameters in a conditional multivariate extreme
#' values model
#'
#' Estimate the dependence parameters in a conditional multivariate extreme
#' values model using the approach of Heffernan and Tawn, 2004.
#'
#' Estimates the extremal dependence structure of the data in \code{x}.  The
#' precise nature of the estimation depends on the value of \code{margins}.  If
#' \code{margins="laplace"} (the default) then dependence parameters a and b
#' are estimated after transformation of the data to Laplace marginal
#' distributions.  These parameters can describe both positive and negative
#' dependence.  If \code{margins="gumbel"} then the parameters a, b, c and d in
#' the dependence structure described by Heffernan and Tawn (2004) are
#' estimated in the following two steps: first, a and b are estimated; then, if
#' a=0 and b is negative, parameters c and d are estimated (this is the case of
#' negative dependence). Otherwise c and d will be fixed at zero (this is the
#' case of positive dependence).
#'
#' If \code{margins="laplace"} then the option of constrained parameter
#' estimation is available by setting argument \code{constrain=TRUE}.  The
#' default is to constrain the values of the parameters
#' (\code{constrain=TRUE}).  This constrained estimation ensures validity of
#' the estimated model, and enforces the consistency of the fitted dependence
#' model with the strength of extremal dependence exhibited by the data.  More
#' details are given in Keef et al. (2013).  The effect of this constraint is
#' to limit the shape of the dependence parameter space so that its boundary is
#' curved rather than following the original box constraints suggested by
#' Heffernan and Tawn (2004).  The constraint brings with it some performance
#' issues for the optimiser used to estimate the dependence parameters, in
#' particular sensitivity to choice of starting value which we describe now.
#'
#' The dependence parameter estimates returned by this function can be
#' particularly sensitive to the choice of starting value used for the
#' optimisation.  This is especially true when \code{margins="laplace"} and
#' \code{constrain=TRUE}, in which case the maximum of the objective function
#' can lie on the edge of the (possibly curved) constrained parameter space.
#' It is therefore up to the user to check that the reported parameter
#' estimates really do correspond to the maximum of the profile lilkelihood
#' surface.  This is easily carried out by using the visual diagnostics invoked
#' by setting \code{PlotLikDo=TRUE} and adjusting the plotting area by using
#' the argument \code{PlotLikRange} to focus on the region containing the
#' surface maximum.  See an example below which illustrates the use of this
#' diagnostic.
#'
#' @usage mexDependence(x, which, dqu, margins="laplace",
#'     constrain=TRUE, v = 10, maxit=1000000, start=c(.01, .01),
#'     marTransform="mixture", referenceMargin = NULL, nOptim = 1,
#'     PlotLikDo=FALSE, PlotLikRange=list(a=c(-1,1),b=c(-3,1)),
#'     PlotLikTitle=NULL)
#' @param x An object of class "migpd" as returned by
#'     \code{\link{migpd}}.
#' @param which The name of the variable on which to condition. This
#'     is the name of a column of the data that was passed into
#'     \code{migpd}.
#' @param dqu See documentation for this argument in
#'     \code{\link{mex}}.
#' @param margins The form of margins to which the data are
#'     transformed for carrying out dependence estimation.  Defaults
#'     to "laplace", with the alternative option being "gumbel".  The
#'     choice of margins has an impact on the interpretation of the
#'     fitted dependence parameters.  Under Gumbel margins, the
#'     estimated parameters a and b describe only positive dependence,
#'     while c and d describe negative dependence in this case.  For
#'     Laplace margins, only parameters a and b are estimated as these
#'     capture both positive and negative dependence.
#' @param constrain Logical value.  Defaults to \code{constrain=TRUE}
#'     although this will subsequently be changed to FALSE if
#'     \code{margins="gumbel"} for which constrained estimation is not
#'     implemented.  If \code{margins="laplace"} and
#'     \code{constrain=TRUE} then the dependence parameter space is
#'     constrained to allow only combinations of parameters which give
#'     the correct stochastic ordering between (positively and
#'     negatively) asymptotically dependent variables and variables
#'     which are asymptotically independent.
#' @param v Scalar. Tuning parameter used to carry out constrained
#'     estimation of dependence structure under
#'     \code{constrain=TRUE}. Takes positive values greater than 1;
#'     values between 2 and 10 are recommended.
#' @param maxit The maximum number of iterations to be used by the
#'     optimizer.  Defaults to \code{maxit = 1000000}.
#' @param start Optional starting value for dependence estimation.
#'     This can be: a vector of length two, with values corresponding
#'     to dependence parameters a and b respectively, and in which
#'     case \code{start} is used as a starting value for numerical
#'     estimation of each of the dependence models to be estimated; a
#'     matrix with two rows corresponding to dependence parameters a
#'     and b respectively and number of columns equal to the number of
#'     dependence models to be estimated (the ordering of the columns
#'     will be as in the original data matrix); or a previously
#'     estimated object of class "mex" whose dependence parameter
#'     estimates are used as a starting point for estimation.  Note
#'     that under \code{constrain=TRUE}, if supplied, \code{start}
#'     must lie within the permitted area of the parameter space.
#' @param marTransform Optional form of transformation to be used for
#'     probability integral transform of data from original to Gumbel
#'     or Laplace margins.  Takes values \code{marTransform="mixture"}
#'     (the default) or \code{marTransform="empirical"}. When
#'     \code{marTransform="mixture"}, the rank transform is used below
#'     the corresponding GPD fitting threshold used in \code{x}, and
#'     the fitted gpd tail model is used above this threshold.  When
#'     \code{marTransform="empirical"} the rank transform is used for
#'     the entire range of each marginal distribution.
#' @param referenceMargin Optional set of reference marginal
#'     distributions to use for marginal transformation if the data's
#'     own marginal distribution is not appropriate (for instance if
#'     only data for which one variable is large is available, the
#'     marginal distributions of the other variables will not be
#'     represented by the available data).  This object can be created
#'     from a combination of datasets and fitted GPDs using the
#'     function \code{makeReferenceMarginalDistribution}.
#' @param nOptim Number of times to run optimiser when estimating
#'     dependence model parameters. Defaults to 1.  In the case of
#'     \code{nOptim > 1} the first call to the optimiser uses the
#'     value \code{start} as a starting point, while subsequent calls
#'     to the optimiser are started at the parameter value to which
#'     the previous call converged.
#' @param PlotLikDo Logical value: whether or not to plot the profile
#'     likelihood surface for dependence model parameters under
#'     constrained estimation.
#' @param PlotLikRange This is used to specify a region of the
#'     parameter space over which to plot the profile log-likelihood
#'     surface.  List of length 2; each item being a vector of length
#'     two corresponding to the plotting ranges for dependence
#'     parameters a and b respectively. If this argument is not
#'     missing, then \code{PlotLikDo} is set equal to TRUE.
#' @param PlotLikTitle Used only if \code{PlotLikDo=TRUE}.  Character
#'     string.  Optional title added to the profile log-likelihood
#'     surface plot.
#' @return An object of class \code{mex} which is a list containing
#'     the following three objects: \item{margins}{An object of class
#'     \code{\link{migpd}}.} \item{dependence}{An object of class
#'     \code{\link{mexDependence}}.} \item{call}{This matches the
#'     original function call.}
#' @author Harry Southworth, Janet E. Heffernan
#' @seealso \code{\link{migpd}}, \code{\link{bootmex}},
#'     \code{\link{predict.mex}}, \code{\link{plot.mex}}
#' @references J. E. Heffernan and J. A. Tawn, A conditional approach
#'     for multivariate extreme values, Journal of the Royal
#'     Statistical society B, 66, 497 -- 546, 2004.
#'
#' C. Keef, I. Papastathopoulos and J. A. Tawn.  Estimation of the conditional
#' distribution of a multivariate variable given that one of its components is
#' large: Additional constraints for the Heffernan and Tawn model, Journal of
#' Multivariate Analysis, 115, 396 -- 404, 2013
#' @keywords models multivariate
#' @examples
#'
#' data(winter)
#' mygpd <- migpd(winter , mqu=.7, penalty="none")
#' mexDependence(mygpd , which = "NO", dqu=.7)
#'
#' # focus on 2-d example with parameter estimates on boundary of constrained parameter space:
#' NO.NO2 <- migpd(winter[,2:3] , mqu=.7, penalty="none")
#'
#' # starting value gives estimate far from true max:
#' mexDependence(NO.NO2, which = "NO",dqu=0.7,start=c(0.01,0.01),
#'               PlotLikDo=TRUE,PlotLikTitle=c("NO2 | NO"))
#'
#' # zoom in on plotting region containing maximum:
#' mexDependence(NO.NO2, which = "NO",dqu=0.7,start=c(0.01,0.01),
#'               PlotLikDo=TRUE,PlotLikTitle=c("NO2 | NO"),
#'               PlotLikRange = list(a=c(0,0.8),b=c(-0.2,0.6)))
#'
#' # try different starting value:
#' mexDependence(NO.NO2, which = "NO",dqu=0.7,start=c(0.1,0.1),
#'               PlotLikDo=TRUE,PlotLikTitle=c("NO2 | NO"),
#'               PlotLikRange = list(a=c(0,0.8),b=c(-0.2,0.6)))
#'
#'
#' @export mexDependence
`mexDependence` <-
    function (x, which, dqu, margins = "laplace",
              constrain=TRUE, v = 10, maxit=1000000,
              start=c(.01, .01), marTransform="mixture",
              referenceMargin=NULL, nOptim = 1,
              PlotLikDo=FALSE, PlotLikRange=list(a=c(-1,1),b=c(-3,1)),
              PlotLikTitle=NULL){
   theCall <- match.call()
   if (!inherits(x, "migpd"))
       stop("you need to use an object created by migpd")

   margins <- list(casefold(margins),
                   p2q = switch(casefold(margins),
                                "gumbel" = function(p)-log(-log(p)),
                                "laplace" = function(p)ifelse(p < .5, log(2 * p), -log(2 * (1 - p)))),
                   q2p = switch(casefold(margins),
                                "gumbel" = function(q)exp(-exp(-q)),
                                "laplace" = function(q)ifelse(q < 0, exp(q)/2, 1- 0.5*exp(-q))))

   x <- mexTransform(x, margins = margins, method = marTransform, r=referenceMargin)
   x$referenceMargin <- referenceMargin

   if (margins[[1]] == "gumbel" & constrain){
     warning("With Gumbel margins, you can't constrain, setting constrain=FALSE")
     constrain <- FALSE
   }

   if (missing(which)) {
       message("Missing 'which'. Conditioning on", dimnames(x$transformed)[[2]][1], "\n")
       which <- 1
   }
   else if (length(which) > 1)
       stop("which must be of length 1")
   else if (is.character(which))
       which <- match(which, dimnames(x$transformed)[[2]])

   if (missing(dqu)) {
       message("Assuming same quantile for dependence thesholding as was used\n     to fit corresponding marginal model...\n")
       dqu <- x$mqu[which]
   }
   dth <- quantile(x$transformed[, which], dqu)

   dependent <- (1:(dim(x$data)[[2]]))[-which]
   if (length(dqu) < length(dependent))
       dqu <- rep(dqu, length = length(dependent))

   # Allowable range of 'a' depends on marginal distributions
   aLow <- ifelse(margins[[1]] == "gumbel", 10^(-10),-1 + 10^(-10))

   if (missing(start)){
     start <- c(.01, .01)
   } else if(inherits(start, "mex")){
     start <- start$dependence$coefficients[1:2,]
   }

   if( length(start) == 2 ){
     start <- matrix(rep(start,length(dependent)),nrow=2)
   }

   if( length(start) != 2*length(dependent)){
     stop("start should be of type 'mex' or be a vector of length 2, or be a matrix with 2 rows and ncol equal to the number of dependence models to be estimated")
   }

   if( ! missing(PlotLikRange) ){
     PlotLikDo <- TRUE
   }

   qfun <- function(X, yex, wh, aLow, margins, constrain, v, maxit, start){
     Qpos <- function(param, yex, ydep, constrain, v, aLow) {

  	   a <- param[1]
       b <- param[2]

       res <- PosGumb.Laplace.negProfileLogLik(yex, ydep, a, b, constrain, v, aLow) # defined in file mexDependenceLowLevelFunctions
       res$profLik
     } # Close Qpos <- function

     o <- try(optim(par=start, fn=Qpos,
              control=list(maxit=maxit),
              yex = yex[wh], ydep = X[wh], constrain=constrain, v=v, aLow=aLow), silent=TRUE)

     if (inherits(o, "try-error")){
        warning("Error in optim call from mexDependence")
        o <- as.list(o)
        o$par <- rep(NA, 6)
        o$value <- NA
     } else if (o$convergence != 0) {
        warning("Non-convergence in mexDependence")
        o <- as.list(o)
        o$par <- rep(NA, 6)

     } else if(nOptim > 1) {

        for( i in 2:nOptim ){
           o <- try(optim(par=o$par, fn=Qpos,
                    control=list(maxit=maxit),
                    yex = yex[wh], ydep = X[wh], constrain=constrain, v=v, aLow=aLow), silent=TRUE)
           if (inherits(o, "try-error")){
             warning("Error in optim call from mexDependence")
             o <- as.list(o)
             o$par <- rep(NA, 6)
             o$value <- NA
             break()
           } else if (o$convergence != 0) {
             warning("Non-convergence in mexDependence")
             o <- as.list(o)
             o$par <- rep(NA, 6)
             break()
           }
        }
     }

     if ( PlotLikDo ){# plot profile likelihood for (a,b)
       nGridPlotLik <- 50
       a.grid <- seq(PlotLikRange$a[1],PlotLikRange$a[2],length=nGridPlotLik)
       b.grid <- seq(PlotLikRange$b[1],PlotLikRange$b[2],length=nGridPlotLik)
       NegProfLik <- matrix(0,nrow=nGridPlotLik,ncol=nGridPlotLik)
       for(i in 1:nGridPlotLik){
         for(j in 1:nGridPlotLik){
           NegProfLik[i,j] <- PosGumb.Laplace.negProfileLogLik(yex=yex[wh], ydep=X[wh],
                                  a = a.grid[i],b=b.grid[j], constrain=constrain,v=v,aLow=aLow)$profLik
         }
       }
       NegProfLik[NegProfLik > 10^10] <- NA
       if(sum(!is.na(NegProfLik))){
          filled.contour(a.grid,b.grid,-NegProfLik,main=paste("Profile likelihood",PlotLikTitle),color.palette = terrain.colors,
                         xlab="a",ylab="b",plot.axes={ axis(1); axis(2); points(o$par[1],o$par[2]) })
       }
     }

     if (!is.na(o$par[1])) { # gumbel margins and negative dependence
        if (margins == "gumbel" & o$par[1] <= 10^(-5) & o$par[2] < 0) {
           lo <- c(10^(-10), -Inf, -Inf, 10^(-10), -Inf, 10^(-10))
           Qneg <- function(yex, ydep, param) {
               param <- param[-1]
               b <- param[1]
               cee <- param[2]
               d <- param[3]
               m <- param[4]
               s <- param[5]

               obj <- function(yex, ydep, b, cee, d, m, s) {
                      mu <- cee - d * log(yex) + m * yex^b
                      sig <- s * yex^b
                      log(sig) + 0.5 * ((ydep - mu)/sig)^2
                      }
               res <- sum(obj(yex, ydep, b, cee, d, m, s))
               res
          }
          o <- try(optim(c(0, 0, 0, 0, 0, 1), Qneg, method = "L-BFGS-B", lower=lo,
                   upper=c(1, 1-10^(-10), Inf, 1-10^(-10), Inf, Inf),
                   yex = yex[wh], ydep = X[wh]), silent=TRUE)

          if (inherits(o, "try-error") || o$convergence != 0) {
             warning("Non-convergence in mexDependence")
             o <- as.list(o)
             o$par <- rep(NA, 6)
          }
        } else { # end if gumbel margins and neg dependence
          Z <- (X[wh] - yex[wh] * o$par[1]) / (yex[wh]^o$par[2])
          o$par <- c(o$par[1:2], 0, 0, mean(Z),sd(Z))
        }
    }
    c(o$par[1:6], o$value) # Parameters and negative loglik
   } # Close qfun <- function(

   yex <- c(x$transformed[, which])
   wh <- yex > unique(dth)

   res <- sapply(1:length(dependent),
                 function(X,dat,yex,wh,aLow,margins,constrain,v,maxit,start)qfun(dat[,X],yex,wh,aLow,margins,constrain,v,maxit,start[,X]),
                 dat=as.matrix(x$transformed[, dependent]), yex=yex, wh=wh, aLow=aLow, margins=margins[[1]],
                 constrain=constrain, v=v, maxit=maxit, start=start)

   loglik <- -res[7,]
   res <- matrix(res[1:6,], nrow=6)

   dimnames(res)[[1]] <- c(letters[1:4],"m","s")
   dimnames(res)[[2]] <- dimnames(x$transformed)[[2]][dependent]
   gdata <- as.matrix(x$transformed[wh, -which])
   tfun <- function(i, data, yex, a, b, cee, d) {
       data <- data[, i]
       a <- a[i]
       b <- b[i]
       cee <- cee[i]
       d <- d[i]
       if (is.na(a))
           rep(NA, length(data))
       else {
           if (a < 10^(-5) & b < 0)
               a <- cee - d * log(yex)
           else a <- a * yex
           (data - a)/(yex^b)
       }
   }
   z <- try(sapply(1:(dim(gdata)[[2]]), tfun, data = gdata,
       yex = yex[wh], a = res[1, ], b = res[2, ], cee = res[3, ], d = res[4, ]))
   if (inherits(z, c("Error", "try-error"))) {
       z <- matrix(nrow = 0, ncol = dim(x$data)[[2]] - 1)
   }
   else if (!is.array(z)) {
       z <- matrix(nrow = 0, ncol = dim(x$data)[[2]] - 1)
   }
   dimnames(z) <- list(NULL,dimnames(x$transformed)[[2]][dependent])
   res2 <- list(coefficients = res, Z = z, dth = unique(dth),
               dqu = unique(dqu), which = which, conditioningVariable= colnames(x$data)[which],
	             loglik=loglik, margins=margins, constrain=constrain, v=v)
   oldClass(res2) <- "mexDependence"

   output <- list(margins=x, dependence=res2, call=theCall)
   oldClass(output) <- "mex"
   output
}

