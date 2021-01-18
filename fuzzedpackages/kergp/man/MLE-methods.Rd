\name{mle-methods}
\docType{methods}
\alias{mle,covAll-method}
\alias{mle-methods}

\title{
  Maximum Likelihood Estimation of Gaussian Process Model Parameters
}

\description{

  Maximum Likelihood estimation of Gaussian Process models which
  covariance structure is described in a covariance kernel object.
  
}
\usage{

\S4method{mle}{covAll}(object, 
    y, X, F = NULL, beta = NULL,
    parCovIni = coef(object),
    parCovLower = coefLower(object), 
    parCovUpper = coefUpper(object),
    noise = TRUE, varNoiseIni = var(y) / 10,
    varNoiseLower = 0, varNoiseUpper = Inf,
    compGrad = hasGrad(object),
    doOptim = TRUE,
    optimFun = c("nloptr::nloptr", "stats::optim"),
    optimMethod = ifelse(compGrad, "NLOPT_LD_LBFGS", "NLOPT_LN_COBYLA"),
    optimCode = NULL,
    multistart = 1,
    parTrack = FALSE, trace  = 0, checkNames = TRUE,
    ...)
}

\arguments{
  \item{object}{

    An object representing a covariance kernel.

  }
  \item{y}{

    Response vector.

  }
  \item{X}{

    Spatial (or input) design matrix.

  } 
  \item{F}{

    Trend matrix.

  }
  \item{beta}{

    Vector of trend coefficients if known/fixed.

  }
  \item{parCovIni}{

    Vector with named elements or matrix giving the initial values for the
    parameters. See \bold{Examples}. When this argument is omitted, the
    vector of covariance parameters given in \code{object} is
    used if \code{multistart == 1}; If \code{multistart > 1},
    a matrix of parameters is simulated by using \code{\link{simulPar}}. 
    Remind that you can use the \code{coef} and \code{coef<-}
    methods to get and set this slot of the covariance object.

  }
  \item{parCovLower}{

    Lower bounds for the parameters. When this argument is omitted, the
    vector of parameters lower bounds in the covariance given in
    \code{object} is used. You can use \code{coefLower} and
    \code{coefLower<-} methods to get and set this slot of the
    covariance object.

  }  
  \item{parCovUpper}{

    Upper bounds for the parameters. When this argument is omitted, the
    vector of parameters lower bounds in the covariance given in
    \code{object} is used. You can use \code{coefUpper} and
    \code{coefUpper<-} methods to get and set this slot of the
    covariance object.

  }
  \item{noise}{

    Logical. Should a noise be added to the error term?

  }  
  \item{varNoiseIni}{

    Initial value for the noise variance.

  }  
  \item{varNoiseLower}{

    Lower bound for the noise variance. Should be \code{<= varNoiseIni}.

  }  
  \item{varNoiseUpper}{

    Upper bound for the noise variance. Should be \code{>= varNoiseIni}.

  }  
  %% \item{parFixed}{

  %%   Vector with named elements giving the fixed parameters (if any).
  %%   Not implemented yet.

  %% }
  \item{compGrad}{

    Logical: compute and use the analytical gradient in optimization?
    This is only possible when \code{object} provides the analytical
    gradient.
    %% Must be TRUE at this stage.

  }
  \item{doOptim}{

    Logical. If \code{FALSE} no optimization is done.

  }
  \item{optimFun}{

    Function used for optimization. The two pre-defined choices are
    \code{nloptr::nloptr} (default) and \code{stats::optim}, both in
    combination with a few specific optimization methods. 
    Ignored if \code{optimCode} is provided.

  }
  \item{optimMethod}{

    Name of the optimization method or algorithm. This is passed as the
    \code{"algorithm"} element of the \code{opts} argument when
    \code{nloptr::nloptr} is used (default), or to the \code{method}
    argument when \code{stats::optim} is used. When another value of
    \code{optimFun} is given, the value of \code{optimMethod} is
    ignored. Ignored if \code{optimCode} is provided. Use
    \code{\link{optimMethods}} to obtain the list of usable values.

  }
  \item{optimCode}{

    An object with class \code{"expression"} or \code{"character"}
    representing a user-written R code to be parsed and performing the
    log-likelihood maximization. Notice that this argument will bypass 
    \code{optimFun} and \code{optimMethod}. The expression must define
    an object named \code{"opt"}, which is either a list containing 
    optimization results, either an object inheriting from \code{"try-error"} 
    to cope with the case where a problem occurred during the optimization. 
    
  }

  \item{multistart}{

    Number of optimizations to perform from different starting points
    (see \code{parCovIni}).  Parallel backend is encouraged.

  }
  %% \item{method}{

  %%   To be passed to \code{\link{optim}}.

  %% }
  %% \item{control}{

  %%   To be passed to \code{\link{optim}}. The element \code{fnscale = -1}
  %%   states that a maximization is required and must not be changed.

  %% }
  \item{parTrack}{

    If \code{TRUE}, the parameter vectors used during the optimization
    are tracked and returned as a matrix.

  }
  \item{trace}{

    Integer level of verbosity.
    
  }
  \item{checkNames}{

    if \code{TRUE} (default), check the compatibility of \code{X} with
    \code{object}, see \code{\link{checkX}}.

  }
  \item{\dots}{

    Further arguments to be passed to the optimization function,
    \code{\link[nloptr]{nloptr}} or
    \code{\link[stats]{optim}}.
    
  }
  
}

\details{

  The choice of optimization method is as follows.
  
  \itemize{

    \item When \code{optimFun} is \code{nloptr:nloptr}, it is assumed
    that we are minimizing the negative log-likelihood \eqn{- \log
    L}{-log L}. Note that both predefined methods \code{"NLOPT_LD_LBFGS"} 
    and \code{"NLOPT_LN_COBYLA"} can cope with a non-finite value
    of the objective, except for the initial value of the parameter.
    Non-finite values of \eqn{- \log L}{-log L} are often met in
    practice during optimization steps. The method
    \code{"NLOPT_LD_LBFGS"} used when \code{compGrad} is \code{TRUE}
    requires that the gradient is provided by/with the covariance
    object. You can try other values of \code{optimMethod} corresponding
    to the possible choice of the \code{"algorithm"} element in the
    \code{opts} argument of \code{nloptr:nloptr}. It may be useful
    to give other options in order to control the optimization and its
    stopping rule.

    \item When \code{optimFun} is \code{"stats:optim"}, it is assumed
    that we are maximizing the log-likelihood \eqn{\log L}{log L}. 
    We suggest to use one of the methods \code{"L-BFGS-B"} or \code{"BFGS"}. 
    Notice that \code{control} can be provided in \code{\dots}, but 
    \code{control$fnscale} is forced to be \code{- 1}, corresponding to maximization. 
    Note that \code{"L-BFGS-B"} uses box constraints, but the optimization stops
    as soon as the log-likelihood is non-finite or \code{NA}. The method
    \code{"BFGS"} does not use constraints but allows the log-likelihood
    to be non-finite or \code{NA}. Both methods can be used without
    gradient or with gradient if \code{object} allows this.

  }
  
  The vectors \code{parCovIni}, \code{parCovLower}, \code{parCovUpper}
  must have elements corresponding to those of the vector of kernel
  parameters given by \code{coef(object)}. These vectors should have
  suitably named elements.
  
}
\value{

  A list with elements hopefully having understandable names. 

  \item{opt}{

    List of optimization results if it was successful, or an error
    object otherwise.

  }
  \item{coef.kernel}{

    The vector of 'kernel' coefficients. This will include one or
    several variance parameters.

  }
  \item{coef.trend}{

    Estimate of the vector \eqn{\boldsymbol{\beta}}{\beta} of the trend
    coefficients.

  }
  \item{parTracked}{

    A matrix with rows giving the successive iterates during
    optimization if the \code{parTrack} argument was set to \code{TRUE}.

  }

}

\author{

  Y. Deville, O. Roustant

}
\note{
  
  The checks concerning the parameter names, dimensions of provided
  objects, \dots are not fully implemented yet.

  Using the \code{optimCode} possibility requires a bit of programming
  effort, although a typical code only contains a few lines.

}

\seealso{
  
  \code{\link{gp}} for various examples, \code{\link{optimMethods}}
  to see the possible values of the argument \code{optimMethod}.

}

\examples{

set.seed(29770)

##=======================================================================
## Example. A 4-dimensional "covMan" kernel
##=======================================================================
d <- 4
myCovMan <- 
      covMan(
         kernel = function(x1, x2, par) { 
         htilde <- (x1 - x2) / par[1]
         SS2 <- sum(htilde^2)
         d2 <- exp(-SS2)
         kern <- par[2] * d2
         d1 <- 2 * kern * SS2 / par[1]            
         attr(kern, "gradient") <- c(theta = d1,  sigma2 = d2)
         return(kern)
      },
      label = "myGauss",
      hasGrad = TRUE,
      d = 4,    
      parLower = c(theta = 0.0, sigma2 = 0.0),
      parUpper = c(theta = +Inf, sigma2 = +Inf),
      parNames = c("theta", "sigma2"),
      par = c(NA, NA)
      )
kernCode <- "
       SEXP kern, dkern;
       int nprotect = 0, d;
       double SS2 = 0.0, d2, z, *rkern, *rdkern;

       d = LENGTH(x1);
       PROTECT(kern = allocVector(REALSXP, 1)); nprotect++;
       PROTECT(dkern = allocVector(REALSXP, 2)); nprotect++;
       rkern = REAL(kern);
       rdkern = REAL(dkern);

       for (int i = 0; i < d; i++) {
         z = ( REAL(x1)[i] - REAL(x2)[i] ) / REAL(par)[0];
         SS2 += z * z; 
       }

       d2 = exp(-SS2);
       rkern[0] = REAL(par)[1] * d2;
       rdkern[1] =  d2; 
       rdkern[0] =  2 * rkern[0] * SS2 / REAL(par)[0];

       SET_ATTR(kern, install(\"gradient\"), dkern);
       UNPROTECT(nprotect);
       return kern;
     "

## inline the C function into an R function: MUCH MORE EFFICIENT!!!
\dontrun{
require(inline)
kernC <- cfunction(sig = signature(x1 = "numeric", x2 = "numeric",
                                   par = "numeric"),
                    body = kernCode)
myCovMan <- covMan(kernel = kernC, hasGrad = TRUE, label = "myGauss", d = 4,
                   parNames = c("theta", "sigma2"),
                   parLower = c("theta" = 0.0, "sigma2" = 0.0),
                   parUpper = c("theta" = Inf, "sigma2" = Inf))
}

##=======================================================================
## Example (continued). Simulate data for covMan and trend
##=======================================================================
n <- 100; 
X <- matrix(runif(n * d), nrow = n)
colnames(X) <- inputNames(myCovMan)

coef(myCovMan) <- myPar <- c(theta = 0.5, sigma2 = 2)
C <- covMat(object = myCovMan, X = X,
            compGrad = FALSE,  index = 1L)

library(MASS)
set.seed(456)
y <- mvrnorm(mu = rep(0, n), Sigma = C)
p <- rpois(1, lambda = 4)
if (p > 0) {
  F <- matrix(runif(n * p), nrow = n, ncol = p)
  beta <- rnorm(p)
  y <- F \%*\% beta + y
} else F <- NULL
par <- parCovIni <- c("theta" = 0.6, "sigma2" = 4)

##=======================================================================
## Example (continued). ML estimation. Note the 'partrack' argument
##=======================================================================           
est <- mle(object = myCovMan,
           parCovIni = parCovIni,
           y = y, X = X, F = F,
           parCovLower = c(0.05, 0.05), parCovUpper = c(10, 100),
           parTrack = TRUE, noise = FALSE, checkNames = FALSE)
est$opt$value

## change the (constrained) optimization  method
\dontrun{
est1 <- mle(object = myCovMan,
            parCovIni = parCovIni,
            optimFun = "stats::optim",
            optimMethod = "L-BFGS-B",
            y = y, X = X, F = F,
            parCovLower = c(0.05, 0.05), parCovUpper = c(10, 100),
            parTrack = TRUE, noise = FALSE, checkNames = FALSE)
est1$opt$value
}

##=======================================================================
## Example (continued). Grid for graphical analysis
##=======================================================================
\dontrun{
    theta.grid <- seq(from = 0.1, to = 0.7, by = 0.2)
    sigma2.grid <- seq(from = 0.3, to = 6, by = 0.4)
    par.grid <- expand.grid(theta = theta.grid, sigma2 = sigma2.grid)
    ll <- apply(as.matrix(par.grid), 1, est$logLikFun)
    llmat <- matrix(ll, nrow = length(theta.grid),
                    ncol = length(sigma2.grid))
}                

##=======================================================================
## Example (continued). Explore the surface ?
##=======================================================================
\dontrun{
   require(rgl)
   persp3d(x = theta.grid, y = sigma2.grid, z = ll,
           xlab = "theta", ylab = "sigma2", zlab = "logLik",
           col = "SpringGreen3", alpha = 0.6)
}

##=======================================================================
## Example (continued). Draw a contour plot for the log-lik 
##                        and show iterates
##=======================================================================
\dontrun{
    contour(x = theta.grid, y = sigma2.grid, z = llmat,
            col = "SpringGreen3", xlab = "theta", ylab = "sigma2",
            main = "log-likelihood contours and iterates",
            xlim = range(theta.grid, est$parTracked[ , 1], na.rm = TRUE),
            ylim = range(sigma2.grid, est$parTracked[ , 2], na.rm = TRUE))
    abline(v = est$coef.kernel[1], h = est$coef.kernel[2], lty = "dotted")
    niter <- nrow(est$parTracked)
    points(est$parTracked[1:niter-1, ],
           col = "orangered", bg = "yellow", pch = 21, lwd = 2, type = "o")
    points(est$parTracked[niter, , drop = FALSE],
           col = "blue", bg = "blue", pch = 21, lwd = 2, type = "o", cex = 1.5)
    ann <- seq(from = 1, to = niter, by = 5)
    text(x = est$parTracked[ann, 1], y = est$parTracked[ann, 2],
         labels = ann - 1L, pos = 4, cex = 0.8, col = "orangered")
    points(x = myPar["theta"], y = myPar["sigma2"],
           bg = "Chartreuse3", col = "ForestGreen",
           pch = 22, lwd = 2, cex = 1.4)

    legend("topright", legend = c("optim", "optim (last)", "true"),
           pch = c(21, 21, 22), lwd = c(2, 2, 2), lty = c(1, 1, NA),
           col = c("orangered", "blue", "ForestGreen"),
           pt.bg = c("yellow", "blue", "Chartreuse3"))
 
}
}


