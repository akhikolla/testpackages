\name{simulate.gp}
\alias{simulate.gp}
\title{

  Simulation of Paths from a \code{gp} Object
}
\usage{
\method{simulate}{gp}(object, nsim = 1L, seed = NULL,
         newdata = NULL,
         cond = TRUE,
         trendKnown = FALSE,
         newVarNoise = NULL,
         checkNames = TRUE,
         output = c("list", "matrix"),
         label = "y", unit = "",
         ...)
}
\arguments{
  \item{object}{

    An object with class \code{"gp"}.

  }
  \item{nsim}{

    Number of paths wanted.

  }
  \item{seed}{

    Not used yet.

  }
  \item{newdata}{

    A data frame containing the inputs values used for simulation as
    well as the required trend covariates, if any. This is similar to
    the \code{newdata} formal in \code{\link{predict.gp}}.
    
  }

  \item{cond}{

    Logical. Should the simulations be conditional on the observations
    used in the object or not?
    
  }
  \item{trendKnown}{

    Logical. If \code{TRUE} the vector of trend coefficients will be
    regarded as known so all simulated paths share the same trend. When
    \code{FALSE}, the trend must have been estimated so that its estimation
    covariance is known. Then each path will have a different trend.
    
  }
  \item{newVarNoise}{
    
    Variance of the noise for the "new" simulated observations. For the
    default \code{NULL}, the noise variance found in \code{object} is
    used. Note that if a very small positive value is used, each
    simulated path is the sum of the trend the smooth GP part and an
    interval containing say \eqn{95}\% of the simulated responses can
    be regarded as a confidence interval rather than a prediction
    interval.
    
  }
  \item{checkNames}{
    
    Logical. It \code{TRUE} the colnames of \code{X} and the input names
    of the covariance in \code{object} as given by
    \code{inputNames(object)} must be identical sets.
    
  }  
  \item{output}{

    The type of output wanted. A simple matrix as in standard simulation
    methods may be quite poor, since interesting intermediate results
    are then lost.
    
  }
  \item{label, unit}{

    A label and unit that will be copied into the output object
    when \code{output} is \code{"list"}.
    
  }
  \item{\dots}{

    Further arguments to be passed to the \code{simulate} method
    of the \code{"covAll"} class.
    
  }
  
}
\value{

  A matrix with the simulated paths as its columns or a more complete
  list with more results. This list which is given the S3 class
  \code{"simulate.gp"} has the following elements.
  
  \itemize{ \item{\code{X}, \code{F}, \code{y}}{ Inputs,
      trend covariates and response. } \item{\code{XNew},
      \code{FNew}}{ New inputs, new trend covariates. }
    \item{\code{sim}}{ Matrix of simulated paths. }
    \item{\code{trend}}{ Matrix of simulated trends. }
    \item{\code{trendKnown}, \code{noise}, \code{newVarNoise}}{
      Values of the formals. } \item{\code{Call}}{ The call. } }
}
\description{

  Simulation of paths from a \code{gp} object.

}
\note{
  
  When \code{betaKnown} is \code{FALSE}, the \emph{trend} and the
  \emph{smooth GP} parts of a simulation are usually correlated, and
  their sum will show less dispersion than each of the two
  components. The covariance of the vector
  \eqn{\widehat{\boldsymbol{\beta}}}{\beta} can be regarded as the
  posterior distribution corresponding to a non-informative prior, the
  distribution from which a new path is drawn being the predictive
  distribution.
  
}
\examples{
set.seed(314159)
n <- 40
x <- sort(runif(n))
y <- 2 + 4 * x  + 2 * x^2 + 3 * sin(6 * pi * x ) + 1.0 * rnorm(n)
df <- data.frame(x = x, y = y)

##-------------------------------------------------------------------------
## use a Matern 3/2 covariance. With model #2, the trend is mispecified,
## so the smooth GP part of the model captures a part of the trend.
##-------------------------------------------------------------------------

myKern <- k1Matern3_2
inputNames(myKern) <- "x"
mygp <- list()
mygp[[1]] <- gp(formula = y ~ x + I(x^2) + sin(6 * pi * x), data = df, 
                parCovLower = c(0.01, 0.01), parCovUpper = c(10, 100),
                cov = myKern, estim = TRUE, noise = TRUE)
mygp[[2]] <- gp(formula = y ~ sin(6 * pi * x), data = df, 
                parCovLower = c(0.01, 0.01), parCovUpper = c(10, 100),
                cov = myKern, estim = TRUE, noise = TRUE)

##-------------------------------------------------------------------------
## New data
##-------------------------------------------------------------------------

nNew <- 150
xNew <- seq(from = -0.2, to= 1.2, length.out = nNew)
dfNew <- data.frame(x = xNew)

opar <- par(mfrow = c(2L, 2L))

nsim <- 40
for (i in 1:2) {

    ##--------------------------------------------------------------------
    ## beta known or not, conditional
    ##--------------------------------------------------------------------

    simTU <- simulate(object = mygp[[i]], newdata = dfNew,  nsim = nsim,
                      trendKnown = FALSE)
    plot(simTU, main = "trend unknown, conditional")

    simTK <- simulate(object = mygp[[i]], newdata = dfNew, nsim = nsim,
                      trendKnown = TRUE)
    plot(simTK, main = "trend known, conditional")

    ##--------------------------------------------------------------------
    ## The same but UNconditional
    ##--------------------------------------------------------------------

    simTU <- simulate(object = mygp[[i]], newdata = dfNew,  nsim = nsim,
                     trendKnown = FALSE, cond = FALSE)
    plot(simTU, main = "trend unknown, unconditional")
    simTK <- simulate(object = mygp[[i]], newdata = dfNew, nsim = nsim,
                      trendKnown = TRUE, cond = FALSE)
    plot(simTK, main = "trend known, unconditional")
}

par(opar)
}
\author{
Yves Deville
}
