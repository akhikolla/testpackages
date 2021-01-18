\name{simulate, covAll-method}
\alias{simulate,covAll-method}
\title{

  Simulation of a \code{covAll} Object

}
\description{
  
    Simulation of a \code{covAll} object.
  
}
\usage{

\S4method{simulate}{covAll}(object, nsim = 1, seed = NULL,
         X, mu = NULL, method = "mvrnorm", checkNames = TRUE,
         ...)  

}
\arguments{
  \item{object}{

    A covariance kernel object.

  }
  \item{nsim}{

    Number of simulated paths.

  }
  \item{seed}{

    Not used yet.
    
  }
  \item{X}{

    A matrix with the needed inputs as its columns.
   
  }
  \item{mu}{

    Optional vector with length \code{nrow(X)} giving the expectation
    \eqn{\mu(\mathbf{x})}{\mu(x)} of the Gaussian Process at the
    simulation locations \eqn{\mathbf{x}}{x}.
   
  }
  \item{method}{

    Character used to choose the simulation method. For now the only
    possible value is \code{"mvrnorm"} corresponding to the function
    with this name in the \bold{MASS} package.
   
  }
  \item{checkNames}{

    Logical. It \code{TRUE} the colnames of \code{X} and the input names
    of \code{object} as given by \code{inputNames(object)} must be
    identical sets.
   
  }  
  \item{\dots}{

    Other arguments for methods.

  }
}
\value{

  A numeric matrix with \code{nrow(X)} rows and \code{nsim} columns.
  Each column is the vector of the simulated path at the simulation
  locations.

}
\note{

  The simulation is unconditional.
  
}
\seealso{

  The \code{\link{mvrnorm}} function.

}
\examples{

## -- as in example(kergp) define an argumentwise invariant kernel --

kernFun <- function(x1, x2, par) {
  h <- (abs(x1) - abs(x2)) / par[1]
  S <- sum(h^2)
  d2 <- exp(-S)
  K <- par[2] * d2
  d1 <- 2 * K * S / par[1]   
  attr(K, "gradient") <- c(theta = d1,  sigma2 = d2)
  return(K)
}

covSymGauss <- covMan(kernel = kernFun,
                      hasGrad = TRUE,
                      label = "argumentwise invariant",
                      d = 2,
                      parNames = c("theta", "sigma2"),
                      par = c(theta = 0.5, sigma2 = 2))

## -- simulate a path from the corresponding GP --

nGrid <- 24; n <- nGrid^2; d <- 2
xGrid <- seq(from = -1, to = 1, length.out = nGrid)
Xgrid <- expand.grid(x1 = xGrid, x2 = xGrid)

ySim <- simulate(covSymGauss, X = Xgrid)
contour(x = xGrid, y = xGrid,
        z = matrix(ySim, nrow = nGrid, ncol = nGrid), 
        nlevels = 15)
}