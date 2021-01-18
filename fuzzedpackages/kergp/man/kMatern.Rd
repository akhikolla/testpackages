\name{kMatern}
\alias{kMatern}
\alias{kExp}
\encoding{UTF8}
\title{
  
  Matérn Kernels

}
\description{
  
  Matérn kernels, obtained by plugging the Euclidian norm into a
  1-dimensional Matérn function.
  
}

\usage{

   kMatern(d, nu = "5/2")

}

\arguments{
  \item{d}{

    Dimension.

  }
  \item{nu}{

    Character corresponding to the smoothness parameter \eqn{\nu} of
    Matérn kernels. At this stage, the possible values are "1/2"
    (exponential kernel), "3/2" or "5/2".  } }

\value{

  An object of class \code{"covMan"} with default parameters: 1 for
  ranges and variance values.

}

\references{

  C.E. Rasmussen and C.K.I. Williams (2006), \emph{Gaussian Processes
    for Machine Learning}, the MIT Press,
  \url{http://www.GaussianProcess.org/gpml}
  
}
\note{

  Notice that these kernels are NOT obtained by tensor product.

}
\examples{
kMatern()  # default: d = 1, nu = 5/2
kMatern(d = 2)
myMatern <- kMatern(d = 5, nu = "3/2")
coef(myMatern) <- c(range = 1:5, sigma2 = 0.1)
myMatern
try(kMatern(nu = 2))  # error
}
