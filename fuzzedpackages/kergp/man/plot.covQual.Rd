\name{plot}

\alias{plot}
\alias{plot.covQual}
\alias{plot,covQual-method}
\alias{plot,covQual,ANY-method} 
\title{
Plot for a qualitative input
}
\description{
Plots of the covariance matrix or the correlation matrix of a qualitative input. 
For an ordinal factor, the warping function can also be plotted. 
}
\usage{
\S4method{plot}{covQual}(x, y, type = c("cov", "cor", "warping"), ...)
}
\arguments{
  \item{x}{ An object of class \code{\link{covQual-class}}.}
  \item{y}{ Not used.}
  \item{type}{ A character indicating the desired type of plot. Type \code{warping} only works for an ordinal input.}
  \item{...}{ Other arguments passed to \code{corrplot::corrplot} or \code{plot}.}
}
\details{
 Covariance / correlation plots are done with package \code{corrplot} if loaded, or \code{lattice} else.
}
\seealso{
 \code{\link{covOrd}}.
}
\examples{
u <- ordered(1:6, levels = letters[1:6])

myCov2 <- covOrd(ordered = u, k1Fun1 = k1Fun1Cos, warpFun = "norm")
coef(myCov2) <- c(mean = 0.5, sd = 0.05, theta = 0.1) 

plot(myCov2, type = "cor", method = "ellipse")
plot(myCov2, type = "warp", col = "blue", lwd = 2)
}
