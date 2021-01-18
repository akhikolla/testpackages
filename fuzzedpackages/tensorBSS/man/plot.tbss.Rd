\name{plot.tbss}
\alias{plot.tbss}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot an Object of the Class tbss}
\description{
Plots the most interesting components (in the sense of extreme kurtosis) obtained by a tensor blind source separation method.
}

\usage{
\method{plot}{tbss}(x, first = 2, last = 2, datatype = NULL, 
     main = "The components with most extreme kurtoses", ...)
}
\arguments{
  \item{x}{Object of class tbss.}
  \item{first}{Number of components with maximal kurtosis to be selected.\cr See \code{\link{selectComponents}} for details.}
  \item{last}{Number of components with minimal kurtosis to be selected.\cr See \code{\link{selectComponents}} for details.}
  \item{main}{The title of the plot.}
  \item{datatype}{Parameter for choosing the type of plot, either \code{NULL}, \code{"iid"} or \code{"ts"}. The default \code{NULL} means the value from the tbss object \code{x} is taken.}
  \item{\dots}{Further arguments to be passed to the plotting functions, see details.}
}
\details{
The function \code{plot.tbss} first selects the most interesting components using \code{\link{selectComponents}} and then plots them either as a matrix of scatter plots using \code{\link{pairs}} (\code{datatype} = "iid") or as a time series plot using \code{\link[stats]{plot.ts}} (\code{datatype} = "ts"). 
Note that for \code{\link{tSOBI}} this criterion might not necessarily be meaningful as the method is based on second moments only.
}
\author{
Joni Virta
}

\examples{
data(zip.train)
x <- zip.train

rows <- which(x[, 1] == 0 | x[, 1] == 1)
x0 <- x[rows, 2:257]
y0 <- x[rows, 1] + 1

x0 <- t(x0)
dim(x0) <- c(16, 16, 2199)

tfobi <- tFOBI(x0)
plot(tfobi, col=y0)

if(require("stochvol")){
  n <- 1000
  S <- t(cbind(svsim(n, mu = -10, phi = 0.98, sigma = 0.2, nu = Inf)$y,
               svsim(n, mu = -5, phi = -0.98, sigma = 0.2, nu = 10)$y,
               svsim(n, mu = -10, phi = 0.70, sigma = 0.7, nu = Inf)$y,
               svsim(n, mu = -5, phi = -0.70, sigma = 0.7, nu = 10)$y,
               svsim(n, mu = -9, phi = 0.20, sigma = 0.01, nu = Inf)$y,
               svsim(n, mu = -9, phi = -0.20, sigma = 0.01, nu = 10)$y))
  dim(S) <- c(3, 2, n)

  A1 <- matrix(rnorm(9), 3, 3)
  A2 <- matrix(rnorm(4), 2, 2)

  X <- tensorTransform(S, A1, 1)
  X <- tensorTransform(X, A2, 2)

  tgfobi <- tgFOBI(X)
  plot(tgfobi, 1, 1)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }

