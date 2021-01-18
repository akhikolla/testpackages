\name{plot.gdpc}
\alias{plot.gdpc}
\title{Plot Generalized Dynamic Principal Components}
\description{
Plots a \code{gdpc} object.
}


\usage{
\method{plot}{gdpc}(x, which = 'Component', which_load = 0, \dots)
}
\arguments{
  \item{x}{An object of class \code{gdpc}, usually the result of \link{gdpc} or one of the entries of the result of \link{auto.gdpc}.}
  \item{which}{String. Indicates what to plot, either 'Component' or 'Loadings'. Default is 'Component'.}
  \item{which_load}{Lag number indicating which loadings should be plotted. Only used if which = 'Loadings'. Default is 0.}
  \item{\dots}{Additional arguments to be passed to the plotting functions.}
}


\author{
Daniel Peña, Ezequiel Smucler, Victor Yohai
}



\seealso{
\code{\link{gdpc}}, \code{\link{auto.gdpc}}, \code{\link{plot.gdpcs}}
}

\examples{
T <- 200 #length of series
m <- 200 #number of series
set.seed(1234)
f <- rnorm(T + 1)
x <- matrix(0, T, m)
u <- matrix(rnorm(T * m), T, m)
for (i in 1:m) {
    x[, i] <- 10 * sin(2 * pi * (i/m)) * f[1:T] + 10 * cos(2 * pi * (i/m)) * f[2:(T + 1)] + u[, i]
}
#Choose number of lags using the LOO type criterion.
#k_max=3 to keep computation time low
autofit <- auto.gdpc(x, k_max = 3) 
plot(autofit[[1]], xlab = '', ylab = '')
}

\keyword{ts}