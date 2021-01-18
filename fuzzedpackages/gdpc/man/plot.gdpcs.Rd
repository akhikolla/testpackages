\name{plot.gdpcs}
\alias{plot.gdpcs}
\title{Plot Generalized Dynamic Principal Components}
\description{
Plots a \code{gdpcs} object.
}


\usage{
\method{plot}{gdpcs}(x, which_comp = 1, plot.type = 'multiple', \dots)
}
\arguments{
  \item{x}{An object of class \code{gdpcs}, usually the result of \link{auto.gdpc}.}
  \item{which_comp}{Numeric vector indicating which components to plot. Default is 1.}
  \item{plot.type}{Argument to be passed to \code{\link[zoo]{plot.zoo}}. Used only when the original data set was stored in an object of class \code{zoo}. Default is 'multiple'.}
  \item{\dots}{Additional arguments to be passed to the plotting functions.}
}


\author{
Daniel Peña, Ezequiel Smucler, Victor Yohai
}

\seealso{
\code{\link{gdpc}}, \code{\link{auto.gdpc}}, \code{\link{plot.gdpc}}
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
#Choose number of lags using the LOO criterion.
#k_max=2 to keep computation time low
autofit <- auto.gdpc(x, k_max = 2, auto_comp = FALSE, num_comp = 2) 
autofit
plot(autofit, which_comp = c(1,2), xlab = '', ylab = '')
}

\keyword{ts}