\name{gdpc}
\alias{gdpc}
\title{Generalized Dynamic Principal Components}
\description{
Computes a single Generalized Dynamic Principal Component with a given number of lags.
}


\usage{gdpc(Z, k, f_ini = NULL, tol = 1e-4, niter_max = 500, crit = 'LOO')}
\arguments{
  \item{Z}{Data matrix. Each column is a different time series.}
  \item{k}{Integer. Number of lags to use.}
  \item{f_ini}{(Optional). Numeric vector. Starting point for the iterations. If no argument is passed the ordinary (non-dynamic) first principal component completed with k lags is used.}
  \item{tol}{Relative precision. Default is 1e-4.}
  \item{niter_max}{Integer. Maximum number of iterations. Default is 500.}
  \item{crit}{A string specifying the criterion to be used to evaluate the fitted model. Options are 'LOO', 'AIC', 'BIC' and 'BNG'. Default is 'LOO'.}
}

\details{

See \code{\link{auto.gdpc}} for the definition of criterion that is part of the output of this function.

}

\value{
  An object of class \code{gdpc}, that is, a list with entries:
   \item{expart}{Proportion of the variance explained.}
 \item{mse}{Mean squared error.}
 \item{crit}{The value of the criterion of the reconstruction, according to what the user specified.}
 \item{k}{Number of lags used.}
 \item{alpha}{Vector of intercepts corresponding to f.}
   \item{beta}{Matrix of loadings corresponding to f. Column number \eqn{k} is the vector of \eqn{k-1} lag loadings.}
 \item{f}{Coordinates of the first dynamic principal component corresponding to the periods \eqn{1,\dots,T}.}
  \item{initial_f}{Coordinates of the first dynamic principal component corresponding to the periods \eqn{-k+1,\dots,0}. Only for the case \eqn{k>0}, otherwise 0.}
       \item{call}{The matched call.}
       \item{conv}{Logical. Did the iterations converge?}
       \item{niter}{Integer. Number of iterations.}
 \code{fitted}, \code{plot} and \code{print} methods are available for this class.
 }
\author{
Daniel Peña, Ezequiel Smucler, Victor Yohai
}


\seealso{
\code{\link{auto.gdpc}}, \code{\link{plot.gdpc}}
}

\examples{
T <- 200 #length of series
m <- 500 #number of series
set.seed(1234)
f <- rnorm(T + 1)
x <- matrix(0, T, m)
u <- matrix(rnorm(T * m), T, m)
for (i in 1:m) {
    x[, i] <- 10 * sin(2 * pi * (i/m)) * f[1:T] + 10 * cos(2 * pi * (i/m)) * f[2:(T + 1)] + u[, i]
}
fit <- gdpc(x, k = 1) #find first DPC with one lag
fit
par(mfrow = c(1, 2)) #plot loadings
plot(fit, which = 'Loadings', which_load = 0, xlab = '', ylab = '') 
plot(fit, which = 'Loadings', which_load = 1, xlab = '', ylab = '') 
}

\keyword{ts}
