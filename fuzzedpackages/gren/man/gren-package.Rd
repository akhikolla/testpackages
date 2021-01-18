\name{gren-package}
\alias{gren-package}
\docType{package}
\title{
Adaptive group-regularized logistic elastic net regression
}
\description{
The package allows the user to incorporate multiple sources of co-data (e.g., previously obtained p-values, published gene lists, and annotation) in the estimation of a logistic regression model to enhance predictive performance.
}
\details{
The main function of the package is \code{\link{gren}}, which estimates a group-regularized elastic net regression model. The following functions are convenience functions:

\describe{
\item{\code{\link{cv.gren}}}{estimates performance measures by efficient cross-validation.}
\item{\code{\link{coef.gren}}}{S3 method to retrieve model parameters from a \code{gren} fit.}
\item{\code{\link{predict.gren}}}{S3 method to get predictions for new data from a \code{gren} fit.}
\item{\code{\link{denet}}}{density function of the elastic net prior distribution.}
\item{\code{\link{renet}}}{generate samples from the elastic net prior distribution.}
}
}

\author{
Magnus M. Münch
Maintainer: Magnus M. Münch <m.munch@vumc.nl>
}
\references{
Münch, M.M., Peeters, C.F.W., van der Vaart, A.W., and van de Wiel, M.A. (2018). Adaptive group-regularized logistic elastic net regression. arXiv:1805.00389v1 [stat.ME].
}
\keyword{ package }
\seealso{
\code{\link{cv.gren}}
}
\examples{
\donttest{
## Create data
p <- 1000
n <- 100
set.seed(2018)
x <- matrix(rnorm(n*p), ncol=p, nrow=n)
beta <- c(rnorm(p/2, 0, 0.1), rnorm(p/2, 0, 1))
m <- rep(1, n)
y <- rbinom(n, m, as.numeric(1/(1 + exp(-x \%*\% as.matrix(beta)))))
partitions <- list(groups=rep(c(1, 2), each=p/2))

## estimate model
fit.gren <- gren(x, y, m, partitions=partitions)
}
}
