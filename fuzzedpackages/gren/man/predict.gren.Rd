\name{predict.gren}

\alias{predict.gren}
\alias{coef.gren}

\title{
Predictions and coefficients from gren model.}

\description{
Create predictions from new data using a fitted gren model/retrieve coefficients from fitted model. Both are S3 methods.}

\usage{
\method{predict}{gren}(object, newx, unpenalized=NULL, s=NULL, 
             type=c("groupreg", "regular"), ...)

\method{coef}{gren}(object, s=NULL, type=c("groupreg", "regular"), ...)
}

\arguments{
\item{object}{
A fitted \code{gren} model.}
\item{newx}{
New data for which to do predictions.}
\item{unpenalized}{
New unpenalized data for which to do predictions.}
\item{s}{
Value of \code{lambda} for which to create predictions/coefficients, may be a \code{vector}.}
\item{type}{
Either \code{groupreg}, which creates predictions/coefficients of group-regularized model, or \code{regular}, for predictions/coefficients from regular model.}
\item{...}{
Further arguments to be passed.}
}

\details{
This are the predict/coefficient functions of the \code{gren} package.}

\value{
\code{predict} returns a \code{numeric} \code{matrix} with predicted probabilities. \code{coef} returns a \code{matrix} with coefficients.}

\references{
Münch, M.M., Peeters, C.F.W., van der Vaart, A.W., and van de Wiel, M.A. (2018). Adaptive group-regularized logistic elastic net regression. arXiv:1805.00389v1 [stat.ME].
}

\author{
Magnus M. Münch <m.munch@vumc.nl>}

\seealso{
\code{\link{gren}},\code{\link{coef.gren}}}

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

## create new data
xnew <- matrix(rnorm(n*p), ncol=p, nrow=n)

## create predictions/coefficients
preds <- predict(fit.gren, xnew, type="groupreg")
coefs <- coef(fit.gren, type="groupreg")
}
}