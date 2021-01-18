\name{cv.gren}
\alias{cv.gren}
\title{
Performance cross-validation of group-regularized logistic elastic net regression
}
\description{
Function that cross-validations for performance estimation of \code{gren} models.
}
\usage{
cv.gren(x, y, m=rep(1, nrow(x)), unpenalized=NULL, partitions=NULL, alpha=0.5, 
        lambda=NULL, intercept=TRUE, monotone=NULL, psel=TRUE, compare=TRUE, 
        posterior=FALSE, nfolds=nrow(x), foldid=NULL, trace=TRUE,
        control=list(epsilon=0.001, maxit=500, maxit.opt=1000, maxit.vb=100),
        keep.pred=TRUE, fix.lambda=FALSE, nfolds.out=nrow(x), foldid.out=NULL,
        type.measure=c("auc", "deviance", "class.error"))
}

\arguments{
\item{x}{
See \code{\link{gren}}.}
\item{y}{
See \code{\link{gren}}.}
\item{m}{
See \code{\link{gren}}.}
\item{unpenalized}{
See \code{\link{gren}}.}
\item{partitions}{
See \code{\link{gren}}.}
\item{alpha}{
See \code{\link{gren}}.}
\item{lambda}{
See \code{\link{gren}}.}
\item{intercept}{
See \code{\link{gren}}.}
\item{monotone}{
See \code{\link{gren}}.}
\item{psel}{
See \code{\link{gren}}.}
\item{compare}{
See \code{\link{gren}}.}
\item{posterior}{
See \code{\link{gren}}.}
\item{nfolds}{
See \code{\link{gren}}.}
\item{foldid}{
See \code{\link{gren}}.}
\item{trace}{
if \code{TRUE}, progress of the cross-validation is printed.}
\item{control}{
See \code{\link{gren}}.}
\item{keep.pred}{
\code{logical}. If \code{TRUE} the cross-validated predictions are saved.}
\item{fix.lambda}{
\code{logical}. If \code{TRUE} \code{lambda} is cross-validated only once and used in every fold.}
\item{nfolds.out}{
\code{numeric} that gives the number of folds to use in the performance cross-validation. Default is \code{nrow(x)}.}
\item{foldid.out}{
optional \code{numeric} vector of length \code{nrow(x)} with the performance cross-validation fold assignments of the observations.}
\item{type.measure}{
\code{numeric} that gives the performance measures to calculate. Currently possible are AUC, deviance, and misclassification error.
}
}

\details{
\code{cv.gren} is a convenience function that gives cross-validated predictions. Performance measures are optionally calculated with these predictions. \code{cv.gren} is more efficient than simply looping over the folds, since it uses the final estimates of the previous fold as starting values for the next fold. This substantially reduces computation time. 
}

\value{
Function returns a \code{list} of length two with the following components:
\item{groupreg}{
\code{list} with the cross-validated predictions of length \code{nrow(x)} and possibly performance metrics of the group-regularized model.
}
\item{regular}{
\code{list} with the cross-validated predictions of length \code{nrow(x)} and possibly performance metrics of the regular model.
}
}

\references{
Münch, M.M., Peeters, C.F.W., van der Vaart, A.W., and van de Wiel, M.A. (2018). Adaptive group-regularized logistic elastic net regression. arXiv:1805.00389v1 [stat.ME].
}

\author{
Magnus M. Münch <m.munch@vumc.nl>
}

\seealso{
\code{\link{gren}}, \code{\link{predict.gren}}, \code{\link{coef.gren}}
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

## calculate cross-validated predictions and performance measures
fit.cv.gren <- cv.gren(x, y, m, partitions=partitions, fix.lambda=TRUE)
}
}