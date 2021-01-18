\name{stepsel}
\alias{stepsel}
\alias{stepsel.MRC}
\alias{stepsel.AIC}
\alias{stepsel.BIC}

\title{Automatic selection of number of jumps}

\description{
Select the number of jumps.
}

\usage{
stepsel(path, y, type = c("MRC", "AIC", "BIC"), \dots)
stepsel.MRC(path, y, q, alpha = 0.05, r = ceiling(50 / min(alpha, 1 - alpha)),
  lengths = if(attr(path$cand, "family") == "gaussKern")
    2^(floor(log2(length(y))):ceiling(log2(length(attr(path$cand, "param")$kern)))) else
    2^(floor(log2(length(y))):0),
  penalty = c("none", "log", "sqrt"), name = if(attr(path$cand, "family") == "gaussKern")
    ".MRC.ktable" else ".MRC.table",
  pos = .GlobalEnv)
stepsel.AIC(path, \dots)
stepsel.BIC(path, \dots)
}

\arguments{
\item{path}{an object of class \code{\link{steppath}}}
\item{y}{for \code{type=MRC} only: a numeric vector containing the serial data}
\item{type}{how to select, dispatches specific method}
\item{\dots}{further argument passed to specific method}
\item{q, alpha, r, lengths, penalty, name, pos}{see \code{\link{bounds}}}
}

\value{
A single integer giving the number of blocks selected, with \code{\link{attr}}ibute \code{crit} containing the values of the criterion (MRC / AIC / BIC) for each fit in the path.
}

\note{
To obtain the threshold described in Boysen et al.~(2009, Theorem~5), set \code{q=(1+delta) * sdrobnorm(y) * sqrt(2*length(y))} for some positive \code{delta} and \code{penalty="none"}.
}

\seealso{\code{\link{steppath}}, \code{\link{stepfit}}, \link{family}, \code{\link{stepbound}}}

\references{
Boysen, L., Kempe, A., Liebscher, V., Munk, A., Wittich, O. (2009) Consistencies and rates of convergence of jump-penalized least squares estimators. \emph{The Annals of Statistics} \bold{37}(1), 157--183.

Yao, Y.-C. (1988) Estimating the number of change-points via Schwarz' criterion. \emph{Statistics & Probability Letters} \bold{6}, 181--189.
}

\examples{
# simulate 5 blocks (4 jumps) within a total of 100 data points
b <- c(sort(sample(1:99, 4)), 100)
f <- rep(rnorm(5, 0, 4), c(b[1], diff(b)))
rbind(b = b, f = unique(f))
# add gaussian noise
y <- f + rnorm(100)
# find 10 candidate jumps
path <- steppath(y, max.cand = 10)
# select number of jumps by simulated MRC with sqrt-penalty
# thresholded with positive delta, and by BIC
sel.MRC <- stepsel(path, y, "MRC", alpha = 0.05, r = 1e2, penalty = "sqrt")
sel.MRC
delta <- .1
sel.delta <- stepsel(path, y, "MRC",
  q = (1 + delta) * sdrobnorm(y) * sqrt(2 * length(y)), penalty = "none")
sel.delta
sel.BIC <- stepsel(path, type="BIC")
sel.BIC
# compare results with truth
fit.MRC <- path[[sel.MRC]]
as.data.frame(fit.MRC)
as.data.frame(path[[sel.delta]])
as.data.frame(path[[sel.BIC]])
}

\keyword{nonparametric}
