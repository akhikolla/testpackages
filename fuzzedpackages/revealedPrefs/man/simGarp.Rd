\name{simGarp}
\alias{simGarp}
\alias{simWarp}
\alias{simSarp}

\title{
Generate random data consistent with rationality axioms (WARP, SARP, GARP).
}
\description{
Functions for generating random data (prices and quantities) consistent with the chosen rationality axiom.
}
\usage{
simWarp(nobs, ngoods, afriat.par= 1, maxit= 10 * nobs, 
        qmin= 0, qmax= 1, pmin= 0, pmax= 1)
simSarp(nobs, ngoods, afriat.par= 1, maxit= 10 * nobs, 
        qmin= 0, qmax= 1, pmin= 0, pmax= 1)
simGarp(nobs, ngoods, afriat.par= 1, maxit= 10 * nobs, 
        qmin= 0, qmax= 1, pmin= 0, pmax= 1)
}
\arguments{
  \item{nobs}{the desired number of observations (number of rows in the quantities and prices matrices),}
  \item{ngoods}{the number of goods in the dataset (number of columns in the quantities and prices matrices),}
  \item{afriat.par}{the Afriat parameter, a real number in [0,1], which allows a certain level of error in the optimization of choices; default is 1, ie. no optimization error allowed,}
  \item{maxit}{maximum number of iterations (default to 10 times \code{nobs}),}
  \item{qmin}{minimum quantities for each good,}
  \item{qmax}{maximum quantities for each good,}
  \item{pmin}{minimum prices for each good,}
  \item{pmax}{maximum prices for each good.}
}
\value{
  \item{\code{x}}{numeric matrix of generated quantities,}
  \item{\code{p}}{numeric matrix of generated prices,} 
  \item{\code{iter}}{number of iterations before the algorithm stopped,}
  \item{\code{nobs}}{number of generated observations.}
}
\details{
The data are iteratively incremented: at each iteration a new random observation (prices and quantities) is generated, and is accepted only if it is consistent with the previously accepted data, in which case it is added to the data. The random observations (price-quantities couples) are independently generated from uniform distributions in the support defined by \code{qmin}, \code{qmax}, and \code{pmin}, \code{pmax}.

For GARP and SARP the depth-first search method is used to check for consistency (a recursive search using only the new candidate observation as starting point), for WARP the candidate observation is pairwise checked against all previously accepted data.

The algorithm stops if the desired number of observations \code{nobs} is reached. If the desired number of observations \code{nobs} is not reached in \code{maxit} iterations, a warning is issued and the function returns the largest dataset attained.
}
\references{
Varian, H. R. (1982) The Nonparametric Approach to Demand Analysis, 
\emph{Econometrica}, 50(4):945-973.

Varian, H. R. (1984) \emph{Microeconomic Analysis}. New York/London:
Norton, 2nd edition, pp 141-143.
}
\author{
Julien Boelaert \email{jubo.stats@gmail.com}
}
\seealso{
See \link{checkGarp} for rationality tests.
}
\examples{
# Generate 100 WARP-consistent observations of 5 goods
simdata <- simWarp(100, 5)
summary(checkWarp(simdata$x, simdata$p))
summary(checkGarp(simdata$x, simdata$p))
}
