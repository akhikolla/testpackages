\name{fecr_probs}
\alias{fecr_probs}
\title{Compute probability of the reduction parameter relative to a given threshold}
\usage{
fecr_probs(stanFit, threshold = 0.95, lessthan = TRUE, 
           plot = TRUE, xlab, ylab, main, verbose = TRUE, ...)
}
\arguments{
\item{stanFit}{a \code{\link[rstan]{stanfit}} object from the output of \code{\link{fecr_stan}()}.}

\item{threshold}{numeric. The default threshold is 0.95 (95\%).}

\item{lessthan}{logical. If TRUE, the probability less than the threshold is computed. Otherwise greater or equal to the threshold is computed. Default is TRUE.}

\item{plot}{logical. If TRUE, the posterior density of the reduction is plotted with region less than the threshold shaded.}

\item{xlab, ylab, main}{strings. Arguments for plotting. Only used if \code{plot = TRUE}.}

\item{verbose}{logical. If TRUE, a statement with computed probability is printed.}

\item{...}{additional plotting arguments}
}
\value{
Returns a numeric value indicating the probability in percentage.
}
\description{
Computes probability of the reduction parameter's marginal posterior density relative to a threshold.
}
\author{
Craig Wang
}

\examples{
\donttest{
## load sample data
data(epgs)

## apply individual efficacy model to the data vectors
model <- fecr_stan(epgs$before, epgs$after, rawCounts = FALSE, preCF = 10, 
                   paired = TRUE, indEfficacy = TRUE)
fecr_probs(model$stan.samples)
}}

