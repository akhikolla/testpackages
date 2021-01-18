\name{compareBlocks}
\alias{compareBlocks}

\title{Compare fit blockwise with ground truth}
\description{Blockwise comparison of a fitted step function with a known ground truth using different criteria.}

\usage{
compareBlocks(truth, estimate, dist = 5e3)
}

\arguments{
\item{truth}{an object of class \code{\link{stepblock}} giving the ground truth, or a list of such objects}
\item{estimate}{corresponding estimated object(s) of class \code{\link{stepblock}}}
\item{dist}{a single \code{\link{numeric}} specifying the distance for at which jumps will be considered as having matched in the qualitative criterion}
}

\value{
A \code{\link{data.frame}}, containing just one row if two single \code{\link{stepblock}} were given, with columns
\item{true.num, est.num}{the true / estimated number of blocks}
\item{true.pos, false.pos, false.neg, sens.rate, prec.rate}{the number of true / false positive, false negatives, as well as the corresponding sensitivity and precision rates, where an estimated block is considered a true positive if it there is a corresponding block in the ground truth with both endpoints within \code{dist} of each other}
\item{fpsle}{false positive sensitive localization error: for each estimated block's midpoint find into which true block it falls, and sum distances of the respective borders}
\item{fnsle}{false negative sensitive localization error: for each true block's mid-point find into which estimated block it falls, and sum distances of the respective borders}
\item{total.le}{total localization error: sum of \code{fpsle} and \code{fnsle}}
}

\note{
No differences between true and fitted parameter \emph{values} are taking into account, only the precision of the detected blocks is considered; also, differing from the criteria in Elhaik et al.~(2010), no blocks are merged in the ground truth if its parameter values are close, as this may punish sensitive estimators. 

Beware that these criteria compare \emph{blockwise}, i.e. they do \emph{not} compare the precision of single jumps but for each block both endpoints have to match well at the same time.
}

\seealso{\code{\link{stepblock}}, \code{\link{stepfit}}, \code{\link{contMC}}}

\encoding{utf8}
\references{
Elhaik, E., Graur, D., \enc{Josić}{Josic}, K. (2010) Comparative testing of DNA segmentation algorithms using benchmark simulations. \emph{Molecular Biology and Evolution} \bold{27}(5), 1015-24.

Futschik, A., Hotz, T., Munk, A. Sieling, H. (2014) Multiresolution DNA partitioning: statistical evidence for segments. \emph{Bioinformatics},  \bold{30}(16), 2255--2262.
} 


\examples{
# simulate two Gaussian hidden Markov models of length 1000 with 2 states each
# with identical transition rates being 0.01 and 0.05, resp, signal-to-noise ratio is 5
sim <- lapply(c(0.01, 0.05), function(rate)
  contMC(1e3, 0:1, matrix(c(0, rate, rate, 0), 2), param=1/5))
plot(sim[[1]]$data)
lines(sim[[1]]$cont, col="red")
# use smuceR to estimate fit
fit <- lapply(sim, function(s) smuceR(s$data$y, s$data$x))
lines(fit[[1]], col="blue")
# compare fit with (discretised) ground truth
compareBlocks(lapply(sim, function(s) s$discr), fit)
}

\keyword{nonparametric}
