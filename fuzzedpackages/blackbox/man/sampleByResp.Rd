\name{sampleByResp}
\alias{sampleByResp}
\title{
Sample predictor points according to predicted response
}
\description{
Assuming that \code{\link{calcPredictorOK}} and \code{\link{maximizeOK}} have been first run:
predictor points can be sampled in several ways: the convex hull of predictor points with predicted response higher than some threshold value can be sampled uniformly. An Expected Improvement (e.g. Bingham et al., 2014) strategy can be used; whereby points with the highest predicted probability of improvement of the response value among a set of candidates sampled uniformly are retained. An expanded convex hull allowing further exploration of predictor space can also be considered. This function performs various combinations of these methods and (if the response was treated as a likelihood surface) can further use information from any previous likelihood ratio test of confidence interval computations.
}
\usage{
sampleByResp(size = blackbox.getOption("nextPointNumber"), outfile = NULL, useEI,
             NextBoundsLevel = 0.001,
             threshold=qchisq(1-NextBoundsLevel, 1)/2,
             rnd.seed = NULL, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{size}{
 sample size
}
  \item{outfile}{
  If not NULL, the name of an ASCII file where to print the result as a table.
}
  \item{useEI}{
  Whether to use an expected improvement criterion
}
  \item{NextBoundsLevel}{
   Controls \code{threshold} in a way meaningful for log-likelihood surfaces
}
  \item{threshold}{
   Controls the threshold for selection of the vertices of the convex hull to be sampled, and for inclusion of candidate predictor points in the sample. This threshold corresponds to a difference between predicted value and maximum predicted value. The actual maximal difference for inclusion of vertices additionally depends on the residual error of the predictor.
}
  \item{rnd.seed}{
  NULL (in which case nothing is done) or an integer (in which case \code{set.seed(seed=rnd.seed)} is called).
}
  \item{verbose}{
  To print information about evaluation, for development purposes.
}
}
\details{
The sampling procedure is designed to balance exploration of new regions of the predictor space and filling the top of a likelihood surface, or accurately locating the maximum and bounds of one-dimensional profile likelihood confidence interval. Details are yet to be documented.
}
\value{
Returns the predictor points invisibly.
}
\references{
D. Bingham, P. Ranjan, and W.J. Welch (2014) Design of Computer Experiments for Optimization, Estimation of Function Contours, and Related Objectives, pp. 109-124 in Statistics in Action: A Canadian Outlook (J.F. Lawless, ed.). Chapman and Hall/CRC.
}
