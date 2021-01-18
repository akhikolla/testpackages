\name{BCEE-package}
\alias{BCEE-package}
\alias{BCEE}
\docType{package}
\title{
The Bayesian Causal Effect Estimation (BCEE) Algorithm}
\description{
A Bayesian model averaging approach to causal effect estimation
based on the BCEE algorithm. Currently supports binary or continuous
exposures and outcomes. For more details, see: 
Talbot et al. (2015) DOI:10.1515/jci-2014-0035,
Talbot and Beaudoin (2020) arXiv:2003.11588.  
}
\details{
\tabular{ll}{
Package: \tab BCEE\cr
Type: \tab Package\cr
Version: \tab 1.3.0\cr
Date: \tab 2020-04-01\cr
License: \tab GPL (>=2) \cr
}
ABCEE(X, Y, U, omega), \cr
NBCEE(X, Y, U, omega), \cr
GBCEE(X, Y, U, omega)
}
\author{
Denis Talbot, Genevieve Lefebvre, Juli Atherton, Yohann Chiu.

Maintainer: Denis Talbot \email{denis.talbot@fmed.ulaval.ca}
}
\references{
Talbot, D., Lefebvre, G., Atherton, J. (2015) \emph{The Bayesian causal effect estimation algorithm}, Journal of Causal Inference, 3(2), 207-236.\\
Talbot, D., Beaudoin, C (2020) \emph{A generalized double robust Bayesian model averaging approach to causal effect estimation with application to the Study of Osteoporotic Fractures} arXiv:2003.11588 
}
\keyword{causal}
\keyword{confounding}
\keyword{model average}
\seealso{
 \code{\link{ABCEE}}, \code{\link{NBCEE}}, \code{\link{GBCEE}}.
}
