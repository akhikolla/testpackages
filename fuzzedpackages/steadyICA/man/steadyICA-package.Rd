\name{steadyICA-package}
\alias{steadyICA-package}
\docType{package}
\title{
ICA via distance covariance, tests of mutual independence, and other ICA functions
}
\description{
Functions related to multivariate measures of independence and ICA: \cr
 -estimate independent components by minimizing distance covariance;\cr
 -conduct a test of mutual independence based on distance covariance; \cr
 -estimate independent components via infomax (a popular method but generally performs poorer than steadyICA or ProDenICA but is useful for comparisons);\cr
 -order independent components by skewness;\cr
 -match independent components from multiple estimates;\cr
 -other functions useful in ICA.
}
\details{
\tabular{ll}{
Package: \tab steadyICA\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2015-11-08\cr
License: \tab GPL (>= 2)\cr
Depends: \tab Rcpp (>= 0.9.13), MASS\cr
Suggests: \tab irlba, JADE, ProDenICA, fastICA
}
}
\author{
Benjamin B. Risk and Nicholas A. James and David S. Matteson.\cr
Maintainer: Benjamin Risk <bbr28@cornell.edu>
}

\references{
Bernaards, C. & Jennrich, R. (2005) Gradient projection algorithms and software for arbitrary rotation criteria in factor analysis. \emph{Educational and Psychological Measurement} 65, 676-696

Matteson, D. S. & Tsay, R. Independent component analysis via U-Statistics. 
    <http://www.stat.cornell.edu/~matteson/#ICA>

Szekely, G., Rizzo, M. & Bakirov, N. Measuring and testing dependence by correlation of distances. (2007) \emph{The Annals of Statistics}, 35, 2769-2794.

Tichavsky, P. & Koldovsky, Z. Optimal pairing of signal components separated by blind techniques. (2004) \emph{Signal Processing Letters} 11, 119-122.
}

\keyword{ ica }
\keyword{ infomax}
\keyword{ matching}
\keyword{ distance }
\keyword{ covariance}
\keyword{ whiten}
\keyword{dcovICA}

\seealso{
\code{\link[fastICA:fastICA]{fastICA}}
\code{\link[ProDenICA:ProDenICA]{ProDenICA::ProDenICA}}
}
\examples{
#see steadyICA
}
