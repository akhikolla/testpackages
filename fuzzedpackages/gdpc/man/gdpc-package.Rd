\name{gdpc-package}
\alias{gdpc-package}
\docType{package}

\title{Generalized Dynamic Principal Components}

\description{
Computes the Generalized Dynamic Principal Components proposed in Peña and Yohai (2016).
}
\details{
\tabular{ll}{
Package: \tab gdpc\cr
Type: \tab Package\cr
Version: \tab 1.1.0\cr
Date: \tab 2018-05-08\cr
Depends: \tab R (>= 3.3.0)\cr
License: \tab GPL (>= 2)\cr
Imports: \tab xts, zoo, methods, Rcpp (>= 0.12.7), parallel, doParallel, foreach\cr
LinkingTo: \tab Rcpp, RcppArmadillo (>= 0.7.500.0.0)\cr
Suggests: \tab testthat\cr
NeedsCompilation: \tab yes\cr
}

Index:
\preformatted{
auto.gdpc   Automatic Fitting of Generalized Dynamic Principal Components.
}
\preformatted{
components.gdpcs    Get Generalized Dynamic Principal Components from a gdpcs object.
}
\preformatted{
fitted.gdpcs    Get reconstructed time series from a gdpcs object.
}
\preformatted{
gdpc    Computes a single Generalized Dynamic Principal Component
        with a given number of lags.
}
\preformatted{
ipi91   Six series corresponding to the Industrial Production
        Index (IPI) of France, Germany, Italy, United Kingdom, USA and 
        Japan. Monthly data from January 1991 to December 2012.
}
\preformatted{
plot.gdpc   Plots a gdpc object.
}
\preformatted{
plot.gdpcs   Plots a gdpcs object.
}
\preformatted{
pricesSP50   Fifty series corresponding to the stock prices of the first 50
             components of the Standard&Poor's 500 index. Five hundred daily
             observations starting 1/1/2010.
}
}
\author{
Daniel Peña, Ezequiel Smucler, Victor Yohai

Maintainer: Ezequiel Smucler <ezequiels.90@gmail.com>
}

\references{
Peña D. and Yohai V.J. (2016). “Generalized Dynamic Principal Components.” Journal of the
American Statistical Association, 111(515), 1121–1131.

Peña D., Smucler E. and Yohai V.J. (2020). “{gdpc}: An {R} Package for Generalized Dynamic Principal Components.” Journal of Statistical Software, to appear.
}


\examples{
data(ipi91)
\dontrun{
  #Compute GDPC, number of components and number of lags is chosen automatically.
  #This might take a bit.
  ipi_autogdpc <- auto.gdpc(ipi91)
  }
}
\keyword{ts}
\keyword{package}