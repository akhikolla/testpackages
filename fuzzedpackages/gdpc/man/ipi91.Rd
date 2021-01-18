\name{ipi91}
\alias{ipi91}
\docType{data}
\title{
Industrial Production Index (IPI) of France, Germany, Italy, United Kingdom, USA and Japan
}
\description{
Six series corresponding to the Industrial Production Index (IPI) of France, Germany, Italy, United Kingdom, USA and Japan. Monthly data from January 1991 to December 2012.
}
\usage{data(ipi91)}
\format{
  A matrix time series with 264 observations on the following 6 variables.
 \describe{
    \item{\code{France}}{IPI of France.}
    \item{\code{Germany}}{IPI of Germany.}
    \item{\code{Italy}}{IPI of Italy.}
    \item{\code{United Kingdom}}{IPI of United Kingdom.}
    \item{\code{USA}}{IPI of USA.}
    \item{\code{Japan}}{IPI of Japan.}
  }
}



\examples{
data(ipi91)
plot(ipi91, plot.type = 'multiple', main = 'Industrial Production Index')
\dontrun{
#Compute first GDPC with nine lags; this may take a bit.
gdpc_ipi <- gdpc(ipi91, 9, niter_max = 1500)
#Plot the component
plot(gdpc_ipi, which = 'Component', ylab = '')
#Get reconstruction of the time series and plot
recons <- fitted(gdpc_ipi)
colnames(recons) <- colnames(ipi91)
plot(recons, main = 'Fitted values')
}
}
\keyword{datasets}