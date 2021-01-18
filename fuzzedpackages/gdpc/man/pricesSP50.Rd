\name{pricesSP50}
\alias{pricesSP50}
\docType{data}
\title{
Stock Prices of the First 50 Components of S&P500
}
\description{
Fifty series corresponding to the stock prices of the first 50 components of the Standard&Poor's 500 index. Five hundred daily observations starting 1/1/2010.
}
\usage{data(pricesSP50)}
\format{
  A matrix time series with 500 observations on the stock prices of the first 50 components of the Standard&Poor's 500 index.
}



\examples{
data(pricesSP50)
\dontrun{
#Plot the first four series
plot(pricesSP50[, 1:4], main = 'Four components of the S&P500 index')
#Compute GDPCs; this may take a bit.
fit_SP <- auto.gdpc(pricesSP50, normalize = 2, niter_max = 1000, ncores= 4)
fit_SP
#Get reconstruction and plot
recons <- fitted(fit_SP, num_comp = 2)
colnames(recons) <- colnames(pricesSP50)
plot(recons[, 1:4], main = 'Reconstruction of four components of the S&P500 index')
}
}
\keyword{datasets}