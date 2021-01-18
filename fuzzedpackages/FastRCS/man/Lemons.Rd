\name{Lemons}
\alias{Lemons}
\title{Sales Data for the Chrysler Town \& Country}
\description{
 Sales data for the Chrysler Town \& Country. 
}
\usage{
Lemons
}
\format{
  \describe{
	\item{\code{VehBCost}}{Acquisition cost paid for the vehicle at time of purchase.}
	\item{\code{MMRAcquisitionAuctionAveragePrice}}{Acquisition price for this vehicle in average condition at time of purchase.}
	\item{\code{MMRAcquisitonRetailCleanPrice}}{Acquisition price for this vehicle in the above Average condition at time of purchase.}	
	\item{\code{MMRAcquisitionRetailAveragePrice}}{Acquisition price for this vehicle in the retail market in average condition at time of purchase.}
	\item{\code{MMRAcquisitonRetailCleanPrice}}{Acquisition price for this vehicle in the retail market in above average condition at time of purchase.}
	\item{\code{MMRCurrentAuctionAveragePrice}}{Acquisition price for this vehicle in average condition as of current day.}
	\item{\code{MMRCurrentAuctionCleanPrice}}{Acquisition price for this vehicle in above condition as of current day.}
	\item{\code{MMRCurrentRetailAveragePrice}}{Acquisition price for this vehicle on the retail market in average condition as of current day.}
	\item{\code{MMRCurrentRetailCleanPrice}}{Acquisition price for this vehicle on the retail market in above average condition as of current day.}
	\item{\code{WarrantyCost}}{Warranty price (term=36month  and millage=36K).}
	\item{\code{VehOdo}}{The vehicle's odometer reading.}
  }
}
\examples{
data(Lemons)
alpha<-0.5
p<-ncol(Lemons)
ns<-FRCSnumStarts(p=p,eps=(1-alpha)*4/5)
Fit<-FastRCS(x=Lemons[,-1],y=Lemons[,1],nSamp=ns,seed=1)
plot(Fit)
}
\keyword{datasets}
