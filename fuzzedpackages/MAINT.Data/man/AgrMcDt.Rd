\name{AgrMcDt}
\alias{AgrMcDt}
\title{Agregate Micro Data}
\description{
AgrMcDt creates \code{\linkS4class{IData}} objects by agregating a Data Frame of Micro Data.
}

\usage{
AgrMcDt(MicDtDF, agrby, agrcrt="minmax")
}

\arguments{
  \item{MicDtDF}{A data frame with the original values of the micro data.}
  \item{agrby}{A factor with categories on which the micro data should be aggregated.}
  \item{agrcrt}{The aggregation criterion. Either the \sQuote{minmax} string, or a two dimensional vector with the prob. value for the left (lower) percentile, followed by the prob. value for the right (upper) percentile, used in the aggregation.}
}

\value{An object of class \code{\linkS4class{IData}} with the data set of Interval-valued variables resulting from the aggregation performed.}

\seealso{
  \code{\linkS4class{IData}} 
}

\keyword{interval data}

\examples{

# Create an Interval-Data object by agregating the microdata consisting 
# of 336776 NYC flights included in the FlightsDF data frame, 
# by the statistical units specified in the FlightsUnits factor.

Flightsminmax <- AgrMcDt(FlightsDF,FlightsUnits)

#Display the first and last observations

head(Flightsminmax)
tail(Flightsminmax)

#Print summary statistics

summary(Flightsminmax)

# Repeat this procedure using now the 10th and 90th percentiles.

Flights1090prcnt <- AgrMcDt(FlightsDF,FlightsUnits,agrcrt=c(0.1,0.9))

#Display the first and last observations

head(Flights1090prcnt)
tail(Flights1090prcnt)

summary(Flights1090prcnt)

}

