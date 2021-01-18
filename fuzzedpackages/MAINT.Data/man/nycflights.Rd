\name{nycflights}
\docType{data}
\alias{nycflights}
\alias{FlightsDF}
\alias{FlightsUnits}
\alias{FlightsIdt}
\title{New York City flights Data Set}
\description{A interval-valued data set containing 142 units and four interval-valued variables (dep_delay, arr_delay, air_time and distance), created from from the flights data set in the R package nycflights13
(on-time data for all flights that departed the JFK, LGA or EWR airports in 2013), after removing all rows with missing observations, and aggregating by month and carrier.}
\usage{data(nycflights)}
\format{FlightsDF: A data frame containing the original 327346 valid (i.e. with non missing values) flights from the nycflights13 package, described by the 4 variables: dep_delay, arr_delay, air_time and distance.\cr
       FlightsUnits: A factor with 327346 observations and 142 levels, indicating the month by carrier combination to which each orginal flight belongs to.\cr
       FlightsIdt: An IData object with 142 observations and 4 interval-valued variables, describing the intervals formed by agregating the FlightsDF microdata by the 0.05 and 0.95 quantiles of the subsamples formed by FlightsUnits factor.}
\keyword{datasets}
\keyword{interval data}


