\docType{data}
\name{dataf.population}
\alias{dataf.population}
\alias{population}
\title{World Historical Population-by-Country Dataset}
\format{
  The functional data as a data structure.
  \describe{
    \item{\code{dataf}}{
      The functional data as a list of objects. Each object is characterized by two coordinates 
      \describe{
        \item{\code{args}}{\bold{year} - a numeric vector of years 1950-2015 (66 years)}
        \item{\code{vals}}{\bold{population} - a numeric vector of the estimated total population in thousands in 233 countries and regions}
       }
    }
    \item{\code{labels}}{The geographic region of the country: Africa, Asia, Europe, Latin America, North America, Oceania}
    \item{\code{identifier}}{The name of country or region}
  }
}
\source{
United Nations, Department of Economic and Social Affairs,
Population Division,
\url{https://esa.un.org/unpd/wpp/Download/Standard/Population/}, 
file \code{Total population - Both sexes}
}
\usage{
dataf.population()
}
\description{
Historical world population data by countries.
}
\details{
World population data by a country, area or region as of 1 July 
of the year indicated. Figures are presented in thousands.
}
\examples{
## load the Population dataset
dataf = dataf.population()

## view the classes
unique(dataf$labels)

## access the 5th point of the 2nd object
dataf$dataf[[2]]$args[5]
dataf$dataf[[2]]$vals[5]

## plot the data
## Not run:
labels = unlist(dataf$labels)
plot(dataf,
  main = "World population data",
  xlab="Year", ylab="Population (in thousands)"
  )
## End(Not run)

## compute the integrated and infimal depths of the data curves
## with respect to the same set of curves
depthf.fd1(dataf$dataf, dataf$dataf)
}
\seealso{
\code{\link{dataf.population2010}}

\code{\link{dataf.*}} for other functional data sets

\code{\link{plot.functional}} for building plots of functional data
}
\keyword{datasets}
\keyword{functional}
