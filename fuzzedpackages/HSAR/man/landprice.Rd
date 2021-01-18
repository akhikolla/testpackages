\name{landprice}
\alias{landprice}
\docType{data}
\title{Leased residential land parcels, from 2003 to 2009 in Beijing, China}
\description{
  The residential land parcel data leased to real estate developers from the government during 2003 to 2009. The data contains price information and a range of locational and neighbourhood characteristics for each land parcel. There are 1117 samples in the data after dropping those without price information.     

}
\usage{data(landprice)}
\format{
  A data frame with 1117 observations on the following 11 variables.
  \describe{
    \item{obs}{An unique identifier for each land parcel.}
    \item{lnprice }{The log of the leasing price per square metre of each residential land parcel (unit: RMB, Chinese yuan) }
    \item{dsubway }{The log of the distance of each land parcel to the nearest railway station (unit:meters)}
    \item{dele }{The log of the distance of each land parcel to the nearest elementary school (unit:meters) }
    \item{dpark }{The log of the distance of each land parcel to the nearest green park (unit:meters) }
    \item{lnarea}{The log of the size of each land parcel (unit: square meters).}
    \item{lndcbd}{The log of the distance of each land parcel to the CBD (centre business district) in Beijing (unit:meters) }
    \item{year}{The year when each land parcel was leased with values of 0,1,2,3,4,5,6 representing year 2003,2004,2005,2006,2007,2008,2009}
    \item{popden}{The population density of each district (unit: 1000 persons per square kilometers) }
    \item{crimerate}{The number of reported serious crimes committed in each district per 1000 persons.}
    \item{district.id}{The identifier of the district where each land parcel is located.}
    }
  }

\seealso{ \code{\link{Beijingdistricts}}, \code{\link{landSPDF}}}

\references{
Harris, R., G.P. Dong, and W.Z. Zhang. (2013).Using Contextualised Geographically Weighted Regression to model the Spatial Heterogeneity of Land Prices in Beijing, China. \emph{Transaction in GIS 17(6)}, 901-19.


Dong G, Harris R, Jones K, Yu J (2015) Multilevel Modelling with Spatial Interaction Effects with Application to an Emerging Land Market in Beijing, China. \emph{PLoS ONE} 10(6): e0130761. doi:10.1371/journal.pone.0130761

}

\examples{
data(landprice)
head(landprice)

}
\keyword{hierarchical spatial data}
