\name{grid_to_spdf}
\alias{grid_to_spdf}
\docType{package}

\title{
Transforms Grid Object into SpatialPolygonsDataFrame

}

\description{
Transforms Grid Object (the output of the gwfa function) into SpatialPolygonsDataFrame.

}

\usage{
  grid_to_spdf(df,epsg,cell_size=NULL, bandwith = NULL, radius = NULL, 
    q = NULL, sample_size = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{df}{the \code{data.frame} obtained from the gwfa function}
  \item{epsg}{EPSG code}
  \item{cell_size}{leave empty}
  \item{bandwith}{leave empty}
  \item{radius}{leave empty}
  \item{q}{leave empty}
  \item{sample_size}{leave empty}
  
  
  
 }
\details{
Returns a \code{spatialPolygonsdataframe}.
}  

\author{ Cecile Tannier, Stephane G. Roux and Francois Semecurbe }

\examples{
%\dontrun{
library(gwfa)

data("mariegalante")

test=gwfa(points=mariegalante,q=0,radius=(20*2^((0:6)/2)),
bandwith=1600,sample_size=500,cell_size=2000)
test=test[test$count>100,]#select the cells with at least 100 points. 


#estimate the fractal dimension on the 7 radius
X=cbind(rep(1,length(test@radius)),log2(test@radius))
fit_frac_dim=(do.call(cbind,test[,4:10]))\%*\%t(solve(t(X)\%*\%X)\%*\%t(X))
test$dimfrac=fit_frac_dim[,2]

#create spatial polygon dataframe
shp=grid_to_spdf(test,"2970")


%}
}
