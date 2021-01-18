\name{gwfa}
\alias{gwfa}
\docType{package}

\title{
Geographically Weighted Fractal Analysis
}

\description{
Estimate local fractal dimension of a 2D set of points
}

\usage{
  gwfa(points, q = 0, radius, bandwith, sample_size, cell_size) 
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{points}{A data.frame containing two rows "x" and "y" corresponding to the longitude and latitude of the set of points analysed}
  
  \item{q}{The vector of order of the Renyi entropy. q=0 allows calculating more robusting box-counting dimension.}

  \item{radius}{the vector of radius of the multi-scale analysis}

  \item{bandwith}{Bandwith of the GWFA. This bandwidth acts as a local-global parameter. A large bandwidth leads to a global analysis. On the contrary, a small bandwidth leads to a local analysis. The unit of measurement is free. It must be the same as the unit of \code{cell_size} variable.}
  
  \item{sample_size}{size of the sample for calculate the average of the masses for the different radius and Renyi entropy orders \eqn{M^q(R)}. The choice of sample_size is a tradeoff between precision and computing speed. In our experience, 3000 is an appropriate choice.   }
  
  \item{cell_size}{Cell size of the grid. The unit of measurement is free. It must be the same as the unit of \code{bandwith} variable.}
 }

\details{
GWFA is a spatial method to analyse the variability of multiscale behavior in space of set of points. The vector of radius sets the scale of the multiscale analysis whereas the bandwith describes the size of the neighborhood of the analysis. The cell_size indicate the spacing between the estimate points. sample_size is a technique parameter indicates the size of the sample of the set of points used in the analysis. 

GWFA is a mix between Sandbox Multifractal analysis and Geographically Weighted Approach. The kernel used is the bisquare \eqn{(1-(\frac{distance}{bandwith})^2)^2} 

As output, we get a Grid object. A Grid class inherits from the class data.frame and it is completed by five slots (cell_size, bandwith, radius,q and sample_size) corresponding of input parameters of the gwfa function.

}


 
 
 
\references{
Vicsek, T. (1990). Mass multifractals. Physica A: Statistical Mechanics and its Applications, 168(1), 490-497.

Brunsdon, C., Fotheringham, S., & Charlton, M. (1998). Geographically weighted regression. Journal of the Royal Statistical Society: Series D (The Statistician), 47(3), 431-443.
}

\author{Cecile Tannier, Stephane G. Roux and Francois Semecurbe}

\examples{
\dontrun{
data("guadeloupe")

test=gwfa(points=guadeloupe,q=0,radius=(20*2^((0:6)/2)),
bandwith=1600,sample_size=3000,cell_size=2000)
test=test[test$count>100,]#select the cells with at least 100 points. 


#estimate the fractal dimension on the 7 radius
X=cbind(rep(1,length(test@radius)),log2(test@radius))
fit_frac_dim=(do.call(cbind,test[,4:10]))\%*\%t(solve(t(X)\%*\%X)\%*\%t(X))
test$dimfrac=fit_frac_dim[,2]

#create spatial polygon dataframe
shp=grid_to_spdf(test,"2970")

#convert to geographic information systems software format 
library(rgdal)
writeOGR(shp,"guade_analysis.shp","guade_analysis",driver="ESRI Shapefile",overwrite_layer = T)

#use the cartography package
library(cartography)
choroLayer(spdf=shp,nclass=5,var="dimfrac",method="fisher-jenks")
}}