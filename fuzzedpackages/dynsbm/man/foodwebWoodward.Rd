\name{foodwebWoodward}
\docType{data}
\alias{foodwebWoodward}
\title{Broadstone Stream seasonal food webs}
\description{
This dataset concerns the aquatic macro-invertebrate community of
Broadstone Stream in south-east England. Six seasonal connectance food  webs 
were  recorded, one every  two months  from May 1996  to April  1997. 
We  restricted here  to simple  presence/absence information on  species (nodes) 
and  binary feeding links  (edges).  This  dataset  forms  a dynamic  trophic  network  
with T=6 snapshots (May, August, October,  December 1996, February, April 1997).
Five species were not sampled each month. Each netowrk is directed.
}
\usage{data(foodwebWoodward)}
\format{An array of size 6x26x26.}
\source{Table 2 of Woordward et al}
\references{
Woodward, G., Speirs, D. C., Hildrew, A. G., \emph{Quantification and Resolution of a Complex, Size-Structured Food Web. In Food Webs: From Connectivity to Energetics}, Vol. 36 of Advances in Ecological Research, pp. 85-135 Academic Press (2005) 
http://dx.doi.org/10.1016/S0065-2504(05)36002-8

Miele, V and Matias, C, \emph{Revealing the hidden structure of dynamic ecological networks}, Royal Society Open Science (2017)
}
\examples{
\dontrun{
data(foodwebWoodward)

## mandatory to use many nstart>>1 starting points
## but estimation can be long;
## better to use nb.cores>1 cores
list.dynsbm <- select.dynsbm(Y=foodwebWoodward$Y, 
				present=foodwebWoodward$present, 
				Qmin=1, Qmax=6, edge.type="binary",
				directed=TRUE, self.loop=TRUE, 
				nstart=200, nb.cores=4) 
				

## selection of Q=4 with the ICL method
dynsbm <- list.dynsbm[[4]]

## plotting intra/inter connectivity patterns
connectivity.plot(dynsbm, foodwebWoodward$Y)

## plotting switches between groups
alluvial.plot(dynsbm, timestep.abbrev="D", only.present=F)
}
}
