\name{alluvial.plot}
\alias{alluvial.plot}
\title{
Alluvial  plot showing  the dynamics  of the  group memberships  associated to a dynamic stochastic block model.
}
\description{
Alluvial  plot showing  the dynamics  of the  group memberships  associated to a dynamic stochastic block model.
}
\usage{
alluvial.plot(dynsbm, timestep.abbrev="T", minimal.flow=1, only.present=FALSE)
}
\arguments{
\item{dynsbm}{An object of class \code{dynsbm} retrieved with the function \code{select.dynsbm}.}
\item{timestep.abbrev}{Time step abbreviation (D fors days, W for weeks,...).}
\item{minimal.flow}{Minimal flow to be displayed, i.e. minimal number of nodes to consider.}
\item{only.present}{Enable/disable the appearence of the group O which gathers absent nodes.}
}
\details{
Alluvial flows between groups. 
Between any two time plots (on the x-axis), lines are drawn to represent flows between groups.
Each flow corresponds to the membership overlapp between one a group to another group at two consecutive time steps  (represented on the y-axis).  
Each group q at each time step t is called Xt_q where X is time step abbreviation set with \code{timestep.abbrev}. 
The thickness of each line is  proportional  to the  corresponding node counts.
}
\references{
Catherine Matias and Vincent Miele, \emph{Statistical clustering of temporal networks through a dynamic stochastic block model}, Journal of the Royal Statistical Society: Series B (2017)
http://dx.doi.org/10.1111/rssb.12200
http://arxiv.org/abs/1506.07464

Vincent Miele and Catherine Matias, \emph{Revealing the hidden structure
  of dynamic ecological networks},  Royal Society Open Science (2017)
http://dx.doi.org/10.1098/rsos.170251
https://arxiv.org/abs/1701.01355
}
\author{
Authors: Catherine Matias, Vincent Miele
  
Maintainer: Vincent Miele <vincent.miele@univ-lyon1.fr>
}
\examples{
####################
## 1 - binary case
data(simdataT5Q4N40binary)

## estimation for Q=1..5 groups
list.dynsbm <- select.dynsbm(simdataT5Q4N40binary, 
				Qmin=1, Qmax=5, edge.type="binary", nstart=1)
				
\dontrun{
## better to use nstart>1 starting points
## but estimation can take 1-2 minutes
list.dynsbm <- select.dynsbm(simdataT5Q4N40binary, 
				Qmin=1, Qmax=5, edge.type="binary", nstart=25)
}

## selection of Q=4
dynsbm <- list.dynsbm[[4]]

## plotting switches between groups
alluvial.plot(dynsbm)

####################
## 2 - continuous case
data(simdataT5Q4N40continuous)

## estimation for Q=1..5 groups
list.dynsbm <- select.dynsbm(simdataT5Q4N40continuous, 
				Qmin=1, Qmax=5, edge.type="continuous", nstart=1)
						
\dontrun{
## better to use nstart>1 starting points
## but estimation can take 1-2 minutes
list.dynsbm <- select.dynsbm(simdataT5Q4N40continuous, 
				Qmin=1, Qmax=5, edge.type="continuous", nstart=25)
}

## selection of Q=4
dynsbm <- list.dynsbm[[4]]

## plotting switches between groups
alluvial.plot(dynsbm)

####################
## 3 - discrete case
data(simdataT5Q4N40discrete)

## estimation for Q=1..5 groups
list.dynsbm <- select.dynsbm(simdataT5Q4N40discrete, 
				Qmin=1, Qmax=5, edge.type="discrete", K=4, nstart=1)
									
\dontrun{
## better to use nstart>1 starting points
## but estimation can take 1-2 minutes
list.dynsbm <- select.dynsbm(simdataT5Q4N40discrete, 
				Qmin=1, Qmax=5, edge.type="discrete", K=4, nstart=25)
}

## selection of Q=4
dynsbm <- list.dynsbm[[4]]

## plotting switches between groups
alluvial.plot(dynsbm)
}