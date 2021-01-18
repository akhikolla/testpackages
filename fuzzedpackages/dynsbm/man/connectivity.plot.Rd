\name{connectivity.plot}
\alias{connectivity.plot}
\title{
Plot the connectivity characteristics between groups associated to a dynamic stochastic block model.
}
\description{  
Plot the connectivity characteristics between groups associated to a dynamic stochastic block model.
}
\usage{
connectivity.plot(dynsbm, Y)
}
\arguments{
\item{dynsbm}{An object of class \code{dynsbm} retrieved with the function \code{select.dynsbm}.}
\item{Y}{An object of class \code{array} of dimension (T x N x N) containing T adjacency matrices of size (N x N), where N is the number of nodes in the network and T is the number of time points.}
}
\details{
Interaction presence and intensity  between  nodes in  any of  the groups to  the others are represented  in a  QxQ matrix. 
The cell in line q/column l deals with the connectivity between groups q/l. 
Each cell contains a curve with T time points on the x-axis corresponding to the T proportions of present edges over all the possible edges, where Q is the number of groups and T is the number of time points, and  
If \code{dynsbm} was estimated with \code{edge.type=="binary"}, the area below the curve is filled in light blue.
If \code{dynsbm} was estimated with \code{edge.type=="discrete"}, the area below the curve is divided into K areas corresponding to the proportion of edges  with value 1 to K (the darker blue, the greater edge intensity).
If \code{dynsbm} was estimated with \code{edge.type=="continuous"}, the area below the curve is filled with a colored gradient representing the mean edge intensity (the darker blue, the greater).
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

## plotting intra/inter connectivity patterns
connectivity.plot(dynsbm, simdataT5Q4N40binary)

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

## plotting intra/inter connectivity patterns
connectivity.plot(dynsbm, simdataT5Q4N40continuous)

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

## plotting intra/inter connectivity patterns
connectivity.plot(dynsbm, simdataT5Q4N40discrete)
}