\name{hfactor}
\alias{hfactor}
\title{Constrained Factor Model
}
\description{Performs factor model analysis with a given 
constrained matrix
}
\usage{
hfactor(X, H, r)
}
\arguments{
  \item{X}{A T-by-k data matrix of an observed k-dimensional time series
}
  \item{H}{The constrained matrix with each column representing a constraint
}
  \item{r}{The number of common factor
}
}
\value{Results of the traditional PCA and constrained factor models are given
}
\references{Tsay (2014, Chapter 6). Tsai and Tsay (2010, JASA)
}
\author{Ruey S. Tsay
}
\examples{
data("mts-examples",package="MTS")
rtn=log(tenstocks[,2:11]+1) # compute log returns
h1=c(1,1,1,1,rep(0,6)) # specify the constraints
h2=c(0,0,0,0,1,1,1,0,0,0)
h3=c(rep(0,7),1,1,1)
H=cbind(h1,h2,h3)
m1=hfactor(rtn,H,3)
}




