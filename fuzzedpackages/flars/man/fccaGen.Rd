\name{fccaGen}
\alias{fccaGen}
\title{
Functional canonical correlation analysis between a scalar variable and a list of mixed scalar and functional variables.
}
\description{
This function carries out the canonical correlation analysis between a scalar variable and a list of mixed scalar and functional variables. There are four choices of the returned values and three representation methods of the functional variables. 
}
\usage{
fccaGen(xL,yVec,type=c('dir','cor','a','all'),method=c('basis',
      'gq','raw'),GCV=TRUE,control=list())
}
\arguments{
  \item{xL}{
The mixed scalar and functional variables. If there is only one functional variable, \code{xL} can be a matrix. If there is only scalar variables, \code{xL} can be a vector or a matrix. If there are more than one functional variables, or there are mixed functional and scalar variables, \code{xL} should be a list. If \code{xL} is a list, each item of the list should correspond to one variable. 
}
  \item{yVec}{
The scalar variable. It should be a matrix. 
}
  \item{type}{
The choice of outcomes. See details for more information.
}
  \item{method}{
The representative methods for the functional coefficients. The method could be one of the 'basis', 'gq' and 'raw' for basis function expression, Gaussian quadrature and representative data points, respectively. 
}
  \item{GCV}{
Use generalized cross validation (GCV) or not to choose the tuning parameter. Logic argument. Currently the only choice is to use GCV.
}
  \item{control}{
List of elements that controls the details of the functional coefficients. See details for more information.
}
}

\details{
There are four choices of \code{type} in the function. 'dir' means that the function only returns the direction coefficients like the one in the traditional Canonical correlation analysis. 'cor' means that the function only returns the correlation coefficients. 'a' means that the function only returns the normalized direction coefficients. With this normalization, the direction coefficients are equivalent to the coefficients from a linear regression with response variable \code{yVec} and covariates \code{xL}. 'all' means that the function returns all three outcomes mentioned above. 

The argument \code{control} is a list. It changes when different representative methods are used for the functional coefficients. If (type=='basis'), the list contains the following items: 
\itemize{
\item{nbasis: }{Number of B-spline basis functions. Default value is 18.}

\item{norder: }{Order of the basis functions. Default value is 6.}

\item{pen1: }{The candidate values of the smoothing parameter. Default values are \item{10^(seq((-20),5,len=41))}}

\item{pen2: }{The candidate values of the ridge tuning parameter. Default value is \item{0.01}}

\item{t: }{IMPORTANT! The time points correspond to the discrete data points of the functional variables. Default to be \code{seq(0,1,len=max(sapply(xL,ncol),na.rm = T))}. Do NOT change the starting and ending point of the sequence.}
}


If (type=='gq'), the list contains the following items: 
\itemize{
\item{nP: }{Number of Gaussian quadrature points. Default value is 18.} 

\item{pen1: }{The candidate values of the smoothing parameter. Default values are \item{10^(seq((-20),5,len=21))}}

\item{pen2: }{The candidate values of the ridge tuning parameter. Default value is \item{0.01}}

\item{t: }{IMPORTANT! The time points correspond to the discrete data points of the functional variables. Default to be \code{seq(-1,1,len=max(sapply(xL,ncol),na.rm = T))}. Do NOT change the starting and ending point of the sequence.}
}


If (type=='raw'), the list contains the following items: 

\itemize{
\item{pen1: }{The candidate values of the smoothing parameter. Default values are \item{10^(seq((-20),5,len=21))}}

\item{pen2: }{The candidate values of the ridge tuning parameter. Default value is \item{0.01}}

\item{t: }{IMPORTANT! The time points correspond to the discrete data points of the functional variables. Default to be \code{seq(0,1,len=max(sapply(xL,ncol),na.rm = T))}. Do NOT change the starting and ending point of the sequence.}
}

The function is designed to be able to handle the situation when different functional variables have different number of discrete data points and the discrete data points could be non-evenly spaced. This would require a list of \code{t} to input in the argument. However, this is not fully tested at the moment. For convenient, especially when we have a large number of functional variables, a universal setting of \code{t} is recommended.
}

\value{
\item{corr}{Correlation coefficient. It is returned when \code{type='corr'} or \code{type='all'}}
\item{a}{Normalized direction coefficients. It is returned when \code{type='a'} or \code{type='all'}}
\item{dir}{Direction coefficients. It is returned when \code{type='dir'}}
\item{K}{Penalized covariance matrix. It is returned when \code{type='all'}.}
\item{gq}{Gaussian quadrature weights. It is returned when \code{type='all'}.}
\item{phiL}{Known part of the functional coefficients. E.g, basis functions. It is returned when \code{type='all'}.}
\item{S}{Hat matrix. It is returned when \code{type='all'}.}
\item{lam1}{The selected smoothing parameter. It is returned when \code{type='all'}.}
\item{lam2}{The selected ridge parameter. It is returned when \code{type='all'}.}
\item{GCV_mat}{The GCV value. It is returned when \code{type='all'}.}
\item{TraceHat}{Trace of the hat matrix. It is returned when \code{type='all'}.}
}

\examples{
library(flars)
## Generate some data.
dataL=data_generation(seed = 1,uncorr = TRUE,nVar = 8,nsamples = 120,
      var_type = 'm',cor_type = 1)

## If there is only one functional variable
# out1=fccaGen(dataL$x[1], dataL$y, type='all', method='basis')

## If there are only a few scalar variables
# x=matrix(rnorm(3*length(dataL$y)),ncol=3)
# out2=fccaGen(x, dataL$y, type='all', method='basis')

## If there are mixed scalar and functional variables
# out3=fccaGen(dataL$x, dataL$y, type='all', method='basis')
}

