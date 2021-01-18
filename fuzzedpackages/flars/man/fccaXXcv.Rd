\name{fccaXXcv}
\alias{fccaXXcv}
\title{
This function finds the best smoothing parameter for the canonical correlation analysis for both groups of variables by using leave-one-out (sample) cross validation. The criterion here is to maximise the first canonical correlation. 
}
\description{
This function carries out the canonical correlation analysis between a scalar variable and a list of mixed scalar and functional variables. There are four choices of the returned values and three representation methods of the functional variables. 
}
\usage{
fccaXXcv(xL1,xL2,method=c('basis','gq','raw'),centre = TRUE,tol=1e-7,
      Control1=list(),Control2=list(),alpha=10^seq(-6,1,len=10))
}
\arguments{
  \item{xL1}{
The mixed scalar and functional variables. For any number and any type of variables, \code{xL1} should be a list. Each item of the list should correspond to one variable. 
}
  \item{xL2}{
Same as \code{xL1}.
}
  \item{method}{
The representative methods for the functional coefficients. The method could be one of the 'basis', 'gq' and 'raw' for basis function expression, Gaussian quadrature and representative data points, respectively. 
}
  \item{centre}{
Logic argument. Default is \code{TRUE}, which means the variables do need to be centred.
}
  \item{tol}{
The threshold to decide whether the correlation is to small to be non-zero.   
}
  \item{Control1}{
List of elements that controls the details of the functional coefficients for \code{xL1}. See details for more information. See the argument \code{control} in function \code{\link{fccaGen}} for details.
}
  \item{Control2}{
Similar to \code{Control1}. 
}
  \item{alpha}{
Candidate tuning parameters for the smoothness of the functional coefficients.
}
}

\details{
Note that the smoothing parameters for both groups of variables are assumed to be the same. This is due to high computational cost of cross validation. See the example in \code{\link{fccaXX}}.
}

\value{
\item{cor}{A vector of the first canonical correlation. Each element of the vector is corresponding to one of the candidate tuning parameters.}
\item{alpha}{The corresponding tuning parameters.}
}


