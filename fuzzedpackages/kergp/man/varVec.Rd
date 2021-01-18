\name{varVec}
\alias{varVec}

\title{

  Generic Function: Variance of Gaussian Process at Specific Locations
  
}
\description{

  Generic function returning a variance vector 
  
}
\usage{

varVec(object, X, ...)

}

\arguments{
  \item{object}{
    
    Covariance kernel object.

  }
  \item{X}{

    A matrix with \eqn{d} columns, where \eqn{d} is the number of inputs
    of the covariance kernel. The \eqn{n} rows define a set of sites or
    locations.

  }
  \item{\dots}{
    
    Other arguments for methods.

  }
}

\value{
  
  A numeric vector with length \code{nrow(X)} containing the variances
  \eqn{K(\mathbf{x}, \mathbf{x})}{K(x, x)} for all the sites
  \eqn{\mathbf{x}}{x}.
  
}

%% \seealso{
%%  
%%  The \code{\link[DiceKriging]{covMat1Mat2}} function.
%%
%% }


