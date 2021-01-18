\name{covMat}
\alias{covMat}

\title{

  Generic Function: Covariance or Cross-Covariance
  Matrix Between two Sets of Locations

}
\description{

  Generic function returning a covariance or a cross-covariance matrix
  between two sets of locations.
  
}
\usage{

covMat(object, X, Xnew, ...)

}

\arguments{
  \item{object}{
    
    Covariance kernel object.

  }
  \item{X}{

    A matrix with \code{d} columns, where \code{d} is the number of inputs
    of the covariance kernel. The \eqn{n_1} rows define a first set of sites or
    locations, typically used for learning.

  }
  \item{Xnew}{

    A matrix with \code{d} columns, where \code{d} is the number of inputs
    of the covariance kernel. The \eqn{n_2} rows define a second set of
    sites or locations, typically used for testing or prediction.
    If \code{Xnew = NULL} the same locations are used: \code{Xnew = X}.

  }
  \item{\dots}{
    
    Other arguments for methods.

  }
}

\value{
  
  A rectangular matrix with \code{nrow(X)} rows and \code{nrow(Xnew)}
  columns containing the covariances \eqn{K(\mathbf{x}_1,
    \mathbf{x}_2)}{K(x1, x2)} for all the couples of sites
  \eqn{\mathbf{x}_1}{x1} and \eqn{\mathbf{x}_2}{x2}.

}

%% \seealso{
%%  
%%  The \code{\link[DiceKriging]{covMat1Mat2}} function.
%%
%% }


