\name{checkX}
\alias{checkX}

%%\alias{checkX}

\title{
  
  Generic function: Check the Compatibility of a Design Matrix with a
  Given Covariance Object

}
\description{
  
  Generic function to check the compatibility of a design matrix with a
  covariance object.

}

\usage{

   checkX(object, X, ...)

}

\arguments{
  \item{object}{

    A covariance kernel object.

  }
  \item{X}{

    A design matrix.

  }
  \item{\dots}{

    Other arguments for methods.

  }
  
}

\value{

  A matrix with columns taken from \code{X} and with column names
  identical to \code{inputNames(object)}.

}

%% \author{Yves Deville}

\seealso{The \code{\link{inputNames}} method.}
