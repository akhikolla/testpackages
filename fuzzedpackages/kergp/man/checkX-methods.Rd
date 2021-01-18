\name{checkX-methods}
\docType{methods}
\alias{checkX,covAll-method}
%%\alias{checkX}

\title{
  
  Check the Compatibility of a Design with a Given
  Covariance Object

}
\description{
  
  Check the compatibility of a design matrix with a covariance object.

}

\usage{

   \S4method{checkX}{covAll}(object, X, strict = FALSE, ...)

}

\arguments{
  \item{object}{

    A covariance kernel object.

  }
  \item{X}{

    A design matrix or data frame.

  }
  \item{strict}{

    Logical. If \code{TRUE}, the character vectors \code{colnames(X)}
    and \code{inputNames(object)} must be the same sets, and hence have
    the same length. If \code{FALSE} the vector
    \code{inputNames(object)} must be a subset of \code{colnames(X)}
    which then can have unused columns.

  }
  \item{\dots}{

    Not used yet.
    
  }
}

\details{

  The matrix \code{X} must have the number of columns expected from the
  covariance kernel object description, and it must have named columns
  conforming to the kernel input names as returned by the
  \code{\link{inputNames}} method. If the two sets of names are identical but
  the names are in a different order, the columns are permuted in order
  to be in the same order as the input names. If the names sets differ,
  an error occurs.

}

\value{

  A matrix with columns names identical to the input names attached with
  the kernel object, i.e. \code{inputNames(object)}. The columns are
  copies of those found under the same names in \code{X}, but are put in
  the order of \code{inputNames(object)}.  When an input name does not
  exist in \code{colnames(X)} an error occurs.

}

%% \author{Yves Deville}

\seealso{The \code{\link{inputNames}} method.}
