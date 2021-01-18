\name{inputNames}
\alias{inputNames}
\alias{inputNames,covAll-method}
\alias{inputNames<-}
\alias{inputNames<-,covAll-method}
\alias{inputNames<-,covAll,character-method}
\title{

  Generic Function: Names of the Inputs of a Covariance Kernel

}
\description{
  
  Generic function returning or setting the names of the inputs attached
  with a covariance kernel.
  
}
\usage{

inputNames(object, ...)

\S4method{inputNames}{covAll}(object, ...) <- value

}
\arguments{
  \item{object}{

    A covariance kernel object.

  }
  \item{value}{

    A suitable character vector.

  }
  \item{\dots}{

    Other arguments for methods.

  }
}

\value{

  A character vector with \emph{distinct} input names that are used
  e.g. in prediction.

}

%%\author{}

\note{

  The input names are usually checked to control that they match the
  colnames of a spatial design matrix. They play an important role since
  in general the inputs are found among other columns of a data frame,
  and their order is not fixed. For instance in a data frame used as
  \code{newdata} formal in the \code{predict} method, the inputs are
  generally found at positions which differ from those in the data frame
  used at the creation of the object.
  
}

\seealso{
  \code{\link{checkX}}
}
