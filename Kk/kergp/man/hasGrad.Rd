\name{hasGrad}
\alias{hasGrad}
\alias{hasGrad,covAll-method}
\title{

  Generic Function: Extract slot hasGrad of a Covariance Kernel

}
\description{
  
  Generic function returning the slot hasGrad of a Covariance Kernel.
  
}
\usage{

hasGrad(object, ...)

\S4method{hasGrad}{covAll}(object, ...)

}
\arguments{
  \item{object}{

    A covariance kernel object.

  }
  \item{\dots}{

    Other arguments for methods.

  }
}

\value{

  A logical indicating whether the gradient is supplied in \code{object} (as indicated in the slot 'hasGrad').

}

