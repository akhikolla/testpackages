\name{checkGrad}
\alias{checkGrad}
\title{

  Check the Gradient Provided in a \code{covMan} Object

}
\usage{
checkGrad(object, sym = TRUE,
          x1 = NULL, n1 = 10,
          x2 = NULL, n2 = NULL,
          XLower = NULL, XUpper = NULL,
          plot = TRUE)
}
\arguments{
  \item{object}{

    A \code{covMan} object.

  }
  \item{sym}{

    Logical. If \code{TRUE}, the check is done assuming that \code{x2}
    is identical to \code{x1}, so the provided values for \code{x2} and
    \code{n2} (if any) will be ignored.
    
  }
  \item{x1}{

    Matrix to be used as the first argument of the kernel.

  }
  \item{n1}{

    Number of rows for the matrix \code{x1}. Used only when \code{x1} is
    not provided.

  }
  \item{x2}{

    Matrix to be used as the second argument of the kernel.

  }
  \item{n2}{

    Number of rows for the matrix \code{x2}. Used only when \code{x2} is
    not provided.

  }
  \item{XLower}{

    Vector of lower bounds to draw \code{x1} and \code{x2} when needed.

  }
  \item{XUpper}{

    Vector of upper bounds to draw \code{x1} and \code{x2} when needed.

  }
  \item{plot}{

    \code{Logical}. If \code{TRUE}, a plot is shown comparing the two
    arrays of gradients.
    
  }
}
\value{
  
  A list of results related to the Jacobians

  \itemize{
    
    \item{\code{test}}{ Max of the absolute difference between
      the gradient obtained by numeric differentiation and the
      gradient provided by the kernel object. }
    
    \item{\code{Jnum}, \code{J}}{ Jacobians (arrays) computed
      with \code{numDeriv::jacobian} and provided by the kernel
      object. }
    
    \item{\code{x1}, \code{x2}, \code{K}}{ The matrices used
      for the check, and the matrix of kernel values with
      dimension \code{c(n1, n2)}. The element \code{x2} can be
      \code{NULL} if the determination of the matrix \code{x2}
      was not necessary. }
    
  }
}
\description{

  Check the gradient provided in a \code{covMan} object.

}
\details{

  Each of the two matrices \code{x1} and \code{x2} with \code{n1} and
  \code{n2} rows can be given or instead be drawn at random. The matrix
  of kernel values with dimension \code{c(n1, n2)} is computed, together
  with its gradient with dimension \code{c(n1, n2, npar)} where
  \code{npar} is the number of parameters of the kernel. A numerical
  differentiation w.r.t. the kernel parameters is performed for the
  kernel value at \code{x1} and \code{x2}, and the result is compared to
  that provided by the kernel function (the function described in the
  slot named \code{"kernel"} of \code{object}). Note that the value of
  the parameter vector is the value provided by \code{coef(object)} and
  it can be changed by using the replacement method \code{`coef<-`} if
  needed.

}
\note{
  
  As a rule of thumb, a gradient coded without error gives a value of
  \code{test} less than \code{1e-4}, and usually the value is much
  smaller than that.

}
\section{Caution}{
  
  For now the function only works when \code{object} has class
  \code{"covMan"}.

}
\author{
  
  Yves Deville

}
