\name{covRadial-class}
\Rdversion{1.1}
\docType{class}
\alias{covRadial-class}
\alias{coef<-,covRadial,numeric-method}
\alias{coefLower<-,covRadial-method}
\alias{coefLower,covRadial-method}
\alias{coef,covRadial-method}
\alias{coefUpper<-,covRadial-method}
\alias{coefUpper,covRadial-method}
\alias{covMat,covRadial-method}
\alias{npar,covRadial-method}
\alias{scores,covRadial-method}
\alias{show,covRadial-method}
\alias{varVec,covRadial-method}

\title{Class \code{"covRadial"}}

\description{

  Class of radial covariance kernels.
  
}
\section{Objects from the Class}{
  
  Objects can be created by calls of the form \code{covRadial(...)}  of
  \code{new("covRadial", ...)}.

  
}
\section{Slots}{
  \describe{
    \item{\code{k1Fun1}:}{

      Object of class \code{"function"} A function of a scalar numeric
      variable. Not that using a one-dimensional kernel here \emph{does
      not} warrant that a positive semi-definite kernel results for any
      dimension \eqn{d}.

    }
    \item{\code{hasGrad}:}{

      Object of class \code{"logical"}. Tells if the value returned by
      the function \code{kern1Fun} has an attribute named \code{"der"}
      giving the derivative(s). 

    }
    \item{\code{cov}:}{

      Object of class \code{"integer"}. The value \code{0L} corresponds
      to a correlation kernel while \code{1L} is for a covariance
      kernel.

    }
    \item{\code{iso}:}{

      Object of class \code{"integer"}. The value \code{1L} corresponds
      to an isotropic covariance, with all the inputs sharing the same
      range value.

    }
    \item{\code{label}:}{

      Object of class \code{"character"}. Short description of the
      object.

    }
    \item{\code{d}:}{

      Object of class \code{"integer"}. Dimension, i.e. number of
      inputs.

    }
    \item{\code{inputNames}:}{

      Object of class \code{"optCharacter"}. Names of the inputs.

    }
    \item{\code{parLower}:}{

      Object of class \code{"numeric"}. Numeric values for the lower
      bounds on the parameters. Can be \code{-Inf}.

    }
    \item{\code{parUpper}:}{

      Object of class \code{"numeric"}. Numeric values for the upper
      bounds on the parameters. Can be \code{Inf}.

    }
    \item{\code{par}:}{

      Object of class \code{"numeric"}. Numeric values for the
      parameters. Can be \code{NA}.

    }
    \item{\code{parN1}:}{

      Object of class \code{"integer"}. Number of parameters of the
      function \code{kern1Fun} (such as a shape).

    } \item{\code{parN}:}{

      Object of class \code{"integer"}? Number of parameters for the
      object. The include: \emph{direct} parameters in the function
      \code{kern1Fun}, ranges, and variance.

    }
    \item{\code{kern1ParNames}:}{

      Object of class \code{"character"}. Names of the \emph{direct}
      parameters.

    }
    \item{\code{kernParNames}:}{

      Object of class \code{"character"}?. Names of the parameters.

    }
  }
}
\section{Extends}{
  
  Class \code{"\linkS4class{covAll}"}, directly.

}
\section{Methods}{
  \describe{
    \item{coef<-}{
      
      \code{signature(object = "covRadial", value = "numeric")}: Set the
      vector of values for the parameters.
      
    }
    \item{coefLower<-}{
      
      \code{signature(object = "covRadial")}: Set the vector of lower
      bounds on the parameters. 

    }
    \item{coefLower}{

      \code{signature(object = "covRadial")}: Get the vector of
      lower bounds on the parameters.

    }
    \item{coef}{
      
      \code{signature(object = "covRadial")}: Get the vector of values
      for the parameters.

    }
    \item{coefUpper<-}{

      \code{signature(object = "covRadial")}: Set the vector of upper
      bounds on the parameters. 

    }
    \item{coefUpper}{

      \code{signature(object = "covRadial")}: Get the vector of upper
      bounds on the parameters. 

    }
    \item{covMat}{

      \code{signature(object = "covRadial")}: Compute the covariance
      matrix for given sites.

    }
    \item{npar}{

      \code{signature(object = "covRadial")}: Get the number of
      parameters.

    }
    \item{scores}{

      \code{signature(object = "covRadial")}: Compute the scores
      i.e. the derivatives w.r.t. the parameters of the contribution of
      the covariance in the log-likelihood of a \code{gp}.

    }
    \item{show}{

      \code{signature(object = "covRadial")}: Print or show the object.

    }
    \item{varVec}{

      \code{signature(object = "covRadial")}: Compute the variance
      vector for given sites.

    }  } }
%% \references{}
%% \note{}
\seealso{

  The creator function \code{\link{covRadial}}, where some details
  are given on the form of kernel.
  \code{\link{covMan}} and \code{\linkS4class{covMan}} for a comparable
  but more general class.

}
\examples{
showClass("covRadial")
}
\keyword{classes}
