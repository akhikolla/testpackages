\name{covQual-class}
\Rdversion{1.1}
\docType{class}
\alias{covQual-class}
\alias{checkX,covQual,matrix-method}
\alias{checkX,covQual,data.frame-method}
\alias{coef,covQual-method}
\alias{coef<-,covQual,numeric-method}
\alias{coefLower<-,covQual-method}
\alias{coefLower,covQual-method}
\alias{coefUpper<-,covQual-method}
\alias{coefUpper,covQual-method}
\alias{coerce,covQual,function-method}
\alias{covMat,covQual-method}
%% \alias{covMat,covMan-method}  % see covMat.covMan.Rd
\alias{npar,covQual-method}
\alias{show,covQual-method}
\alias{scores,covQual-method}
\alias{simulate,covQual-method}
\alias{varVec,covQual-method}
\title{
  Class \code{"covQual"}
}
\description{

  Covariance kernel for qualitative inputs.

}
\section{Objects from the Class}{
  
  Objects can be created by calls of the form \code{new("covQual", ...)}.
  
  
}
\section{Slots}{
  
  \describe{
    \item{\code{covLevels}:}{

      Object of class \code{"function"}. This function has
      arguments \code{'par'} and optional arguments \code{lowerSQRT} and
      \code{compGrad}. It returns the covariance matrix for an input 
      corresponding to all the levels. 
      
    }
    \item{\code{covLevMat}:}{

      Object of class \code{"matrix"}. This is the result returned by the 
      function \code{covLevels} (former slot) with \code{lowerSQRT = FALSE}
      and \code{gradient = FALSE}.

    }
    \item{\code{hasGrad}:}{

      Object of class \code{"logical"}.  When \code{TRUE}, the
      covariance matrix returned by the function in slot
      \code{covLevels} must compute the gradients. The returned 
      covariance matrix must have a \code{"gradient"} attribute;
      this must be an array with dimension \code{c(m, m, np)} where
      \code{m} stands for the number of levels and \eqn{np} is the
      number of parameters.

    }
    \item{\code{acceptLowerSQRT}:}{

      Object of class \code{"logical"}. When \code{TRUE}, the function
      in slot \code{covLevels} must have a formal \code{lowerSQRT}
      which can receive a logical value. When the value is \code{TRUE}
      the Cholesky (lower) root of the covariance is returned instead of
      the covariance.

    }
    \item{\code{label}:}{

      Object of class \code{"character"}. A description of the kernel
      which will remained attached with it.

    }
    \item{\code{d}:}{

      Object of class \code{"integer"}. The dimension or number of
      (qualitative) inputs of the kernel.

    }
    \item{\code{inputNames}:}{

      Object of class \code{"character"}. The names of the (qualitative)
      inputs. These will be matched against the columns of a data frame
      when the kernel will be evaluated.

    }
    \item{\code{nlevels}:}{

      Object of class \code{"integer"}. A vector with length
      \code{d} giving the number of levels for each of the \code{d}
      inputs.

    }
    \item{\code{levels}:}{

      Object of class \code{"list"}. A list of length \code{d}
      containing the \code{d} character vectors of levels for
      the \code{d} (qualitative) inputs.

    }
    \item{\code{parLower}:}{

      Object of class \code{"numeric"}. Vector of \code{parN} lower
      values for the parameters of the structure. The value
      \code{-Inf} can be used when needed.

    }

    \item{\code{parUpper}:}{

      Object of class \code{"numeric"}. Vector of \code{parN} upper
      values for the parameters of the structure. The value
      \code{Inf} can be used when needed.
      
    }
    \item{\code{par}:}{

      Object of class \code{"numeric"}. Vector of \code{parN}
      current values for the structure.

    }
    \item{\code{parN}:}{

      Object of class \code{"integer"}. Number of parameters for the
      structure, as returned by the \code{npar} method.

    }
    \item{\code{kernParNames}:}{

      Object of class \code{"character"}. Vector of length \code{parN}
      giving the names of the parameters. E.g. \code{"range"},
      \code{"var"}, \code{"sigma2"} are popular names.

    }
    \item{\code{ordered}:}{

      Vector of class \code{"logical"} indicating whether the factors
      are ordered or not.
    
    }
    \item{\code{intAsChar}:}{

      Object of class \code{"logical"} indicating how to cope with an
      integer input. When \code{intAsChar} is \code{TRUE} the input is
      coerced into a character; the values taken by this character
      vector should then match the levels in the \code{covQual} object
      as given by \code{levels(object)[[1]]}. If instead
      \code{intAsChar} is \code{FALSE}, the integer values are assumed
      to correspond to the levels of the \code{covQual} object in the
      same order.
      
    }
  }
}
\section{Methods}{
  
  \describe{
    \item{checkX}{

      \code{signature(object = "covQual", X = "data.frame")}: check that
      the inputs exist with suitable column names and suitable \emph{factor}
      content. The levels should match the prescribed levels. Returns a
      matrix with the input columns in the order prescribed by
      \code{object}.
      
      \code{signature(object = "covQual", X = "matrix")}: check that the
      inputs exist with suitable column names and suitable \emph{numeric}
      content for coercion into a factor with the prescribed levels.
      Returns a data frame with the input columns in the order
      prescribed by \code{object}.

    }
    \item{coef<-}{
      
      \code{signature(object = "covQual")}: replace the whole vector of
      coefficients, as required during ML estimation.
      
    }
    \item{coefLower<-}{

      \code{signature(object = "covQual")}: replacement method for lower
      bounds on covQual coefficients.

    }
    \item{coefLower}{
      
      \code{signature(object = "covQual")}: extracts the numeric values of
      the lower bounds.
      
    }
    \item{coef}{
      
      \code{signature(object = "covQual")}: extracts the numeric values
      of the covariance parameters.
     

    }
    \item{coefUpper<-}{

      \code{signature(object = "covQual")}: replacement method for upper
      bounds on \code{covQual} coefficients.

    }
    \item{coefUpper}{

      \code{signature(object = "covQual")}: ...

    }
    \item{covMat}{

      \code{signature(object = "covQual")}: build the covariance matrix
      or the cross covariance matrix between two sets of locations for a
      \code{covQual} object.

    }
    \item{npar}{

      \code{signature(object = "covQual")}: returns the number of
      parameters.

    }
    
    \item{plot}{\code{signature(x = "covQual")}: see \code{\link{plot,covQual-method}}. }

    \item{scores}{

      \code{signature(object = "covQual")}: return the vector of
      scores, i.e. the derivative of the log-likelihood w.r.t. the
      parameter vector at the current parameter values.

    }
    \item{simulate}{

      \code{signature(object = "covQual")}: simulate \code{nsim} paths
      from a Gaussian Process having the covariance structure. The paths
      are indexed by the finite set of levels of factor inputs, and they
      are returned as columns of a matrix.

    }
    \item{varVec}{

      \code{signature(object = "covQual")}: build the variance vector
      corresponding to a set locations for a \code{covQual} object.

    }
  }

}

%% \references{}
%%\author{}
\note{

  This class is to be regarded as experimental. The slot names or list
  may be changed in the future. The methods \code{npar},
  \code{inputNames} or \code{`inputNames<-`} should provide a more
  robust access to some slot values.
  
}
\seealso{
  
  See \code{\linkS4class{covMan}} for a comparable structure dedicated
  to kernels with continuous inputs.

}
\examples{

showClass("covQual")

}
\keyword{classes}
