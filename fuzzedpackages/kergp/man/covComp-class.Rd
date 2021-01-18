\name{covComp-class}
\Rdversion{1.1}
\docType{class}
\alias{covComp-class}
\alias{checkX,covComp,data.frame-method}
\alias{coef,covComp-method}
\alias{coef<-,covComp,numeric-method}
%% \alias{coefLower<-,covComp-method}
\alias{coefLower,covComp-method}
%% \alias{coefUpper<-,covComp-method}
\alias{coefUpper,covComp-method}
%% \alias{coerce,covQual,function-method}
\alias{covMat,covComp-method}
\alias{inputNames<-,covComp,character-method}
\alias{show,covComp-method}
\alias{scores,covComp-method}
\alias{varVec,covComp-method}
\alias{as.list,covComp-method}

\title{

  Class \code{"covComp"}

}
\description{
  
  Class \code{"covComp"} representing a composite kernel combining
    several kernels objects inheriting from the class \code{"covAll"}.

}
\section{Objects from the Class}{

  Objects can be created by calls of the form \code{new("covComp",
    ...)} or by using \code{covComp}.
  
  
}
\section{Slots}{
  
  \describe{
    \item{\code{def}:}{

      Object of class \code{"expression"} defining the
      This is a parsed and cleaned version of the value
      of the \code{formula} formal in \code{\link{covComp}}.

    }
    \item{\code{covAlls}:}{

      Object of class \code{"list"} containing the kernel
      objects used by the formula. The coefficients of these
      kernels can be changed.

    }
    \item{\code{hasGrad}:}{

      Object of class \code{"logical"}: can we differentiate
      the kernel w.r.t. all its parameters?

    }
    \item{\code{label}:}{

      Object of class \code{"character"} A label attached to the kernel
      to describe it.

    }
    \item{\code{d}:}{

      Object of class \code{"integer"}: dimension (or number of inputs).

    }
    \item{\code{parN}:}{

      Object of class \code{"integer"}: number of parameters.

    }
    \item{\code{parNames}:}{

      Object of class \code{"character"}: vector of parameter names. Its
      length is in slot \code{parN}.

    }
    \item{\code{inputNames}:}{

      Object of class \code{"character"}: names of the inputs used by
      the kernel.

    }
    \item{\code{topParN}:}{

      Object of class \code{"integer"}: number of \emph{top} parameters.

    }
    \item{\code{topParNames}:}{

      Object of class \code{"character"}. Names of the top parameters.

    }
    \item{\code{topPar}:}{

      Object of class \code{"numeric"}. Values of the top parameters.

    }
    \item{\code{topParLower}:}{

      Object of class \code{"numeric"}. Lower bounds for the top
      parameters.

    }
    \item{\code{topParUpper}:}{

      Object of class \code{"numeric"}. Upper bounds for the top
      parameters.

    }
    \item{\code{parsedFormula}:}{

      Object of class \code{"list"}. Ugly draft for some slots to be
      added in the next versions.

    }

  }

}
\section{Extends}{
  
  Class \code{"\linkS4class{covAll}"}, directly.
  
}
\section{Methods}{
  
  \describe{
    
    \item{as.list}{
      
      \code{signature(object = "covComp")}: coerce \code{object} into a
      list of covariance kernels, each inheriting from the virual class
      \code{"covAll"}. This is useful e.g., to extract the coefficients
      or to plot a covariance component.
      
    }
    \item{checkX}{
      
      \code{signature(object = "covComp", X = "data.frame")}: check that
      the inputs exist with suitable column names and suitable \emph{factor}
      content. The levels should match the prescribed levels. Returns a
      matrix with the input columns in the order prescribed by
      \code{object}.
      
    }
    \item{coef, coef<-}{
      
      \code{signature(object = "covComp")}: extract or replace the
      vector of coefficients.
      
    }
    \item{coefLower, coefUpper}{
      
      \code{signature(object = "covComp")}: extract the vector of Lower
      or Upper bounds on the coefficients.
      
    }
    \item{scores}{

      \code{signature(object = "covComp")}: return the vector of
      scores, i.e. the derivative of the log-likelihood w.r.t. the
      parameter vector at the current parameter values.

    }

  }
  

  
}
% \references{}
%% \author{}
%\note{}

\seealso{

  The \code{\link{covComp}} creator.
  
}
\examples{

showClass("covComp")

}
\keyword{classes}
