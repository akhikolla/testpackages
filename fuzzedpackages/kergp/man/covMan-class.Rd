\name{covMan-class}
\Rdversion{1.1}
\docType{class}
\alias{covMan-class}
%\alias{checkX,covMan,matrix-method}
%\alias{coef,covMan-method}
\alias{coef<-,covMan,numeric-method}
\alias{coefLower<-,covMan-method}
\alias{coefLower,covMan-method}
%% \alias{coef,covMan-method}  % see coef.covMan.Rd
\alias{coefUpper<-,covMan-method}
\alias{coefUpper,covMan-method}
\alias{coerce,covMan,function-method}
%% \alias{covMat,covMan-method}  % see covMat.covMan.Rd
%% \alias{npar,covMan-method}
\alias{show,covMan-method}
\alias{scores,covMan-method}

%\alias{kernelName,covMan-method}
%\alias{npar,covMan-method}

\title{Class \code{"covMan"}}
\description{

  S4 class representing a covariance kernel defined manually by a
  (semi-)positive definite function.

}
\section{Objects from the Class}{
  
  Objects can be created by calling \code{new("covMan", ...)}
  or by using the \code{\link{covMan}} function.

}
\section{Slots}{
  \describe{
    \item{\code{kernel}:}{

      object of class \code{"function"} defining the kernel (supposed to
      be (semi-)positive definite).

    }
    \item{\code{hasGrad}:}{

      logical indicating whether \code{kernel} returns the gradient
      (w.r.t. the vector of parameters) as \code{"gradient"} attribute
      of the result.

    }
    
    \item{\code{acceptMatrix}:}{

      logical indicating whether \code{kernel} admits matrix
      arguments. Default is \code{FALSE}.

    }
    \item{\code{label}:}{

      object of class character, typically one or two words, used to
      describe the kernel.

    }
    \item{\code{d}:}{

      object of class \code{"integer"}, the spatial dimension or number
      of inputs of the covariance.

    }
    \item{\code{inputNames}:}{

      object of class \code{"character"}, vector of input names. Length
      \code{d}.

    }
    \item{\code{parLower}:}{, }
    \item{\code{parUpper}:}{

      object of class \code{"numeric"}, vector of (possibly infinite)
      lower/upper bounds on parameters.

    }
    \item{\code{par}:}{

      object of class \code{"numeric"}, numeric vector of parameter
      values.

    }
    \item{\code{parN}:}{

      object of class \code{"integer"}, total number of parameters.

    }
    \item{\code{kernParNames}:}{

      object of class \code{"character"}, name of the kernel
      parameters. 
      %(among the allowed ones).
      
    }
  }
}
\section{Methods}{

  \describe{
    %    \item{checkX}{\code{signature(object = "covMan", X =
    %      "matrix")}: check the compatibility of a \code{covMan} object
    %      and a design matrix. }
    \item{coef<-}{
      
      \code{signature(object = "covMan")}: replace the whole vector of
      coefficients, as required during ML estimation.

    }
    \item{coefLower<-}{

      \code{signature(object = "covMan")}: replacement method for lower
      bounds on covMan coefficients.

    }
    \item{coefLower}{

      \code{signature(object = "covMan")}: extracts the numeric values of
      the lower bounds.

    }
    \item{coef}{

      \code{signature(object = "covMan")}: extracts the numeric values
      of the covariance parameters.
      %an ersatz of the \code{covMatrix} method, with slightly different
      %parameters and use.

    }
    \item{coefUpper<-}{

      \code{signature(object = "covMan")}: replacement method for upper
      bounds on covMan coefficients.

    }
    \item{coefUpper}{

      \code{signature(object = "covMan")}: ...

    }
    \item{covMat}{

      \code{signature(object = "covMan")}: builds the covariance matrix
      or the cross covariance matrix between two sets of locations for a
      \code{covMan} object.

    }
    \item{scores}{

      \code{signature(object = "covMan")}: computes the scores
      (derivatives of the log-likelihood w.r.t. the covariance
      parameters.

    }
    \item{show}{

      \code{signature(object = "covMan")}: prints in a custom format.

    }
%    \item{kernelName}{\code{signature(object = "covMan")}: return the
%      character value of the kernel name.}
  }
}

%%\references{}

\author{

  Y. Deville, O. Roustant, D. Ginsbourger and N. Durrande.

}

\note{
  
  While the \code{coef<-} replacement method is typically intended for
  internal use during likelihood maximization, the \code{coefLower<-}
  and \code{coefUpper<-} replacement methods can be used when some
  information is available on the possible values of the parameters.

}
\seealso{

  The \code{\link{covMan}} function providing a creator.
  
}
\examples{
showClass("covMan")
}
\keyword{classes}
