\name{covTS-class}
\Rdversion{1.1}
\docType{class}
\alias{covTS-class}
%% mind the methods written for the class union covAll
%% \alias{checkX,covTS,matrix-method}
%% \alias{coef,covTS-method}  % see coef.covTS.Rd
\alias{coef<-,covTS,numeric-method}
\alias{coefLower,covTS-method}
\alias{coefLower<-,covTS-method}
\alias{coefUpper,covTS-method}
\alias{coefUpper<-,covTS-method}
%% \alias{covMat,covTS-method}  % see covMat.covTS.Rd
%\alias{inputNames,covTS-method}
\alias{kernelName,covTS-method}
%% \alias{npar,covTS-method}
%% \alias{parMap,covTS-method}
%% \alias{simulPar,covTS-method}
\alias{scores,covTS-method}
\alias{show,covTS-method}

\title{Class \code{"covTS"}}
\description{
  
  S4 class representing a Tensor Sum (TS) covariance kernel.

}
\section{Objects from the Class}{

  Objects can be created by call of the form \code{new("covTS", ...)} or
  by using the \code{\link{covTS}} function.
  %%  ~~ describe objects here ~~ 
  
}
\section{Slots}{
  \describe{
    \item{\code{d}:}{

      Object of class \code{"integer"}, the spatial dimension or number
      of inputs of the covariance.

    }
    \item{\code{inputNames}:}{


      Object of class \code{"character"}, vector of input names. Length
      \code{d}.

    }
    \item{\code{kernel}:}{

      Object of class \code{"covMan"} representing a 1d kernel.

    }
    \item{\code{kernParNames}:}{

      Object of class \code{"character"}, name of the kernel (among the
      allowed ones).

    }
    \item{\code{kernParCodes}:}{

      Object of class \code{"integer"}, an integer code stating the
      dependence of the parameter to the input.

    }
    \item{\code{par}:}{

      Object of class \code{"numeric"}, numeric vector of parameter
      values.

    }
    \item{\code{parN}:}{
      
      Object of class \code{"integer"}, total number of parameters.

    }
    \item{\code{parInput}:}{

      Object of class \code{"integer"}, the number of the inputs for
      each parameter. Same length as \code{par}, values between \code{1}
      and \code{d}.

    }
    \item{\code{parLower}:}{, }
    \item{\code{parUpper}:}{

      Object of class \code{"numeric"} numeric, vector of (possibly
      infinite) lower/upper bounds on parameters.

    }
    \item{\code{parBlock}:}{

      Object of class \code{"integer"}

    }
  }
}
\section{Methods}{
  \describe{
    %    \item{checkX}{\code{signature(object = "covTS", X = "matrix")}:
    %      check the compatibility of a \code{covTS} object and a design matrix. }
    \item{coef}{
      
      \code{signature(object = "covTS")}: extracts the numeric values of
      the covariance parameters.
      %an ersatz of the \code{covMatrix} method, with slightly different
      %parameters and use.

    }
    \item{coef<-}{

      \code{signature(object = "covTS")}: replaces the whole vector of
      coefficients, as required during ML estimation.

    }
    \item{coefLower}{

      \code{signature(object = "covTS")}: extracts the numeric values of
      the lower bounds.

    }
    \item{coefLower<-}{

      \code{signature(object = "covTS")}: replacement method for lower
      bounds on covTS coefficients.

    }
    \item{coefUpper}{

      \code{signature(object = "covTS")}: ...

    }
    \item{coefUpper<-}{

      \code{signature(object = "covTS")}: replacement method for upper
      bounds on covTS coefficients.

    }
    \item{covMat}{

      \code{signature(object = "covTS")}: builds the covariance matrix,
      or the cross covariance matrix between two sets of locations for a \code{covTS} object.

    }
    \item{kernelName}{

      \code{signature(object = "covTS")}: return the character value of the kernel name.

    }
    \item{parMap}{

      \code{signature(object = "covTS")}: an integer matrix used to map
      the \code{covTS} parameters on the inputs and kernel parameters
      during the computations.

    }
    \item{scores}{

      \code{signature(object = "covTS")}: computes the scores.

    }
    \item{show}{

      \code{signature(object = "covTS")}: prints in a custom format.

    }
    \item{simulPar}{

      \code{signature(object = "covTS")}: simulates random values for
      the covariance parameters.

    }
  }
}
%%\references{}
\author{
  
  Y. Deville, O. Roustant, D. Ginsbourger.

}
\note{

  The names of the methods strive to respect a
  \href{http://en.wikipedia.org/wiki/CamelCase}{camelCase} naming
  convention.

  While the \code{coef<-} replacement method is typically intended
  for internal use during likelihood maximization, the \code{coefLower<-}
  and \code{coefUpper<-} replacement methods can be used when some
  rough information exists on the possible values of the parameters.

}
\seealso{

  The \code{\link{covTS}} function providing  a creator. 

}
\examples{
showClass("covTS")
}
\keyword{classes}
