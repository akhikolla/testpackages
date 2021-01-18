\name{covQualNested-class}
\Rdversion{1.1}
\docType{class}
\alias{covQualNested-class}

\title{Class \code{"covQualNested"}}
\description{
  
  Correlation or covariance structure for qualitative inputs
  (i.e. factors) obtained by nesting.

}
\section{Objects from the Class}{
  
  Objects can be created by calls of the form \code{new("covQualNested",
    ...)}.
  
}
\section{Slots}{
  \describe{
    \item{\code{covLevels}:}{

      Object of class \code{"function"} computing the covariance matrix
      for the set of all levels.

    }
    \item{\code{covLevMat}:}{

      Object of class \code{"matrix"}. The matrix returned by the
      function in slot \code{covLevels}. Since this matrix is often
      needed, it can be stored rather than recomputed.

    }
    \item{\code{hasGrad}:}{

      Object of class \code{"logical"}. If \code{TRUE}, the analytical
      gradient can be computed.

    }
    \item{\code{acceptLowerSQRT}:}{

      Object of class \code{"logical"}. If \code{TRUE}, the lower square
      root of the matrix can be returned

    }
    \item{\code{label}:}{

      Object of class \code{"character"}. A label to describe the
      kernel.

    }
    \item{\code{d}:}{

      Object of class \code{"integer"}. The number of inputs.

    }
    \item{\code{inputNames}:}{
      
      Object of class \code{"character"} Names of the inputs.

    }
    \item{\code{nlevels}:}{
      
      Object of class \code{"integer"} with length \code{d} give the
      number of levels for the factors.

    }
    \item{\code{levels}:}{
      
      Object of class \code{"list"} with length \code{d}. Gives the
      levels for the inputs.

    }
    \item{\code{parLower}:}{
      
      Object of class \code{"numeric"}. Lower bounds on the (hyper)
      parameters.

    }
    \item{\code{parUpper}:}{

      Object of class \code{"numeric"}. Upper bounds on the (hyper)
      parameters.

    }
    \item{\code{par}:}{

      Object of class \code{"numeric"}. Value of the (hyper) parameters.

    }
    \item{\code{parN}:}{

      Object of class \code{"integer"}. Number of (hyper) parameters.

    }
    \item{\code{kernParNames}:}{

      Object of class \code{"character"}. Name of the parameters.

    }
    \item{\code{group}:}{

      Object of class \code{"integer"}. Group numbers: one for each final
      level.

    }
    \item{\code{groupLevels}:}{
      
      Object of class \code{"character"}. Vector of labels for the
      groups.

    }
    \item{\code{between}:}{

      Object of class \code{"covQual"}. A covariance or correlation
      structure that can be used between groups.

    }
    \item{\code{within}:}{

      Object of class \code{"list"}. A list of covariance or correlation
      structures that are used within the groups. Each item has class
      \code{"covQual"}.

    }
    \item{\code{parNCum}:}{

      Object of class \code{"integer"}. Cumulated number of
      parameters. Used for technical computations.

    }
    \item{\code{contrasts}:}{

      Object of class \code{"function"}. A contrast function like
      \code{contr.helmod}. This function must return a contrast matrix
      with columns having unit norm.

    }
  }
}
\section{Extends}{
  Class \code{"\linkS4class{covQual}"}, directly.
  Class \code{"\linkS4class{covAll}"}, by class "covQual", distance 2.
}
\section{Methods}{
  
  No methods defined with class "covQualNested" in the signature.

}
%% \references{}
%% \note{}
%% \seealso{}
\examples{
showClass("covQualNested")
}
\keyword{classes}
