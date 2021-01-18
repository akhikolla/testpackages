\name{covOrd-class}
\Rdversion{1.1}
\docType{class}
\alias{covOrd-class}
\alias{checkX,covOrd,matrix-method}
\alias{checkX,covOrd,data.frame-method}
\alias{coef,covOrd-method}
\alias{coef<-,covOrd,numeric-method}
\alias{coefLower<-,covOrd-method}
\alias{coefLower,covOrd-method}
\alias{coefUpper<-,covOrd-method}
\alias{coefUpper,covOrd-method}
\alias{coerce,covOrd,function-method}
\alias{covMat,covOrd-method}
%% \alias{covMat,covMan-method}  % see covMat.covMan.Rd
\alias{npar,covOrd-method}
\alias{show,covOrd-method}
\alias{scores,covOrd-method}
\alias{simulate,covOrd-method}
\alias{varVec,covOrd-method}
\title{
  Class \code{"covOrd"}
}
\description{

Covariance kernel for qualitative ordered inputs obtained by warping. 

Let \eqn{u} be an ordered factor with levels \eqn{u_1, \dots, u_L}{u[1], ..., u[L]}. 
Let \eqn{k_1}{k1} be a 1-dimensional stationary kernel (with no or fixed parameters), \eqn{F} a warping function i.e. an increasing function on the interval \eqn{[0,1]} and \eqn{\theta}{theta} a scale parameter. Then \eqn{k} is defined by:
\deqn{k(u_i, u_j) = k_1([F(z_i) - F(z_{j})]/\theta)}{k(u[i], u[j]) = k1([F(z[i]) - F(z[j])]/theta)}
where \eqn{z_1, \dots, z_L}{z[1], ..., z[L]} form a regular sequence from \eqn{0} to \eqn{1} (included). Notice that an example of warping is a distribution function (cdf) restricted to \eqn{[0,1]}.
}

\section{Objects from the Class}{
  
  Objects can be created by calls of the form \code{new("covOrd", ...)}.
  
  
}
\section{Slots}{
  
  \describe{
    \item{\code{covLevels}:}{

      Same as for \code{\link{covQual-class}}. 
      
    }
    \item{\code{covLevMat}:}{

      Same as for \code{\link{covQual-class}}.

    }
    \item{\code{hasGrad}:}{
      Same as for \code{\link{covQual-class}}.
    }

    \item{\code{acceptLowerSQRT}:}{
      Same as for \code{\link{covQual-class}}.
    }
    \item{\code{label}:}{
      Same as for \code{\link{covQual-class}}.
    }
    \item{\code{d}:}{

      Same as for \code{\link{covQual-class}}. Here equal to 1.

    }
    \item{\code{inputNames}:}{

      Same as for \code{\link{covQual-class}}.

    }
    \item{\code{nlevels}:}{

      Same as for \code{\link{covQual-class}}.

    }
    \item{\code{levels}:}{

      Same as for \code{\link{covQual-class}}.

    }
    \item{\code{parLower}:}{

      Same as for \code{\link{covQual-class}}.
    }

    \item{\code{parUpper}:}{

      Same as for \code{\link{covQual-class}}.
      
    }
    \item{\code{par}:}{

      Same as for \code{\link{covQual-class}}.

    }
    \item{\code{parN}:}{

      Same as for \code{\link{covQual-class}}.

    }
    \item{\code{kernParNames}:}{

      Same as for \code{\link{covQual-class}}.

    }
    
    \item{\code{k1Fun1}:}{
      A function representing a 1-dimensional stationary kernel function, with no or fixed parameters.
    }
    \item{\code{warpFun}:}{
      A cumulative density function representing a warping.
    }
    \item{\code{cov}:}{
      Object of class \code{"integer"}. The value \code{0L} corresponds
      to a correlation kernel while \code{1L} is for a covariance
      kernel.
    }
    \item{\code{parNk1}:}{
      Object of class \code{"integer"}. Number of parameters of \code{k1Fun1}. Equal to \code{0} at this stage.
    }
    \item{\code{parNwarp}:}{
      Object of class \code{"integer"}. Number of parameters of \code{warpFun}.
    }
    \item{\code{k1ParNames}:}{
      Object of class \code{"character"}. Parameter names of \code{k1Fun1}.
    }
    \item{\code{warpParNames}:}{
      Object of class \code{"character"}. Parameter names of \code{warpFun}.
    }
  }
}

\section{Methods}{
  
  \describe{
    \item{checkX}{

      \code{signature(object = "covOrd", X = "data.frame")}: check that
      the inputs exist with suitable column names and suitable \emph{factor}
      content. The levels should match the prescribed levels. Returns a
      matrix with the input columns in the order prescribed by
      \code{object}.
      
      \code{signature(object = "covOrd", X = "matrix")}: check that the
      inputs exist with suitable column names and suitable \emph{numeric}
      content for coercion into a factor with the prescribed levels.
      Returns a data frame with the input columns in the order
      prescribed by \code{object}.

    }
    \item{coef<-}{
      
      \code{signature(object = "covOrd")}: replace the whole vector of
      coefficients, as required during ML estimation.
      
    }
    \item{coefLower<-}{

      \code{signature(object = "covOrd")}: replacement method for lower
      bounds on covOrd coefficients.

    }
    \item{coefLower}{
      
      \code{signature(object = "covOrd")}: extracts the numeric values of
      the lower bounds.
      
    }
    \item{coef}{
      
      \code{signature(object = "covOrd")}: extracts the numeric values
      of the covariance parameters.
     

    }
    \item{coefUpper<-}{

      \code{signature(object = "covOrd")}: replacement method for upper
      bounds on \code{covOrd} coefficients.

    }
    \item{coefUpper}{

      \code{signature(object = "covOrd")}: ...

    }
    \item{covMat}{

      \code{signature(object = "covOrd")}: build the covariance matrix
      or the cross covariance matrix between two sets of locations for a
      \code{covOrd} object.

    }
    \item{npar}{

      \code{signature(object = "covOrd")}: returns the number of
      parameters.

    }
    \item{scores}{

      \code{signature(object = "covOrd")}: return the vector of
      scores, i.e. the derivative of the log-likelihood w.r.t. the
      parameter vector at the current parameter values.

    }
    \item{simulate}{

      \code{signature(object = "covOrd")}: simulate \code{nsim} paths
      from a Gaussian Process having the covariance structure. The paths
      are indexed by the finite set of levels of factor inputs, and they
      are returned as columns of a matrix.

    }
    \item{varVec}{

      \code{signature(object = "covOrd")}: build the variance vector
      corresponding to a set locations for a \code{covOrd} object.

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

showClass("covOrd")

}
\keyword{classes}
