\name{varVec-methods}
\docType{methods}
%% \alias{varVec-methods}
\alias{varVec,covMan-method}
%% \alias{varVec,covManVec-method}
\alias{varVec,covTS-method}
\title{
  Covariance Matrix for a Covariance Kernel Object
}
\description{
  
  Covariance matrix for a covariance kernel object.

}
\usage{
\S4method{varVec}{covMan}(object, X, compGrad = FALSE, 
       checkNames = NULL, index = -1L, ...)

\S4method{varVec}{covTS}(object, X, compGrad = FALSE, 
       checkNames = TRUE, index = -1L, ...)
}
\arguments{
  \item{object}{

    An object with S4 class corresponding to a covariance kernel.
    
  }
  \item{X}{

    The usual matrix of spatial design points, with \eqn{n} rows and
    \eqn{d} cols where \eqn{n} is the number of spatial points and
    \eqn{d} is the 'spatial' dimension.

  }
  \item{compGrad}{

    Logical. If \code{TRUE} a derivative with respect to the vector of
    parameters will be computed and returned as an attribute of the
    result. For the \code{covMan} class, this is possible only when the
    gradient of the kernel is computed and returned as a
    \code{"gradient"} attribute of the result.

  }
  \item{checkNames}{

    Logical. If \code{TRUE} (default), check the compatibility of
    \code{X} with \code{object}, see \code{\link{checkX}}.

  }
  \item{index}{

    Integer giving the index of the derivation parameter in the official
    order.

  }
  \item{...}{

    Not used yet.

  } 
}
\details{
  
   The variance vector is computed in a C program using the \code{.Call}
   interface. The R kernel function is evaluated within the C code using
   \code{eval}.

 }
 \value{
   
  A vector of length \code{nrow(X)} with general element \eqn{V_{i} :=
  K(\mathbf{x}_{i},\,\mathbf{x}_{i};\,\boldsymbol{\theta})}{V[i]
  = K(X[i, ], X[i, ], \theta)} where
  \eqn{K(\mathbf{x}_1,\,\mathbf{x}_2;\,\boldsymbol{\theta})}{K(x1, x2,
  \theta)} is the covariance kernel function.

}
%%\references{}
%% \author{
  
%%   Y. Deville, O. Roustant, D. Ginsbourger, N. Durrande.

%% }
\note{

  The value of the parameter \eqn{\boldsymbol{\theta}}{\theta} can be
  extracted from the object with the \code{coef} method.

}
%%\section{Caution}{}
\seealso{

  \code{\link{coef}} method

}
\examples{
myCov <- covTS(inputs = c("Temp", "Humid", "Press"),
               kernel = "k1PowExp",
               dep = c(range = "cst", shape = "cst"),
               value = c(shape = 1.8, range = 1.1))
n <- 100; X <- matrix(runif(n*3), nrow = n, ncol = 3)
try(V1 <- varVec(myCov, X = X)) ## bad colnames
colnames(X) <- inputNames(myCov)
V2 <- varVec(myCov, X = X)

Xnew <- matrix(runif(n * 3), nrow = n, ncol = 3)
colnames(Xnew) <- inputNames(myCov)
V2 <- varVec(myCov, X = X)

}
%%\keyword{ ~kwd1 }
%%\keyword{ ~kwd2 }
