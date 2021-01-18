\name{covMat-methods}
\docType{methods}
%% \alias{covMat-methods}
\alias{covMat,covMan-method}
%% \alias{covMat,covManVec-method}
\alias{covMat,covTS-method}
\title{
  Covariance Matrix for a Covariance Kernel Object
}
\description{
  
  Covariance matrix for a covariance kernel object.

}
\usage{
\S4method{covMat}{covMan}(object, X, Xnew, compGrad = hasGrad(object), 
       checkNames = NULL, index = 1L, ...)

\S4method{covMat}{covTS}(object, X, Xnew, compGrad = FALSE, 
       checkNames = TRUE, index = 1L, ...)
}
\arguments{
  \item{object}{

    An object with S4 class corresponding to a covariance kernel.
    
  }
  \item{X}{

    The matrix (or data.frame) of design points, with \eqn{n} rows and
    \eqn{d} cols where \eqn{n} is the number of spatial points and
    \eqn{d} is the 'spatial' dimension.

  }
  \item{Xnew}{

    An optional new matrix of spatial design points. If missing, the
    same matrix is used: \code{Xnew = X}.

  }
  \item{compGrad}{

    Logical. If \code{TRUE} a derivative with respect to a parameter
    will be computed and returned as an attribute of the result. For
    the \code{covMan} class, this is possible only when the gradient
    of the kernel is computed and returned as a \code{"gradient"}
    attribute of the result.

  }
  \item{checkNames}{

    Logical. If \code{TRUE} (default), check the compatibility of
    \code{X} with \code{object}, see \code{\link{checkX}}.

  }
  \item{index}{

    Integer giving the index of the derivation parameter in the official
    order. Ignored if \code{compGrad = FALSE}.

  }
  \item{...}{

    not used yet.

  } 
}
\details{
  
   The covariance matrix is computed in a C program using the
   \code{.Call} interface. The R kernel function is evaluated within the
   C code using \code{eval}.

 }
 \value{
   
  A \eqn{n_1 \times n_2}{(n_1, n_2)} matrix with general element
  \eqn{C_{ij} :=
  K(\mathbf{x}_{1,i},\,\mathbf{x}_{2,j};\,\boldsymbol{\theta})}{C[i, j] = K(X1[i,
  ], X2[j, ], \theta)} where
  \eqn{K(\mathbf{x}_1,\,\mathbf{x}_2;\,\boldsymbol{\theta})}{K(x1, x2,
  \theta)} is the covariance kernel function.

}
%%\references{}
\author{
  
  Y. Deville, O. Roustant, D. Ginsbourger, N. Durrande.

}
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
try(C1 <- covMat(myCov, X)) ## bad colnames
colnames(X) <- inputNames(myCov)
C2 <- covMat(myCov, X)

Xnew <- matrix(runif(n * 3), nrow = n, ncol = 3)
colnames(Xnew) <- inputNames(myCov)
C2 <- covMat(myCov, X, Xnew)

## check with the same matrix in 'X' and 'Xnew'
CMM <- covMat(myCov, X, X)
CM <- covMat(myCov, X)
max(abs(CM - CMM))
}
%%\keyword{ ~kwd1 }
%%\keyword{ ~kwd2 }
