\name{covMan}
\alias{covMan}

\title{ Creator Function for \code{covMan} Objects }
 
\description{
  
  Creator function for \code{covMan} objects representing a covariance
  kernel entered manually.

}

\usage{
   covMan(kernel, hasGrad = FALSE, acceptMatrix = FALSE,
          inputs = paste("x", 1:d, sep = ""), 
          d = length(inputs), parNames,
          par = NULL, parLower = NULL, parUpper = NULL,
          label = "covMan", ...) 
} 

\arguments{
  \item{kernel}{
    
    A (semi-)positive definite function. This must be an object of class
    \code{"function"} with formal arguments named \code{"x1"}, \code{"x2"}
    and \code{"par"}. The first two formal arguments are locations vectors
    or matrices. The third formal is for the vector
    \eqn{\boldsymbol{\theta}}{\theta} of \emph{all} covariance
    parameters. An analytical gradient can be computed and returned as an
    attribute of the result with name \code{"gradient"}. See
    \bold{Details}.

  }
  \item{hasGrad}{
    
    Logical indicating whether the \code{kernel} function returns the
    gradient w.r.t. the vector of parameters as a \code{"gradient"}
    attribute of the result. See \bold{Details}
    
  }
  \item{acceptMatrix}{
    
    Logical indicating whether \code{kernel} admits matrices as
    arguments. Default is \code{FALSE}. See \bold{Examples} below.

  }
  \item{inputs}{

    Character vector giving the names of the inputs used as arguments
    of \code{kernel}. Optional if \code{d} is given.

  }
  \item{d}{

    Integer specifying the spatial dimension (equal to the number of
    inputs). Optional if \code{inputs} is given.

  }
  \item{parNames}{

    Vector of character strings containing the parameter names.

  }
  \item{par, parLower, parUpper}{

    Optional numeric vectors containing the parameter values, lower
    bounds and upper bounds.

  }
  \item{label}{

    Optional character string describing the kernel.

  }
  \item{...}{

    Not used at this stage.

  }

}
\details{
  
  The formals and the returned value of the \code{kernel} function
  must be in accordance with the value of \code{acceptMatrix}. 

  \itemize{
    \item{

      When \code{acceptMatrix} is \code{FALSE}, the formal arguments
      \code{x1} and \code{x2} of \code{kernel} are numeric vectors with
      length \eqn{d}. The returned result is a numeric vector of length
      \eqn{1}.  The attribute named \code{"gradient"} of the returned
      value (if provided in accordance with the value of \code{hasGrad})
      must then be a numeric vector with length equal to the number of
      covariance parameters. It must contain the derivative of the
      kernel value
      \eqn{K(\mathbf{x}_1,\,\mathbf{x}_2;\,\boldsymbol{\theta})}{K(x1,
	x2, \theta)} with respect to the parameter vector
      \eqn{\boldsymbol{\theta}}{\theta}.
      
    }
    \item{
      
      When \code{acceptMatrix} is \code{TRUE}, the formals \code{x1} and
      \code{x2} are matrices with \eqn{d} columns and with \eqn{n_1}{n1}
      and \eqn{n_2}{n2} rows. The result is then a covariance matrix
      with \eqn{n_1}{n1} rows and \eqn{n_2}{n2} columns. The gradient
      attribute (if provided in accordance with the value of
      \code{hasGrad}) must be a list with length equal to the number of
      covariance parameters. The list element \eqn{\ell}{l} must contain
      a numeric matrix with dimension \eqn{(n_1, n_2)}{c(n_1, n_2)}
      which is the derivative of the covariance matrix w.r.t. the
      covariance parameter \eqn{\theta_\ell}{\theta[l]}. }
    
  }

}
\note{

  The kernel function must be symmetric with respect to its
  first two arguments, and it must be positive definite,
  which is not checked. If the function returns an object
  with a \code{"gradient"} attribute but \code{hasGrad} was
  set to \code{FALSE}, the gradient will \emph{not} be used
  in optimization.

  The name of the class was motivated by earlier stages in the
  development.

}

\examples{
myCovMan <- 
      covMan(
         kernel = function(x1, x2, par) { 
         htilde <- (x1 - x2) / par[1]
         SS2 <- sum(htilde^2)
         d2 <- exp(-SS2)
         kern <- par[2] * d2
         d1 <- 2 * kern * SS2 / par[1]            
         attr(kern, "gradient") <- c(theta = d1,  sigma2 = d2)
         return(kern)
      },
      hasGrad = TRUE,    
      d = 1,
      label = "myGauss",
      parLower = c(theta = 0.0, sigma2 = 0.0),
      parUpper = c(theta = Inf, sigma2 = Inf),
      parNames = c("theta", "sigma2"),
      par = c(NA, NA)
      )
      
# Let us now code the same kernel in C
kernCode <- "
       SEXP kern, dkern;
       int nprotect = 0, d;
       double SS2 = 0.0, d2, z, *rkern, *rdkern;

       d = LENGTH(x1);
       PROTECT(kern = allocVector(REALSXP, 1)); nprotect++;
       PROTECT(dkern = allocVector(REALSXP, 2)); nprotect++;
       rkern = REAL(kern);
       rdkern = REAL(dkern);

       for (int i = 0; i < d; i++) {
         z = ( REAL(x1)[i] - REAL(x2)[i] ) / REAL(par)[0];
         SS2 += z * z; 
       }

       d2 = exp(-SS2);
       rkern[0] = REAL(par)[1] * d2;
       rdkern[1] =  d2; 
       rdkern[0] =  2 * rkern[0] * SS2 / REAL(par)[0];

       SET_ATTR(kern, install(\"gradient\"), dkern);
       UNPROTECT(nprotect);
       return kern;
     "
     
myCovMan  

## "inline" the C function into an R function: much more efficient! 

\dontrun{
require(inline)
kernC <- cfunction(sig = signature(x1 = "numeric", x2 = "numeric",
                                   par = "numeric"),
                    body = kernCode)
myCovMan <- covMan(kernel = kernC, hasGrad = TRUE, d = 1,
                   parNames = c("theta", "sigma2"))
myCovMan
}

## A kernel admitting matricial arguments
myCov <- covMan(
    
    kernel = function(x1, x2, par) { 
      # x1 : matrix of size n1xd
      # x2 : matrix of size n2xd
            
      d <- ncol(x1)
            
      SS2 <- 0
      for (j in 1:d){
        Aj <- outer(x1[, j], x2[, j], "-")
        Aj2 <- Aj^2
        SS2 <- SS2 + Aj2 / par[j]^2
      }
      D2 <- exp(-SS2)
      kern <- par[d + 1] * D2
    },
    acceptMatrix = TRUE,
    d = 2,
    label = "myGauss",
    parLower = c(theta1 = 0.0, theta2 = 0.0, sigma2 = 0.0),
    parUpper = c(theta1 = Inf, theta2 = Inf, sigma2 = Inf),
    parNames = c("theta1", "theta2", "sigma2"),
    par = c(NA, NA, NA)

)
      
coef(myCov) <- c(0.5, 1, 4)
show(myCov)

## computing the covariance kernel between two points
X <- matrix(c(0, 0), ncol = 2)
Xnew <- matrix(c(0.5, 1), ncol = 2)
colnames(X) <- colnames(Xnew) <- inputNames(myCov)
covMat(myCov, X)            ## same points
covMat(myCov, X, Xnew)      ## two different points

# computing covariances between sets of given locations
X <- matrix(c(0, 0.5, 0.7, 0, 0.5, 1), ncol = 2)
t <- seq(0, 1, length.out = 3)
Xnew <- as.matrix(expand.grid(t, t))
colnames(X) <- colnames(Xnew) <- inputNames(myCov)
covMat(myCov, X)         ## covariance matrix
covMat(myCov, X, Xnew)   ## covariances between design and new data

}