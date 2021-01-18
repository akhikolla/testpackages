\name{covRadial}
\alias{covRadial}
\title{Creator for the Class \code{"covRadial"}}
\usage{

   covRadial(k1Fun1 = k1Fun1Gauss,
             cov = c("corr", "homo"), 
             iso = 0, hasGrad = TRUE,
             inputs = NULL, d = NULL,
             parNames, par = NULL,
             parLower = NULL, parUpper = NULL,
             label = "Radial kernel",
             ...)
}
\arguments{

  \item{k1Fun1}{

    A function of a \emph{scalar} numeric variable, and possibly of an
    extra "shape" parameter. This function should return the first-order
    derivative or the two-first order derivatives as an attibute with
    name \code{"der"} and with a matrix content. When an extra shape
    parameter exists, the gradient should also be returned as an
    attribute with name \code{"gradient"}, see \bold{Examples}
    later. The name of the function can be given as a character string.

  }
  \item{cov}{

    A character string specifying the kind of covariance kernel:
    correlation kernel ("corr") or kernel of a homoscedastic GP ("homo").
    Partial matching is allowed.

  }
  \item{iso}{

    Integer. The value \code{1L} corresponds to an isotropic covariance,
    with all the inputs sharing the same range value.
    
  }
  \item{hasGrad}{

    Integer or logical. Tells if the value returned by the function
    \code{k1Fun1} has an attribute named \code{"der"} giving the
    derivative(s).
    
  }
  \item{inputs}{

    Character. Names of the inputs.

  }
  \item{d}{

    Integer. Number of inputs.

  }
  \item{par, parLower, parUpper}{

    Optional numeric values for the lower bounds on the parameters. Can be
    \code{NA} for \code{par}, can be \code{-Inf} for \code{parLower} and
    \code{Inf} for \code{parUpper}.

  }
  \item{parNames}{

    Names of the parameters. By default, ranges are prefixed
    \code{"theta_"} in the non-iso case and the range is names
    \code{"theta"} in the iso case.
    
  }
  \item{label}{

    A short description of the kernel object.
    
  }
  \item{...}{

    Other arguments passed to the method \code{new}.

  }
}
\value{

  An object with class \code{"covRadial"}.

}
\description{

  Creator for the class \code{"covRadial"}, which describles \emph{radial
    kernels}.

}
\details{

  A radial kernel on the \eqn{d}-dimensional euclidean space
  takes the form 
  \deqn{K(\mathbf{x},\,\mathbf{x}') = \sigma^2 k_1(r)}{K(x, x') =
    sigma^2 * k1(r)}
  where \eqn{k_1(r)}{k1(r)} is a suitable correlation kernel for a
  one-dimensional input, and \eqn{r} is given by
  \deqn{r = \left\{\sum_{\ell = 1}^d [x_\ell - x'_\ell]^2 / \theta_\ell^2
    \right\}^{1/2}.}{r = sigma2 * sqrt(sum((x[l] - x'[l])^2 / theta[l]^2))}


  In this default form, the radial kernel depends on \eqn{d + 1} parameters:
  the \emph{ranges} \eqn{\theta_\ell >0}{theta[l] > 0}  and the
  variance \eqn{\sigma^2}{sigma2}.

  An \emph{isotropic} form uses the same range \eqn{\theta}{theta} for
    all inputs, i.e. sets \eqn{\theta_\ell =
      \theta}{theta[l]=theta} for all \eqn{\ell}{l}. This is obtained
    by using \code{iso = TRUE}.

  A \emph{correlation} version uses  \eqn{\sigma^2 = 1}{sigma2 = 1}. This
  is obtained by using \eqn{cov = "corr"}.

  Finally, the correlation kernel \eqn{k_1(r)}{k1(r)} can depend on a
  "shape" parameter, e.g. have the form \eqn{k_1(r;\,\alpha)}{k1(r,
  alpha)}. The extra shape parameter \eqn{\alpha}{alpha} will be
  considered then as a parameter of the resulting radial kernel, making
  it possible to estimate it by ML along with the range(s) and the
  variance.
  

}
\note{

  When \code{k1Fun1} has more than one formal argument, its arguments
  with position \code{> 1} are assumed to be "shape" parameters of the
  model. Examples are functions with formals
  \code{function(x, shape = 1.0)} or
  \code{function(x, alpha = 2.0, beta = 3.0)},
  corresponding to
  vector of parameter names \code{c("shape")} and
  \code{c("alpha", "beta")}.
  For now only more thant one shape parameter has not been
  tested.
  
  Remind that using a one-dimensional correlation kernel
   \eqn{k_1(r)}{k1(r)} here \emph{does not} warrant that a positive
   semi-definite kernel will result for \emph{any} dimension
   \eqn{d}. This question relates to Schoenberg's theorem and the concept
   of completely monotone functions.
   
} 
\seealso{ 

  \code{\link{k1Fun1Exp}}, \code{\link{k1Fun1Matern3_2}},
  \code{\link{k1Fun1Matern5_2}} or \code{\link{k1Fun1Gauss}} for
  examples of functions that can be used as values for the \code{k1Fun1}
  formal.

}
\references{

  Gregory Fassauher and Michael McCourt (2016) \emph{Kernel-based
  Approximation Methods using MATLAB}. World Scientific.  

}
\examples{
set.seed(123)
d <- 2; ng <- 20
xg <- seq(from = 0, to = 1, length.out = ng)
X <- as.matrix(expand.grid(x1 = xg, x2 = xg))

## ============================================================================
## A radial kernel using the power-exponential one-dimensional
## function
## ============================================================================

d <- 2
myCovRadial <- covRadial(k1Fun1 = k1Fun1PowExp, d = 2, cov = "homo", iso = 1)
coef(myCovRadial)
inputNames(myCovRadial) <- colnames(X)
coef(myCovRadial) <- c(alpha = 1.8, theta = 2.0, sigma2 = 4.0)
y <- simulate(myCovRadial, X = X, nsim = 1)
persp(x = xg, y = xg, z = matrix(y, nrow = ng))

## ============================================================================
## Define the inverse multiquadric kernel function. We return the first two
## derivatives and the gradient as attributes of the result.
## ============================================================================

myk1Fun <- function(x, beta = 2) {
    prov <- 1 + x * x
    res <- prov^(-beta)
    der <- matrix(NA, nrow = length(x), ncol = 2)
    der[ , 1] <- - beta * 2 * x * res / prov
    der[ , 2] <- -2 * beta * (1 - (1 + 2 * beta) * x * x) * res / prov / prov
    grad <- -log(prov) * res
    attr(res, "gradient") <- grad
    attr(res, "der") <- der
    res
}

myCovRadial1 <- covRadial(k1Fun1 = myk1Fun, d = 2, cov = "homo", iso = 1)
coef(myCovRadial1)
inputNames(myCovRadial1) <- colnames(X)
coef(myCovRadial1) <- c(beta = 0.2, theta = 0.4, sigma2 = 4.0)
y1 <- simulate(myCovRadial1, X = X, nsim = 1)
persp(x = xg, y = xg, z = matrix(y1, nrow = ng))


}