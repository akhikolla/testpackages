\name{covTS}
\alias{covTS}
\encoding{UTF8}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Creator Function for \code{covTS} Objects
}
\description{
  
  Creator function for \code{covTS} objects representing a Tensor Sum
  covariance kernel.

}
\usage{
covTS(inputs = paste("x", 1:d, sep = ""),
      d = length(inputs), kernel = "k1Matern5_2",
      dep = NULL, value = NULL, var = 1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{inputs}{

    Character vector giving the names of the inputs used as arguments of
    \code{kernel}. Optional if \code{d} is given.

  }
  \item{d}{ 

    Integer specifying the spatial dimension (equal to the number of
    inputs). Optional if \code{inputs} is given.

  }
  \item{kernel}{
    
    Character, name of the one-dimensional kernel.

  }
  \item{dep}{

    Character vector with elements \code{"cst"} or \code{"input"}
    usually built using the concatenation \code{\link{c}}. The names
    must correspond to parameters of the kernel specified with
    \code{kernel}. When an element is \code{"cst"}, the corresponding
    parameter of the 1d kernel will be the same for all inputs.  When
    the element is \code{"input"}, the corresponding parameter of the 1d
    kernel gives birth to \code{d} parameters in the \code{covTS}
    object, one by input.

  }
  \item{value}{
    
    Named numeric vector. The names must correspond to the 1d kernel
    parameters.

  }
  \item{var}{

    Numeric vector giving the variances \eqn{\sigma^2_i}{\sigma^2[i]}
    that weight the \eqn{d} components.

  }
  \item{\dots}{

    Not used at this stage.

  }
}

\details{

  A \code{covTS} object represents a \eqn{d}-dimensional kernel object
  \eqn{K} of the form \deqn{K(\mathbf{x}, \mathbf{x}';
  \boldsymbol{\theta}) = \sum_{i=1}^d k(x_i,
  x_i';\boldsymbol{\theta}_{\mathbf{s}_i})}{ K(x, x'; \theta) =
  \sum_{i=1}^d k(x_i, x'_i; \theta[s[i, ]])} where \eqn{k}{k} is
  the covariance kernel for a Gaussian Process \eqn{Y_x} indexed by a scalar
  \eqn{x}. The \eqn{d} numbers \eqn{x_i} stand for the components of the
  \eqn{d}-dimensional location vector \eqn{\mathbf{x}}{x}. The length
  \eqn{p} of all the vectors \eqn{\mathbf{s}_i}{s[i, ]} is the number of
  parameters of the one-dimensional kernel \eqn{k}, i.e. 2 or 3 for
  classical covariance kernels.

  The package comes with the following covariance kernels which can
  be given as \code{kernel} argument.
  
  \tabular{llll}{

    \emph{name} \tab \emph{description} \tab \eqn{p} \tab \emph{par. names} \cr

    \code{k1Exp} \tab exponential \tab \eqn{2} \tab \code{range},
    \code{var}\cr

    \code{k1Matern3_2} \tab Matérn \eqn{\nu = 3/2} \tab \eqn{2} \tab
    \code{range}, \code{var} \cr

    \code{k1Matern5_2} \tab Matérn \eqn{\nu = 5/2}\tab \eqn{2} \tab
    \code{range}, \code{var} \cr 

    \code{k1PowExp} \tab power exponential \tab \eqn{3} \tab \code{range},
    \code{shape}, \code{var} \cr

    \code{k1Gauss} \tab gaussian or "square exponential" \tab \eqn{2} \tab
    \code{range}, \code{var} \cr

  }

  Note that the exponential kernel of \code{k1Exp} is identical to the
  Matérn kernel for \eqn{\nu = 1/2}, and that the three Matérns kernels
  provided here for \eqn{\nu = 1/2}, \eqn{\nu = 3/2} and \eqn{\nu = 5/2}
  are special cases of Continuous AutoRegressive (CAR) process
  covariances, with respective order \eqn{1}, \eqn{2} and \eqn{3}.
 
}
\value{

  An object with S4 class \code{"covTS"}.

}
\references{
  
  N. Durrande, D. Ginsbourger, and O. Roustant (2012) Additive
  "Covariance kernels for high-dimensional Gaussian Process modeling",
  \emph{Annales de la Faculté des Sciences de Toulouse} 21(3),
  pp. 481–499.

}
\author{

  Y. Deville, O. Roustant D. Ginsbourger

}
\note{

  Most arguments receive default values or are recycled if necessary.

}

\section{Caution}{
  
  The \eqn{1d} kernel \eqn{k} as given in \code{kernel} is always
  assumed to have a variance parameter with name \code{var}. This
  assumption may be relaxed in future versions.

}
%%\seealso{}
\examples{
myCov1 <- covTS(kernel = "k1Exp", inputs = c("v1", "v2", "v3"),
                dep = c(range = "input"))
coef(myCov1) <- c(range = c(0.3, 0.7, 0.9), sigma2 = c(2, 2, 8))

myCov1
coef(myCov1)
coef(myCov1, as = "matrix")
coef(myCov1, as = "list")
coef(myCov1, as = "matrix", type = "range")

# with a common range parameter
myCov2 <- covTS(kernel = "k1Exp", inputs = c("v1", "v2", "v3"), 
                dep = c(range = "cst"), value = c(range = 0.7),
                var = c(2, 2, 8))
myCov2

myCov3 <- covTS(d = 3, kernel = "k1PowExp",
                dep = c(range = "cst", shape = "cst"),
                value = c(shape = 1.8, range = 1.1),
                var = c(2, 2, 8))
myCov3
}
