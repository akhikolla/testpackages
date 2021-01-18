\name{optim_nm}
\alias{optim_nm}
\title{
  Optimization with Nelder-Mead
}
\description{
  This function contains a direct search algorithm, to minimize or maximize an objective function with respect to their
  input parameters.
}
\usage{
  optim_nm(fun, k = 0, start, maximum = FALSE, trace = FALSE,
           alpha = 1, beta = 2, gamma = 1/2, delta = 1/2,
           tol = 0.00001, exit = 500, edge = 1)
}

\arguments{
  \item{fun}{
  Function to minimize or maximize. It should return a single scalar value.
}
  \item{k}{
  Number of parameters of the objective function.
  }
  \item{start}{
  Optional vector with starting values. Number of values must be equal to \code{k}. The initial simplex is constructed around this start vector.
  }
  \item{maximum}{
  Logical. The default is FALSE.
  }
  \item{trace}{
  Logical. If TRUE, interim results are stored. Necessary for the plot function. Default is FALSE.
  }
  \item{alpha}{
  A positive scalar which indicates the size of the reflected simplex. The value 1 leads to a reflected simplex of the same size as the former iteration.
  }
  \item{beta}{
  A positive scalar which indicates the size of the expended simplex. It is usually twice as high as \code{alpha}. It must be higher than \code{alpha}.
   }
  \item{gamma}{
  A positive scalar which indicates the size of either the outside contracted simplex or inside contracted simplex. It
  is usually half as high as \code{alpha}. It must be smaller than \code{alpha}.
  }
  \item{delta}{
  A positive scalar which indicates the size of the shrinked simplex. It is usually half as high as alpha. It must be smaller than \code{alpha}.
  }
  \item{tol}{
  A positive scaler describing the tolerance at which the distances in between the function responses of the simplex vertices are close enough to zero to terminate the algorithm.
  }
  \item{exit}{
  A positive scalar giving the maximum number of iterations the algorithm is allowed to take. It is used to prevent infinite loops. In case of optimizing functions with higher dimensions it is quite likely that the algorithm needs more than 500 iterations. The value should therefore be adjusted to the specific optimization problem.
  }
  \item{edge}{
  A positive scalar providing the edge length of the initial simplex. It is useful to adjust the edge length if the initial guess is close to the global optimum or if the parameter space of the loss function is relatively small.
  }
}
  \details{
  The Nelder-Mead method is a comparatively simple heuristic optimization algorithm. It is, However, useful for relatively simple optimization problems without many local minima and low dimensions(n < 10). Nevertheless, the speed and accuracy are rather useful for simple problems. Moreover, the Nelder-Mead is able to optimize functions without derivatives. The handling of the optimization function is quite easy, because there are only few parameters to adjust.
  }
  \value{
  The output is a nmsa_optim object with following entries:
  \describe{
    \item{\code{par}}{
      Function parameters after optimization.
    }
    \item{\code{function_value}}{
      Function response after optimization.
    }
    \item{\code{trace}}{
      Matrix with interim results. NULL if \code{trace} was not activated.
    }
    \item{\code{fun}}{
      The loss function.
    }
    \item{\code{start}}{
      The initial function parameters.
    }
    \item{\code{lower}}{
      The lower boundaries of the function parameters.
    }
    \item{\code{upper}}{
      The upper boundaries of the function parameters.
    }
    \item{\code{control}}{
      The number of parameters and iterations of the algorithm.
    }

  }
}

\references{
Gao, F. and Han, L. (2012). Implementing the nelder-mead simplex algorithm
with adaptive parameters. Computational Optimization and Applications, 51(1):259
277.

Geiger, C. and Kanzow, C. (1999). Das nelder-mead-verfahren. Numerische
Verfahren zur Loesung unregestrierter Optimierungsaufgaben.

Nelder, J. and Mead, R. (1965). A simplex method for function minimization.
Computer Journal, 7(4).
}

\author{
Alexander Lange
}
\seealso{
\code{\link{optim_sa}}, \code{\link{optim}}, \code{\link{plot.optim_nmsa}}
}
\examples{

#####  Rosenbrock function
# minimum at f(1,1) = 0
   B <- function(x){
    100*(x[2]-x[1]^2)^2+(1-x[1])^2
  }

##### Minimization with an initial guess at c(-2.048, 2.048)
  optim_nm(B, start = c(-2.048, 2.048))

#####  Himmelblau's function
# minimum at f(3,2) = 0
# f(-2.805, -3.1313) = 0
# f(-3.779, -3.283) = 0
#f(3.5844, -1.848) = 0
  H <- function(x){
    (x[1]^2+x[2]-11)^2+(x[1]+x[2]^2-7)^2
  }

##### Minimization with defined number of parameters
  optim_nm(fun = H, k = 2)

##### Colville function with 4 parameters
  co <- function(x){
    x1 <- x[1]
    x2 <- x[2]
    x3 <- x[3]
    x4 <- x[4]

    term1 <- 100 * (x1^2 - x2)^2
    term2 <- (x1 - 1)^2
    term3 <- (x3-1)^2
    term4 <- 90 * (x3^2 - x4)^2
    term5 <- 10.1 * ((x2 - 1)^2 + (x4 - 1)^2)
    term6 <- 19.8 * (x2 - 1)*(x4-1)

     y <- term1 + term2 + term3 + term4 + term5 + term6
  }

  optim_nm(co, k = 4)

#### Minimization with trace
  Output <- optim_nm(H, k = 2, trace = TRUE)
  plot(Output)
  plot(Output, 'contour')

}


