\name{plot.optim_nmsa}
\alias{plot.optim_nmsa}

\title{
Plot an optim_nmsa Object
}
\description{
Creates convergence or contour plots for visual inspection of the optimization result. Note that 'trace' must be activated for this function. \cr
In case of a bivariate optimization, the 'contour' plot gives an overview of the parameter development over time in the entire state space. This is useful for the evaluation of the algorithm settings and therefore helps improving the performance. The development of the response can be visualized via the 'convergence' plot.
}
\usage{
\method{plot}{optim_nmsa}(x, type = 'convergence', lower = NA, upper = NA, ...)
}
\arguments{
  \item{x}{
  Object of type 'optim_nmsa' to be plotted. The 'trace' entry must not be empty.
  }
  \item{type}{
  Character string which determines the plot type. Either 'convergence' or 'contour' is possible.
  }
  \item{lower}{
  Vector containing the lower limits of the variables in the plot. Only useful for 'contour' plots.
  }
  \item{upper}{
  Vector containing the upper limits of the variables in the plot. Only useful for 'contour' plots.
  }
  \item{...}{
  Further arguments for the generic plot function.
  }
}

\author{
Kai Husmann, Alexander Lange
}
\seealso{
\code{\link{optim_nm}}, \code{\link{optim_sa}}
}

\examples{
# S3 method for class 'optim_nlme'

# Himmelblau's function
hi <- function(x){(x[1]**2 + x[2] - 11)**2 + (x[1] + x[2]**2 -7)**2}

out_nm <- optim_nm(hi, k = 2, trace = TRUE)
out_sa <- optim_sa(fun = hi, start = c(runif(2, min = -1, max = 1)),
                   trace = TRUE, lower = c(-4, -4) ,upper=c(4, 4),
                   control = list(t0 = 1000, nlimit = 1500,r = 0.8))

# Examples for optimization results via 'Nelder-Mead' method.
plot(out_nm)
plot(out_nm, type = "contour", lower = c(-4, -4), upper = c(4, 4))

# Examples for optimization results via 'Simulated Annealing' method.
plot(out_sa)
plot(out_sa, type = "contour")
}
