\name{optim_sa}
\alias{optim_sa}
\title{
Flexible Optimization with Simulated Annealing
}
\description{
Random search optimization method with systematic component that searches the global optimum. The loss function is allowed to be non-linear, non-differentiable and multimodal. Undefined responses are allowed as well.
}
\usage{
optim_sa(fun, start, maximization = FALSE, trace = FALSE,
         lower, upper, control = list())
}

\arguments{
  \item{fun}{
  Loss function to be optimized. It must return a scalar value. The variables must be assigned as a vector. See 'details'.
  }
  \item{start}{
  Vector of initial values for the function variables. Must be of same length as the variables vector of the loss function. The response of the initial variables combination must be defined (NA or NaN responses are not allowed).
  }
  \item{maximization}{
  Logical. Default is FALSE.
  }
  \item{trace}{
  Logical. If TRUE, interim results are stored. Necessary for the plot function. Default is FALSE.
  }
  \item{lower}{
  Vector of lower boundaries for the function variables. Must be of same length as the variables vector of the function.
  }
  \item{upper}{
  Vector of upper boundaries for the function variables. Must be of same length as the variables vector of the function.
  }
  \item{control}{
  List with optional further arguments to modify the optimization specifically to the loss function:
    \describe{
      \item{\code{vf}}{
        Function that determines the variation of the function variables for the next iteration. The variation function is allowed to depend on the vector of variables of the current iteration, the vector of random factors \code{rf} and the temperature of the current iteration. Default is a uniform distributed random number with relative range \code{rf}.
      }
      \item{\code{rf}}{
        Numeric vector. Random factor vector that determines the variation of the random number of \code{vf} in relation to the dimension of the function variables for the following iteration. Default is 1. If \code{dyn_rf} is enabled, the \code{rf} change dynamically over time.
      }
      \item{\code{dyn_rf}}{
        Logical. \code{rf} change dynamically over time to ensure increasing precision with increasing number of iterations. Default is TRUE, see 'details'.
      }
      \item{\code{t0}}{
        Numeric. Initial temperature. Default is 1000.
      }
      \item{\code{nlimit}}{
        Integer. Maximum number of iterations of the inner loop. Default is 100.
      }
      \item{\code{r}}{
        Numeric. Temperature reduction in the outer loop. Default is 0.6.
      }
      \item{\code{k}}{
        Numeric. Constant for the Metropolis function. Default is 1.
      }
      \item{\code{t_min}}{
        Numeric. Temperature where outer loop stops. Default is 0.1.
      }
      \item{\code{maxgood}}{
        Integer. Break criterion to improve the algorithm performance. Maximum number of loss function improvements in the inner loop. Breaks the inner loop. Default is 100.
      }
      \item{\code{stopac}}{
        Integer. Break criterion to improve the algorithm performance. Maximum number of repetitions where the loss improvement is lower than \code{ac_acc}. Breaks the inner loop. Default is 30.
      }
      \item{\code{ac_acc}}{
        Numeric. Accuracy of the \code{stopac} break criterion in relation to the response. Default is 1/10000 of the function value at initial variables combination.
      }
    }
  }

}
\details{
Simulated Annealing is an optimization algorithm for solving complex functions that may have several optima. The method is composed of a random and a systematic component. Basically, it randomly modifies the variables combination \code{n_limit} times to compare their response values. Depending on the temperature and the constant \code{k}, there is also a likelihood of choosing variables combinations with worse response. There is thus a time-decreasing likelihood of leaving local opima. The Simulated Annealing Optimization method is therefore advantageous for multimodal functions. Undefined response values (NA) are allowed as well. This can be useful for loss functions with variables restrictions. The high number of parameters allows a very flexible parameterization. \code{optim_sa} is able to solve mathematical formulas as well as complex rule sets.

The performance therefore highly depends on the settings. It is indispensable to parameterize the algorithm carefully. The control list is pre-parameterized for loss functions of medium complexity. To improve the performance, the settings should be changed when solving relatively simple functions (e. g. three dimensional multimodal functions). For complex functions the settings should be changed to improve the accuracy. Most important parameters are \code{nlimit},  \code{r} and \code{t0}.

The dynamic \code{rf} adjustment depends on the number of loss function calls which are out of the variables boundaries as well as the temperature of the current iteration. The obligatory decreasing \code{rf} ensures a relatively wide search grid at the beginning of the optimization process that shrinks over time. It thus automatically adjusts for the trade-off between range of the search grid and accuracy. See Pronzato (1984) for more details. It is sometimes useful to disable the dynamic \code{rf} changing when the most performant \code{rf} are known. As \code{dyn_rf} usually improves the performance as well as the accuracy, the default is TRUE.
}
\value{
The output is a nmsa_optim list object with following entries:
  \describe{
    \item{\code{par}}{
      Function variables after optimization.
    }
    \item{\code{function_value}}{
      Loss function response after optimization.
    }
    \item{\code{trace}}{
      Matrix with interim results. NULL if \code{trace} was not activated.
    }
    \item{\code{fun}}{
      The loss function.
    }
    \item{\code{start}}{
      The initial function variables.
    }
    \item{\code{lower}}{
      The lower boundaries of the function variables.
    }
    \item{\code{upper}}{
      The upper boundaries of the function variables.
    }
    \item{\code{control}}{
      Control arguments, see 'details'.
    }


  }
}

\references{
Corana, A., Marchesi, M., Martini, C. and Ridella, S. (1987), Minimizing Multimodal Functions of Continuous Variables with the 'Simulated Annealing' Algorithm. ACM Transactions on Mathematical Software, 13(3):262-280.

Kirkpatrick, S., Gelatt, C. D. and Vecchi, M. P. (1983). Optimization by Simulated Annealing. Science, 220(4598):671-680.

Pronzato, L., Walter, E., Venot, A. and Lebruchec, J.-F. (1984). A general-purpose global optimizer: Implimentation and applications. Mathematics and Computers in Simulation, 26(5):412-422.
}
\author{
Kai Husmann
}

\seealso{
\code{\link{optim_nm}}, \code{\link{optim}}, \code{\link{plot.optim_nmsa}}
}

\examples{
##### Rosenbrock function
# minimum at f(1,1) = 0
ro <- function(x){
  100*(x[2]-x[1]^2)^2+(1-x[1])^2
}

# Random start values. Example arguments for the relatively simple Rosenbrock function.
ro_sa <- optim_sa(fun = ro,
                  start = c(runif(2, min = -1, max = 1)),
                  lower = c(-5, -5),
                  upper = c(5, 5),
                  trace = TRUE,
                  control = list(t0 = 100,
                            nlimit = 550,
                            t_min = 0.1,
                            dyn_rf = FALSE,
                            rf = 1,
                            r = 0.7
                  )
         )


# Visual inspection.
plot(ro_sa)
plot(ro_sa, type = "contour")


##### Holder table function

# 4 minima at
  #f(8.055, 9.665) = -19.2085
  #f(-8.055, 9.665) = -19.2085
  #f(8.055, -9.665) = -19.2085
  #f(-8.055, -9.665) = -19.2085

ho <- function(x){
  x1 <- x[1]
  x2 <- x[2]

  fact1 <- sin(x1) * cos(x2)
  fact2 <- exp(abs(1 - sqrt(x1^2 + x2^2) / pi))
  y <- -abs(fact1 * fact2)
}

# Random start values. Example arguments for the relatively complex Holder table function.
optim_sa(fun = ho,
         start = c(1, 1),
         lower = c(-10, -10),
         upper = c(10, 10),
         trace = TRUE,
         control = list(dyn_rf = FALSE,
                        rf = 1.6,
                        t0 = 10,
                        nlimit = 200,
                        r = 0.6,
                        t_min = 0.1
         )
)

#### Himmelblau's function
# 4 minima at
  # f(3, 2) = 0
  # f(-2.804, -3.131) = 0
  # f(-3.779, -3.283) = 0
  # f( 3.584, -1.848) = 0

hi <- function(x){
  (x[1]**2 + x[2] - 11)**2 + (x[1] + x[2]**2 -7)**2
}

# Random start values. Example arguments for integer programming.
# Only the integer solution will be found.

var_func_int <- function(para_0, fun_length, rf, temp = NA){
  ret_var_func <- para_0 + sample.int(rf, fun_length, replace = TRUE) *
  ((rbinom(fun_length, 1, 0.5) * -2) + 1)
  return (ret_var_func)
}

optim_sa(fun = hi,
         start = round(c(runif(2, min = -1, max = 1))),
                trace = TRUE,
                lower = c(-4, -4),
                upper=c(4, 4),
                control = list(t0 = 1000,
                               nlimit = 1500,
                               r = 0.8,
                               vf = var_func_int,
                               rf = 3
                )
)
}
