\name{k1Matern3_2}
\alias{k1FunExp}
\alias{k1FunPowExp}
\alias{k1FunGauss}
\alias{k1FunMatern3_2}
\alias{k1FunMatern5_2}
%%
\alias{k1Fun1Cos}
\alias{k1Fun1Exp}
\alias{k1Fun1PowExp}
\alias{k1Fun1Gauss}
\alias{k1Fun1Matern3_2}
\alias{k1Fun1Matern5_2}
%%
\title{
  
  One-Dimensional Classical Covariance Kernel Functions

}
\description{

  One-dimensional classical covariance kernel Functions.

}
\usage{
k1FunExp(x1, x2, par)
k1FunGauss(x1, x2, par)
k1FunPowExp(x1, x2, par)
k1FunMatern3_2(x1, x2, par)
k1FunMatern5_2(x1, x2, par)

k1Fun1Cos(x)
k1Fun1Exp(x)
k1Fun1Gauss(x)
k1Fun1PowExp(x, alpha = 1.5)
k1Fun1Matern3_2(x)
k1Fun1Matern5_2(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x1}{
    
    First location vector.

  }
  \item{x2}{

    Second location vector. Must have the same length as \code{x1}. 

  }
  
  \item{x}{
    
    For stationary covariance functions, the vector containing difference
    of positions: \code{x = x1 - x2}.
    
  }  
  \item{alpha}{
    
    Regularity parameter in \eqn{(0, 2]} for Power Exponential
    covariance function.

  }
  \item{par}{
    
    Vector of parameters. The length and the meaning of the elements in
    this vector depend on the chosen kernel. The first parameter is the
    range parameter (if there is one), the last is the variance. So the
    shape parameter of \code{k1FunPowExp} is the second one out of the
    three parameters.

  }
}
\details{

  These kernel functions are described in the Roustant et al (2012),
  table 1 p. 8. More details are given in chap. 4 of Rasmussen et al
  (2006).

}

\value{

  A matrix with a \code{"gradient"} attribute. This matrix has \eqn{n_1}
  rows and \eqn{n_2} columns where \eqn{n_1} and \eqn{n_2} are the
  length of \code{x1} and \code{x2}.  If \code{x1} and \code{x2} have
  length 1, the attribute is a vector of the same length \eqn{p} as
  \code{par} and gives the derivative of the kernel with respect to the
  parameters in the same order. If \code{x1} or \code{x2} have length
  \eqn{> 1}, the attribute is an array with dimension \eqn{(n_1, n_2,
  p)}.

}
\references{

  C.E. Rasmussen and C.K.I. Williams (2006), \emph{Gaussian Processes
    for Machine Learning}, the MIT Press,
  \url{http://www.GaussianProcess.org/gpml}
  
  O. Roustant, D. Ginsbourger, Y. Deville (2012).
  "DiceKriging, DiceOptim: Two R Packages for the Analysis of
  Computer Experiments by Kriging-Based Metamodeling and Optimization."
  \emph{Journal of Statistical Software}, 51(1), 1-55.
  \url{http://www.jstatsoft.org/v51/i01/}
  %%M.L. Stein (1999), \emph{Interpolation of spatial data, some theory
  %%  for kriging}, Springer.

}
\author{

  Oivier Roustant, David Ginsbourger, Yves Deville
  
}
\note{

  The kernel functions are coded in C through the \code{.Call} interface
  and are mainly intended for internal use. They are used by the
  \code{covTS} class.
  
  Be aware that very few checks are done (length of objects, order of
  the parameters, ...).

}
%%\seealso{}
\examples{
## show the functions
n <- 300
x0 <- 0
x <- seq(from = 0, to = 3, length.out = n)
kExpVal <- k1FunExp(x0, x, par = c(range = 1, var = 2))
kGaussVal <- k1FunGauss(x0, x, par = c(range = 1, var = 2))
kPowExpVal <- k1FunPowExp(x0, x, par = c(range = 1, shape = 1.5, var = 2))
kMatern3_2Val <- k1FunMatern3_2(x0, x, par = c(range = 1, var = 2))
kMatern5_2Val <- k1FunMatern5_2(x0, x, par = c(range = 1, var = 2))
kerns <- cbind(as.vector(kExpVal), as.vector(kGaussVal), as.vector(kPowExpVal),
               as.vector(kMatern3_2Val), as.vector(kMatern5_2Val))
matplot(x, kerns, type = "l", main = "five 'kergp' 1d-kernels", lwd = 2)

## extract gradient
head(attr(kPowExpVal, "gradient"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%%\keyword{ ~kwd1 }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
