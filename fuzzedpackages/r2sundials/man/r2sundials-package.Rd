\name{r2sundials-package}
\alias{r2sundials-package}
\alias{r2sundials}
\docType{package}
\title{
\packageTitle{r2sundials}
}
\description{
\packageDescription{r2sundials}
Interface for commonly used SUNDIALS library for solving ODE system. It is made R friendly by allowing to pass custom parameters to user's callback functions. They can be both written in R and in C++ (Rcpp flavor). In case of C++, performance is greatly improved so this option is highly advisable when performance matters. Jacobian matrix can be specified as either dense or sparse matrix. In the latter case \pkg{rmumps} package is used to solve corresponding linear systems. Root finding and corresponding event management are optional and can be too specified as R or C++ functions which makes them a very flexible tool for controlling the ODE system during the time course simulation.
}
\details{
The DESCRIPTION file:
\packageDESCRIPTION{r2sundials}
\packageIndices{r2sundials}
The main function of the package is r2cvodes() which wraps and converts its arguments in data structures and parameters convenient for calling \code{cvodes()} from SUNDIALS library.

When using r2sundials::r2cvodes(), some callback functions have to be written by the user either in R or in C++ (RcppArmadillo flavor). One of them is mandatory and defines the rhs of the ODE system to solve. Others are optional and can be used to
\itemize{
 \item calculate Jacobian matrix in sparse or dense format;
 \item calculate root vector for tracking particular events and
 \item handle them in a way defined by the user himself;
 \item calculate sensitivity rhs if sensitivity to some parameters is required and user does not want to rely on internal procedures for estimating such rhs.
}

Note that if 'SUNDIALS' 'CVODES' is compiled with parameter SUNDIALS_INDEX_SIZE set to 32, some memory copying is spared if a C++ function calculating sparse Jacobian is provided by user.

The version number of this package if made of CVODES original version, e.g. 5.0.0 followed by one digit for R wrapper release. This can give 5.0.0-1.
}
\author{
\packageAuthor{r2sundials}
Maintainer: \packageMaintainer{r2sundials}
}
\references{
\href{https://computation.llnl.gov/sites/default/files/public/cvs_guide.pdf}{Original SUNDIALS CVODES user documentation}
}
\keyword{ ODE }
\keyword{ sensitivity }
\seealso{
 \code{\link[deSolve:deSolve-package]{deSolve}}
}
\examples{
  # a very simple ODE for exponential transition from 0 to 1: y'=1-y, y(0)=0
  frhs=function(t, y, param, psens) 1-y
  ti=seq(0, 5, length.out=101)
  y0=0
  res_exp=r2sundials::r2cvodes(y0, ti, frhs)
  plot(ti, res_exp[1,], t="l", xlab="Time", ylab="Y")
}
