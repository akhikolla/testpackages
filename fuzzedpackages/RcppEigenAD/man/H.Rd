

\name{H}
\alias{H}
\title{Construct a function to calculate the Hessian of a function.}
\usage{
H(f)
}
\arguments{
\item{f}{A function created using either \link{sourceCppAD} or  the
  composition operator \link{\%.\%}.}
}
\value{
A function which computes the Hessian of the function.
}
\description{

Constructs a function to calculate the Hessian of a function produced either using \link{sourceCppAD} or the composition operator \link{\%.\%} .
The returned function has the same argument signature as f but returns a matrix representing the
blocked Hessians of \code{f} evaluated at the functions arguments. The partial
derivatives are formed with respect to the argument specified when f was created with sourceCppAD.

In what follows it is assumed that f has a single matrix argument (the
one with which f is differentiated with respect to). When this is not
the case, the other arguments will be considered constant at the point the
Hessian is evaluated at. Consequently, the structure of the output of
the function produced by \code{H} is unchanged by the additional arguments.
The blocked Hessian \eqn{{\bf H}} is organised as follows.

If \eqn{f:{\bf R}^{n} \rightarrow {\bf R}^{m}} where \eqn{{\bf Y}_{n_{Y}
    \times m_{Y}} = f({\bf X}_{n_{X} \times m_{X}})} and
\eqn{n=n_{X}m_{X}}, \eqn{m=n_{Y}m_{Y}} then by numbering the elements of
    the matrices row-wise so that,
\deqn{
  {\bf Y} =
  \left[
    \begin{array}{ccc}
      y_{1} & \dots & y_{m_{Y}} \\
      y_{m_{Y}+1} & \dots & y_{2m_{Y}} \\  
      \vdots & \ddots & \vdots \\
       y_{(n_{Y}-1)m_{Y}+1} & \dots & y_{n_{Y}m_{Y}}
    \end{array}
    \right]
}
and
\deqn{
  {\bf X} =
  \left[
    \begin{array}{ccc}
      x_{1} & \dots & x_{m_{X}} \\
      x_{m_{X}+1} & \dots & x_{2m_{X}} \\  
      \vdots & \ddots & \vdots \\
       x_{(n_{X}-1)m_{X}+1} & \dots & x_{n_{X}m_{X}}
    \end{array}
    \right]
}
then the \eqn{n \times nm} blocked Hessian matrix \eqn{{\bf H}} is structured as 
\deqn{
  {\bf H} =
  \left[
    \begin{array}{cccc}
        {\bf H}_{1} & {\bf H}_{2} & \dots & {\bf H}_{m} 
      \end{array}
    \right]
}
with \eqn{{\bf H}_{k}} the Hessian matrix for \eqn{y_{k}} i.e.,
\deqn{
\begingroup
\renewcommand*{\arraystretch}{1.5}
 {\bf H}_{k} = 
  \left[ \begin{array}{ccc}
      \frac{\partial^2 y_{k}}{\partial x_1 \partial x_1} & \dots & \frac{\partial^2 y_{k}}{\partial x_1 \partial x_n} \\
      \vdots & \ddots & \vdots \\
      \frac{\partial^2 y_{k}}{\partial x_n \partial x_1} & \dots  & \frac{\partial^2 y_{k}}{\partial x_n \partial x_n} 
      \end{array} \right]
\endgroup
}
}
\examples{
\donttest{
library(RcppEigenAD)
# define f as the eigen vectors of its argument x
# calculated using the Eigen library
f<-sourceCppAD('
ADmat f(const ADmat& X)
{
   Eigen::EigenSolver<ADmat> es(X);
   return es.pseudoEigenvectors();
}
')
Hf<-H(f)
X<-matrix(c(1,2,3,4),2,2)
Hmat<-Hf(X) 
Hmat # the Hessian matrices of second derivatives stacked column wise
Hmat[3,9] # the second derivative of f(X)[2,1] with respect to X[2,1] and X[1,2]
}
}