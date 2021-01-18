
\name{J}
\alias{J}
\title{Construct a function to calculate the Jacobian of a function.}
\usage{
J(f)
}
\arguments{
\item{f}{A function created using either \link{sourceCppAD} or the
  composition operator \link{\%.\%}.}
}
\value{
A function which computes the Jacobian of the function
}
\description{
Constructs a function to calculate the Jacobian of a function produced
either using \link{sourceCppAD} or the composition operator \link{\%.\%}.
The returned function has the same argument signature as f but returns a matrix representing the
Jacobian of \code{f} evaluated at the functions arguments. The partial
derivatives are formed with respect to the arguments specified when f was created with sourceCppAD.


In what follows it is assumed that f has a single matrix argument (the
one with which f is differentiated with respect to). When this is not
the case, the other arguments will be considered constant at the point the
Jacobian is evaluated at. Consequently, the structure of the output of
the function produced by \code{J} is unchanged by the additional arguments.
The Jacobian matrix \eqn{{\bf J}} is organised as follows.

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
then the \eqn{m \times n} Jacobian matrix is given by 
\deqn{
\begingroup
\renewcommand*{\arraystretch}{1.5}
 {\bf J} = 
  \left[ \begin{array}{ccc}
      \frac{\partial y_{1}}{\partial x_1} & \dots & \frac{\partial y_{1}}{\partial x_n} \\
      \vdots & \ddots & \vdots \\
      \frac{\partial y_{m}}{\partial x_1} & \dots  & \frac{\partial y_{m}}{\partial x_n} 
      \end{array} \right]
\endgroup
}
}

\examples{
\donttest{
# define f as the eigen vectors of its argument X
# calculated using the Eigen library
library(RcppEigenAD)
f<-sourceCppAD('
ADmat f(const ADmat& X)
{
   Eigen::EigenSolver<ADmat> es(X);
   return es.pseudoEigenvectors();
}
')
Jf<-J(f)
X<-matrix(c(1,2,3,4),2,2)
Jmat<-Jf(X)
Jmat # the Jacobian matrix of first derivatives
Jmat[2,3] # the derivative of f(X)[1,2] with respect to X[2,1]
}
}