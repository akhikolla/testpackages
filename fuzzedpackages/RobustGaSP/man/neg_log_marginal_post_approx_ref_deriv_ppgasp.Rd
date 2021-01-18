\name{neg_log_marginal_post_approx_ref_deriv_ppgasp}
\alias{neg_log_marginal_post_approx_ref_deriv_ppgasp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Derivative of negative natural logarithm of approximate reference marginal posterior density of the PP GaSP model
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
The function computes the derivative (with regard to log of inverse range parameter) of natural logarithm of marginal posterior density of the PP GaSP model with the jointly robust  prior after marginalizing out the mean (trend) and variance parameters  by the location-scale prior.
}
\usage{
neg_log_marginal_post_approx_ref_deriv_ppgasp(param, nugget, nugget.est, 
                                       R0, X, zero_mean,output, CL, a, b, 
                                      kernel_type, alpha)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{param}{
%%     ~~Describe \code{param} here~~
A vector of natural logarithm of inverse-range parameters and natural logarithm of the nugget-variance ratio parameter.
}
  \item{nugget}{
%%     ~~Describe \code{nugget} here~~
The nugget-variance ratio parameter if this parameter is fixed.
}
  \item{nugget.est}{
%%     ~~Describe \code{nugget.est} here~~
Boolean value of whether the nugget is estimated or fixed. 
}
  \item{R0}{
%%     ~~Describe \code{R0} here~~
A List of matrix where the j-th matrix is an absolute difference matrix of the j-th input vector.
}
  \item{X}{
%%     ~~Describe \code{X} here~~
The mean basis function i.e. the trend function.
}
  \item{zero_mean}{
%%     ~~Describe \code{zero_mean} here~~
The mean basis function is zero or not.
}

  \item{output}{
%%     ~~Describe \code{output} here~~
The output matrix.
}
  \item{CL}{
%%     ~~Describe \code{CL} here~~
Pseudoparameter in the approximate reference prior.
}
  \item{a}{
%%     ~~Describe \code{a} here~~
Pseudoparameter in the approximate reference prior.
}
  \item{b}{
%%     ~~Describe \code{b} here~~
Pseudoparameter in the approximate reference prior.
}
  \item{kernel_type}{
%%     ~~Describe \code{kernel_type} here~~
A vector of \code{integer} specifying the type of kernels of each coordinate of the input. 
In each coordinate of the vector, 1 means the \code{pow_exp} kernel with roughness parameter specified in alpha; 2 means  \code{matern_3_2}  kernel; 3 means \code{matern_5_2}  kernel; 5 means \code{periodic_gauss} kernel; 5 means \code{periodic_exp} kernel. 

%Type of kernel. \code{matern_3_2} and \code{matern_5_2} are \code{Matern kernel} with roughness parameter 3/2 and 5/2 respectively. \code{pow_exp} is power exponential kernel with roughness parameter alpha. If \code{pow_exp} is to be used, one needs to specify its roughness parameter alpha.
}
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
Roughness parameters in the kernel functions.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
The derivative of natural logarithm of marginal posterior density with the jointly robust prior.
}
\references{
M. Gu. and J.O. Berger (2016). Parallel partial Gaussian process emulation for computer models with massive output. \emph{Annals of Applied Statistics}, 10(3), 1317-1347.

M. Gu. (2016). Robust uncertainty quantification and scalable computation for computer models with massive output. Ph.D. thesis. Duke University.

M. Gu (2018), \emph{Jointly robust prior for Gaussian stochastic process in emulation, calibration and variable selection}, arXiv:1804.09329.

}
\author{
%%  ~~who you are~~
\packageAuthor{RobustGaSP}

Maintainer: \packageMaintainer{RobustGaSP}
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{neg_log_marginal_post_approx_ref_ppgasp}}.
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{internal} % use one of  RShowDoc("KEYWORDS")
