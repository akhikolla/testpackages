\name{neg_log_marginal_post_approx_ref_deriv}
\alias{neg_log_marginal_post_approx_ref_deriv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Derivative of negative natural logarithm of approximate reference marginal posterior density
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
The function computes the derivative (with regard to log of inverse range parameter) of natural logarithm of marginal posterior density with the jointly robust prior prior after marginalizing out the mean (trend) and variance parameters  by the location-scale prior.
}
\usage{
neg_log_marginal_post_approx_ref_deriv(param, nugget, nugget.est, 
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
The output vector.
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
Type of kernel. \code{matern_3_2} and \code{matern_5_2} are \code{Matern kernel} with roughness parameter 3/2 and 5/2 respectively. \code{pow_exp} is power exponential kernel with roughness parameter alpha. If \code{pow_exp} is to be used, one needs to specify its roughness parameter alpha.
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
The derivative of natural logarithm of marginal posterior density with jointly robust prior prior.
}
\references{
Mengyang Gu. (2016). Robust Uncertainty Quantification and Scalable Computation for Computer Models with Massive Output. Ph.D. thesis. Duke University.

M. Gu (2018), \emph{Jointly Robust Prior for Gaussian Stochastic Process in Emulation, Calibration and Variable Selection}, arXiv:1804.09329.

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
\code{\link{neg_log_marginal_post_approx_ref}}.
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{internal} % use one of  RShowDoc("KEYWORDS")
