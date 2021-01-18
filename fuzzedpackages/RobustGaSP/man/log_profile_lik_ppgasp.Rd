\name{log_profile_lik_ppgasp}
\alias{log_profile_lik_ppgasp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Natural logarithm of profile likelihood of the PP GaSP model
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function computes the natural logarithm of profile likelihood of the PP GaSP model after plugging in the maximum likelihood estimator of the mean (trend) and variance parameters.
}
\usage{
log_profile_lik_ppgasp(param, nugget, nugget_est, R0, X, zero_mean
,output,kernel_type, alpha)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{param}{
%%     ~~Describe \code{param} here~~
a vector of natural logarithm of inverse-range parameters and natural logarithm of the nugget-variance ratio parameter.
}
  \item{nugget}{
%%     ~~Describe \code{nugget} here~~
the nugget-variance ratio parameter if this parameter is fixed.
}
  \item{nugget_est}{
%%     ~~Describe \code{nugget.est} here~~
Boolean value of whether the nugget is estimated or fixed.
}
  \item{R0}{
%%     ~~Describe \code{R0} here~~
a List of matrix where the j-th matrix is an absolute difference matrix of the j-th input vector.
}
  \item{X}{
%%     ~~Describe \code{X} here~~
the mean basis function i.e. the trend function.
}
  \item{zero_mean}{
%%     ~~Describe \code{zero_mean} here~~
the mean basis function is zero or not.
}

  \item{output}{
%%     ~~Describe \code{output} here~~
a matrix where each row is one runs of the computer model output.
}
  \item{kernel_type}{
%%     ~~Describe \code{kernel_type} here~~
A vector of \code{integer} specifying the type of kernels of each coordinate of the input. 
In each coordinate of the vector, 1 means the \code{pow_exp} kernel with roughness parameter specified in alpha; 2 means  \code{matern_3_2}  kernel; 3 means \code{matern_5_2}  kernel; 5 means \code{periodic_gauss} kernel; 5 means \code{periodic_exp} kernel. 

%A vector specifying the type of kernels of each coordinate of the input. \code{matern_3_2} and \code{matern_5_2} are \code{Matern correlation} with roughness parameter 3/2 and 5/2 respectively. \code{pow_exp} is power exponential correlation with roughness parameter alpha. If \code{pow_exp} is to be used, one needs to specify its roughness parameter alpha. The default choice is \code{matern_5_2}. The \code{periodic_gauss} means the Gaussian kernel with periodic folding method with be used. The \code{periodic_exp} means the exponential kernel with periodic folding method will be used.
}
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
roughness parameters in the kernel functions.
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
The numerical value of natural logarithm of the profile likelihood.
}
\references{
M. Gu. and J.O. Berger (2016). Parallel partial Gaussian process emulation for computer models with massive output. \emph{Annals of Applied Statistics}, 10(3), 1317-1347.

M. Gu. (2016). Robust Uncertainty Quantification and Scalable Computation for Computer Models with Massive Output. Ph.D. thesis. Duke University.
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

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")

\keyword{internal}