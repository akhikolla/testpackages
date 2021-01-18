\name{neg_log_marginal_post_approx_ref}
\alias{neg_log_marginal_post_approx_ref}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Natural logarithm of approximate reference marginal posterior density of the robust GaSP model
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Natural logarithm of marginal posterior density with the approximate reference prior of inverse range parameter after marginalizing out the mean (trend) and variance parameters  by the location-scale prior.
}
\usage{
neg_log_marginal_post_approx_ref(param, nugget, nugget.est, 
                                  R0, X, zero_mean,output, CL, a, b, 
                                  kernel_type, alpha)
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
  \item{nugget.est}{
%%     ~~Describe \code{nugget.est} here~~
Boolean value of whether the nugget is estimated or fixed.
}
  \item{R0}{
%%     ~~Describe \code{R0} here~~
a list of matrix where the j-th matrix is an absolute difference matrix of the j-th input vector.
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
the output vector.
}
  \item{CL}{
%%     ~~Describe \code{CL} here~~
prior parameters in the approximate reference prior.
}
  \item{a}{
%%     ~~Describe \code{a} here~~
prior parameter in the approximate reference prior.
}
  \item{b}{
%%     ~~Describe \code{b} here~~
prior parameter in the approximate reference prior.
}
  \item{kernel_type}{
%%     ~~Describe \code{kernel_type} here~~
type of kernel. \code{matern_3_2} and \code{matern_5_2} are \code{Matern kernel} with roughness parameter 3/2 and 5/2 respectively. \code{pow_exp} is power exponential kernel with roughness parameter alpha. If \code{pow_exp} is to be used, one needs to specify its roughness parameter alpha.
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
The natural logarithm of the marginal posterior density with approximate reference prior of inverse range parameter (beta parameterization). 
}

\references{
Mengyang Gu. (2016). Robust Uncertainty Quantification and Scalable Computation for Computer Models with Massive Output. Ph.D. thesis. Duke University.
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
\keyword{internal}% use one of  RShowDoc("KEYWORDS")
