\name{construct_ppgasp}
\alias{construct_ppgasp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{PP GaSP constructor after estimating the parameters}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function constructs the PP GaSP model if the range and noise-variance ratio parameters are given or have been estimated.
}
\usage{
construct_ppgasp(beta, nu, R0, X, zero_mean, output, kernel_type, alpha)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{beta}{
%%     ~~Describe \code{beta} here~~
inverse-range parameters.
}
  \item{nu}{
%%     ~~Describe \code{nu} here~~
noise-variance ratio parameter.
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
the output matrix.
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
A list. The first element is a lower triangular matrix \code{L}, a cholesky decomposition of \code{R}, i.e. \ifelse{html}{\out{LL^t=R}}{\eqn{LL^t=R}{LL^t=R}}
and \code{R} is the correlation matrix. The second element is  lower triangular matrix \code{LX} a cholesky decomposition of \ifelse{html}{\out{X^tR^{-1}X}}{\eqn{X^tR^{-1}X}{X^tR^{-1}X}}. The third element is a matrix of \code{theta_hat}, the mean (trend) parameters. The last element is a vector of \code{sigma2_hat}, which is the estimated variance parameter on each function.  
}
%\references{
%% ~put references to the literature/web site here ~
%}
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

\keyword{internal}% __ONLY ONE__ keyword per line

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
