\name{separable_kernel}
\alias{separable_kernel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Product correlation matrix with the product form
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Function to construct the product correlation matrix with the product form.
}
\usage{
separable_kernel(R0, beta, kernel_type, alpha)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{R0}{
%%     ~~Describe \code{R0} here~~
A List of matrix where the j-th matrix is an absolute difference matrix of the j-th input vector.
}
  \item{beta}{
%%     ~~Describe \code{beta} here~~
The range parameters.
}
  \item{kernel_type}{
%%     ~~Describe \code{kernel_type} here~~
A vector specifying the type of kernels of each coordinate of the input. \code{matern_3_2} and \code{matern_5_2} are \code{Matern correlation} with roughness parameter 3/2 and 5/2 respectively. \code{pow_exp} is power exponential correlation with roughness parameter alpha. If \code{pow_exp} is to be used, one needs to specify its roughness parameter alpha. The default choice is \code{matern_5_2}. The \code{periodic_gauss} means the Gaussian kernel with periodic folding method with be used. The code{periodic_exp} means the exponential kernel with periodic folding method will be used.
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
The product correlation matrix with the product form.
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

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{internal}