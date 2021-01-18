\name{separable_multi_kernel}
\alias{separable_multi_kernel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Product correlation matrix with the product form
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Function to construct the product correlation matrix with the product form. The kernel can be different for each coordinate of the input. 
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
A vector of \code{integer} specifying the type of kernels of each coordinate of the input. 
In each coordinate of the vector, 1 means the \code{pow_exp} kernel with roughness parameter specified in alpha; 2 means  \code{matern_3_2}  kernel; 3 means \code{matern_5_2}  kernel; 5 means \code{periodic_gauss} kernel; 5 means \code{periodic_exp} kernel. 
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