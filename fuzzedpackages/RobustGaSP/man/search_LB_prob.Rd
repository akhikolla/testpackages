\name{search_LB_prob}
\alias{search_LB_prob}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Search for the default lower bound of range parameters. 
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Function to find the values to construct the default lower bound of range parameters.
}
\usage{
search_LB_prob(param, R0, COND_NUM_UB, p, kernel_type, alpha, nugget)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{param}{
%%     ~~Describe \code{param} here~~
A vector of natural logarithm of inverse-range parameters and natural logarithm of the nugget-variance ratio parameter.
}
  \item{R0}{
%%     ~~Describe \code{R0} here~~
A List of matrix where the j-th matrix is an absolute difference matrix of the j-th input vector.
}
  \item{COND_NUM_UB}{
%%     ~~Describe \code{COND_NUM_UB} here~~
The maximum condition number of the correlation matrix.
}
  \item{p}{
%%     ~~Describe \code{p} here~~
}
  \item{kernel_type}{
%%     ~~Describe \code{kernel_type} here~~
Type of kernel. \code{matern_3_2} and \code{matern_5_2} are \code{Matern kernel} with roughness parameter 3/2 and 5/2 respectively. \code{pow_exp} is power exponential kernel with roughness parameter alpha. If \code{pow_exp} is to be used, one needs to specify its roughness parameter alpha.
}
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
Roughness parameters in the kernel functions.
}
  \item{nugget}{
%%     ~~Describe \code{nugget} here~~
The nugget-variance ratio parameter if this parameter is fixed.
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
A \code{vector} of values used in constructing the default lower bound of range parameters.
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

\keyword{internal}