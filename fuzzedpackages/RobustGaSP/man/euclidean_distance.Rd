\name{euclidean_distance}
\alias{euclidean_distance}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Euclidean distance matrix between two input matrices
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Function to construct the euclidean distance matrix with the two input matrices.
}
\usage{
euclidean_distance(input1, input2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input1}{
%%     ~~Describe \code{R0} here~~
A matrix of input with the number of rows being the number of observations and the number of columns being the number of variables
}
  \item{input2}{
%%     ~~Describe \code{beta} here~~
A matrix of input with the number of rows being the number of observations and the number of columns being the number of variables
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
The  euclidean distance matrix with the number of rows and the number of columns being the number of rows in the first and second input matrices.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
%%  ~~who you are~~
\packageAuthor{RobustGaSP}

Maintainer: \packageMaintainer{RobustGaSP}
}


\examples{
  # dimensional of the inputs
  dim_inputs <- 8    
  
  # number of the inputs
  num_obs <- 30       
  
  # uniform samples of design
  input <-matrix(runif(num_obs*dim_inputs), num_obs,dim_inputs) 

  # the Euclidean distance matrix 
  R0=euclidean_distance(input, input)
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