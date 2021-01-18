\name{higdon.1.data}
\alias{higdon.1.data}
\alias{higdon}
\alias{hig02}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Higdon (2002) Function
}
\description{

Higdon (2002) 1-dimensional test function.  
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
%This function is a simple example used to test \code{s} method of space modeling. The values of the input variable \code{s} are spatial locations. A small error term \ifelse{html}{\out{epsilon ~ N(0, 0.12)}}{\eqn{\varepsilon \sim N(0, 0.12)}{epsilon ~ N(0, 0.12)}} is added to the response.
}
\usage{
higdon.1.data(s)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{s}{ a 1-dimensional vector specifying the location where the function is to be evaluated.  
%%     ~~Describe \code{s} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
For more details, see Higdon D. (2002) Space and Space-Time Modeling using Process Convolutions. In: Anderson C.W., Barnett V., Chatwin P.C., El-Shaarawi A.H. (eds) Quantitative Methods for Current Environmental Issues. Springer, London. DOI:\url{https://doi.org/10.1007/978-1-4471-0657-9_2}.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...

 A real number equal to the Higdon (2002) function value at \code{s}.
 
%\ifelse{html}{\out{f(s)}}{\eqn{f(s)=\sin(2\pi s/10)+0.2\sin(2\pi s/2.5)}{f(s)}}
}
\references{
%% ~put references to the literature/web site here ~

Higdon, D. (2002). Space and space-time modeling using process convolutions. In Quantitative methods for current environmental issues (pp. 37-56). Springer London.
\url{https://doi.org/10.1007/978-1-4471-0657-9_2}

S. Surjanovic, D. Bingham, Virtual Library of Simulation Experiments:  Test Functions
and Datasets, retrieved March 29, 2016, from \url{https://www.sfu.ca/~ssurjano/hig02.html}.

}
\author{
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
\examples{

s <- seq(0,10,0.01)
y <- higdon.1.data(s)
plot(s,y, xlab='s',ylab='y',type='l',main='Higdon (2002) function')

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{internal}