\name{Borehole}
\alias{borehole}
\alias{Borehole}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Borehole function
%%  ~~function to do ... ~~
}
\description{

It is an 8-dimensional test function that models water flow through a borehole. Its simplicity and quick evaluation makes it a commonly used function for testing a wide variety of methods in computer experiments. 

}
\usage{
borehole(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ an 8-dimensional vector specifying the location where the function is to be evaluated. }
%%     ~~Describe \code{x} here~~
}

\details{
%%  ~~ If necessary, more details than the description above ~~
For more details, see Worley, B. A. (1987). Deterministic uncertainty analysis (No. CONF-871101-30). Oak Ridge National Lab., TN (USA).
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
The response is given by
\deqn{
f(x) = \frac{2\pi T_u (H_u - H_l)}{\ln(\frac{r}{r_w})\Big(1+\frac{2LT_u}{\ln(\frac{r}{r_w})r^2_w K_w}+\frac{T_u}{T_l}\Big)}}{Borehole Function (r_w,r,T_u,H_u,T_l,H_l,L,K_w)}

where \code{x} is an 8-dimensional vector with \code{rw <- x[1]}, \code{r  <- x[2]}, \code{Tu <- x[3]}, \code{Hu <- x[4]}, \code{Tl <- x[5]}, \code{Hl <- x[6]}, \code{L  <- x[7]}, \code{Kw <- x[8]}.

}
\references{
%% ~put references to the literature/web site here ~

Harper, W. V., & Gupta, S. K. (1983). Sensitivity/uncertainty analysis of a borehole scenario comparing Latin Hypercube Sampling and deterministic sensitivity approaches (No. BMI/ONWI-516). Battelle Memorial Inst., Columbus, OH (USA). Office of Nuclear Waste Isolation.

Worley, B. A. (1987). Deterministic uncertainty analysis (No. CONF-871101-30). Oak Ridge National Lab., TN (USA).

S. Surjanovic, D. Bingham, Virtual Library of Simulation Experiments:  Test Functions
and Datasets, retrieved March 29, 2016, from \url{http://www.sfu.ca/~ssurjano/borehole.html}.

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

%\examples{
%%##---- Should be DIRECTLY executable !! ----
%%##-- ==>  Define data, use random,
%%##--	or do  help(data=index)  for the standard data sets.
%%
%%## The function is currently defined as
%%function (x) 
%%{
%%  }
%}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{optimize}
\keyword{internal}% __ONLY ONE__ keyword per line
