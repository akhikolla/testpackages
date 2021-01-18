\name{plot}
\alias{plot}
\alias{plot.rgasp}
\alias{plot,rgasp-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Plot for Robust GaSP model
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Function to make plots on Robust GaSP models after the Robust GaSP model has been constructed.
}
\usage{
\S4method{plot}{rgasp}(x,y, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ an object of  class \code{rgasp}.}
  \item{y}{ not used.  }
  \item{...}{ additional arguments not implemented yet.}

}



%\details{
%%  ~~ If necessary, more details than the description above ~~
%Provide here some details.
%}


\value{

Three plots: the leave-one-out fitted values versus exact values, standardized residuals and QQ plot.
}
\references{
M. Gu. (2016). Robust uncertainty quantification and scalable computation for computer models with massive output. Ph.D. thesis. Duke University.
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

\examples{
 library(RobustGaSP)
  #------------------------
  # a 3 dimensional example
  #------------------------
  # dimensional of the inputs
  dim_inputs <- 3    
  # number of the inputs
  num_obs <- 30       
  # uniform samples of design
  input <- matrix(runif(num_obs*dim_inputs), num_obs,dim_inputs) 
  
  # Following codes use maximin Latin Hypercube Design, which is typically better than uniform
  # library(lhs)
  # input <- maximinLHS(n=num_obs, k=dim_inputs)  ##maximin lhd sample
  
  # outputs from the 3 dim dettepepel.3.data function
  
  output = matrix(0,num_obs,1)
  for(i in 1:num_obs){
    output[i]<-dettepepel.3.data (input[i,])
  }
  
  # use constant mean basis, with no constraint on optimization
  m1<- rgasp(design = input, response = output, lower_bound=FALSE)
 
  # plot
  plot(m1)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
