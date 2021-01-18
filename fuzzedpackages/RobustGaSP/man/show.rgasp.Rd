\name{show}
\alias{show}
\alias{show.rgasp}
\alias{show.rgasp-class}
\alias{show,rgasp-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Show Robust GaSP object
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Function to print Robust GaSP models after the Robust GaSP model has been constructed.
}
\usage{
\S4method{show}{rgasp}(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ an object of  class \code{rgasp}.}
}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%Provide here some details.
%}

%\references{
%% ~put references to the literature/web site here ~
%}

\author{
%%  ~~who you are~~
\packageAuthor{RobustGaSP}

Maintainer: \packageMaintainer{RobustGaSP}
}

\examples{
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
  
  ####
  # outputs from the 3 dim dettepepel.3.data function
  
  output = matrix(0,num_obs,1)
  for(i in 1:num_obs){
    output[i]<-dettepepel.3.data (input[i,])
  }
  
  # use constant mean basis, with no constraint on optimization
  m1<- rgasp(design = input, response = output, lower_bound=FALSE)
  
  # the following use constraints on optimization
  # m1<- rgasp(design = input, response = output, lower_bound=TRUE)
  
  # the following use a single start on optimization
  # m1<- rgasp(design = input, response = output, lower_bound=FALSE)
  
  
  show(m1)
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
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
