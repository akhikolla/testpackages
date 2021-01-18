\name{leave_one_out_rgasp}
\alias{leave_one_out_rgasp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
leave-one-out fitted values and standard deviation for robust GaSP model
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
A function to calculate leave-one-out fitted values and the standard deviation of the prediction on robust GaSP models after the robust GaSP model has been constructed.
}
\usage{
leave_one_out_rgasp(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
%%     ~~Describe \code{beta} here~~
  \item{object}{ an object of  class \code{rgasp}.}
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
A list of 2 elements with 
\item{mean }{leave one out fitted values.}
\item{sd }{standard deviation of each prediction.}


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

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{rgasp}}
}



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
  
  ####
  # outputs from the 3 dim dettepepel.3.data function
  
  output = matrix(0,num_obs,1)
  for(i in 1:num_obs){
    output[i]<-dettepepel.3.data (input[i,])
  }
  
  # use constant mean basis, with no constraint on optimization
  m1<- rgasp(design = input, response = output, lower_bound=FALSE)
  
  ##leave one out predict
  leave_one_out_m1=leave_one_out_rgasp(m1)
  
  ##predictive mean 
  leave_one_out_m1$mean
  ##standard deviation
  leave_one_out_m1$sd
  ##standardized error
  (leave_one_out_m1$mean-output)/leave_one_out_m1$sd

  
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
