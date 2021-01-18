\name{findInertInputs}
\alias{findInertInputs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
find inert inputs with the posterior mode
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
The function tests for inert inputs (inputs that barely affect the outputs) using the posterior mode. 
}
\usage{
findInertInputs(object,threshold=0.1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object }{
%%     ~~Describe \code{param} here~~
an object of  class \code{rgasp} or the \code{ppgasp}.
}
   \item{threshold}{
   a threshold between 0 to 1. If the normalized inverse parameter of an input is smaller this value, it is classified as inert inputs. 
   
   }
}

\details{
%  ~~ If necessary, more details than the description above ~~expand here the details.

This function utilizes the following quantity 

object@p*object@beta_hat*object@CL/sum(object@beta_hat*object@CL)

for each input to identify the inert outputs. The average estimated normalized inverse range parameters will be 1. If the estimated normalized inverse range parameters of an input is close to 0, it means this input might be an inert input. 

In this method, a prior that has shrinkage effects is suggested to use, .e.g the jointly robust prior (i.e. one should set prior_choice='ref_approx' in rgasp() to obtain the use code{rgasp} object before using this function). Moreover, one may not add a lower bound of the range parameters to perform this method (i.e. one should set lower_bound=F in rgasp()).  For more details see Chapter 4 in 

Mengyang Gu. (2016). Robust Uncertainty Quantification and Scalable Computation for Computer Models with Massive Output. Ph.D. thesis. Duke University.

}


\value{
A vector that has the same dimension of the number of inputs indicating how likely the inputs are inerts. The average value is 1. When a value is very close to zero, it tends to be an inert inputs. 
}
\references{
Mengyang Gu. (2016). Robust Uncertainty Quantification and Scalable Computation for Computer Models with Massive Output. Ph.D. thesis. Duke University.
}
\author{
%%  ~~who you are~~
\packageAuthor{RobustGaSP}

Maintainer: \packageMaintainer{RobustGaSP}
}

\examples{
  #-----------------------------------------------
  # test for inert inputs in the Borehole function
  #-----------------------------------------------
# dimensional of the inputs
dim_inputs <- 8    
# number of the inputs
num_obs <- 40       

# uniform samples of design
set.seed(0)
input <-matrix(runif(num_obs*dim_inputs), num_obs,dim_inputs) 
# Following codes use maximin Latin Hypercube Design, which is typically better than uniform
# library(lhs)
# input <- maximinLHS(n=num_obs, k=dim_inputs)  # maximin lhd sample

# rescale the design to the domain
input[,1]<-0.05+(0.15-0.05)*input[,1];
input[,2]<-100+(50000-100)*input[,2];
input[,3]<-63070+(115600-63070)*input[,3];
input[,4]<-990+(1110-990)*input[,4];
input[,5]<-63.1+(116-63.1)*input[,5];
input[,6]<-700+(820-700)*input[,6];
input[,7]<-1120+(1680-1120)*input[,7];
input[,8]<-9855+(12045-9855)*input[,8];

# outputs from the 8 dim Borehole function

output=matrix(0,num_obs,1)
for(i in 1:num_obs){
  output[i]=borehole(input[i,])
}





# use constant mean basis with trend, with no constraint on optimization
m3<- rgasp(design = input, response = output, lower_bound=FALSE)

P=findInertInputs(m3)


}

%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%\code{\link{neg_log_marginal_post_approx_ref}}.
%}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{internal} % use one of  RShowDoc("KEYWORDS")
