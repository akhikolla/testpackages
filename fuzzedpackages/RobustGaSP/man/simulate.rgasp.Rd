\name{simulate}
\alias{simulate}
\alias{simulate.rgasp}
\alias{simulate.rgasp-class}
\alias{simulate,rgasp-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Sample for Robust GaSP model
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Function to sample Robust GaSP after the Robust GaSP model has been constructed.
}
\usage{
\S4method{simulate}{rgasp}(object, testing_input, num_sample=1,
testing_trend= matrix(1,dim(testing_input)[1],1),
r0=NA,rr0=NA,sample_data=T,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ an object of  class \code{rgasp}.}
  \item{testing_input}{a matrix containing the inputs where the \code{rgasp} is to sample.}
    \item{num_sample}{number of samples one wants. }
  \item{testing_trend}{a matrix of mean/trend for prediction.}
    \item{r0}{
  the distance between input and testing input. If the value is \code{NA}, it will be calculated later. It can also be specified by the user. If specified by user, it is either a \code{matrix} or \code{list}. The default value is \code{NA}.
  }
      \item{rr0}{
  the distance between testing input and testing input. If the value is \code{NA}, it will be calculated later. It can also be specified by the user. If specified by user, it is either a \code{matrix} or \code{list}. The default value is \code{NA}.
  }

\item{sample_data}{
   a boolean value. If \code{T}, the interval of the data will be calculated. Otherwise, the interval of the mean of the data will be calculted.
}
  \item{...}{Extra arguments to be passed to the function (not implemented yet).}
}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%Provide here some details.
%}


\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
The returned value is a \code{matrix} where each column is a sample on the prespecified inputs. 
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
  #------------------------
  # a 1 dimensional example
  #------------------------
  
###########1dim higdon.1.data 
p1 = 1     ###dimensional of the inputs
dim_inputs1 <- p1
n1 = 15   ###sample size or number of training computer runs you have 
num_obs1 <- n1
input1 = 10*matrix(runif(num_obs1*dim_inputs1), num_obs1,dim_inputs1) ##uniform
#####lhs is better
#library(lhs)
#input1 = 10*maximinLHS(n=num_obs1, k=dim_inputs1)  ##maximin lhd sample
output1 = matrix(0,num_obs1,1)
for(i in 1:num_obs1){
  output1[i]=higdon.1.data (input1[i])
}




m1<- rgasp(design = input1, response = output1, lower_bound=FALSE)

#####locations to samples
testing_input1 = seq(0,10,1/50) 
testing_input1=as.matrix(testing_input1)
#####draw 10 samples
m1_sample=simulate(m1,testing_input1,num_sample=10)

#####plot these samples
matplot(testing_input1,m1_sample, type='l',xlab='input',ylab='output')
lines(input1,output1,type='p')


}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
