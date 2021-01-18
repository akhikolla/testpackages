\name{as.S4prediction}
\alias{as.S4prediction}
\alias{as.S4prediction.predict}
\alias{as.S4prediction.predrgasp-class}
\alias{as.S4prediction,predrgasp-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Convert a \code{RobustGaSP} \code{S3 object} prediction into a \code{S4 object}
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function converts a \code{S3 object} prediction into the default \code{S4 object} prediction
}
\usage{
as.S4prediction(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{  an object of type list obtained by \code{\link{predict.rgasp}} contains the prediction and it will be converted into a \code{S4 object}
%%     ~~Describe \code{object} here~~
}
  \item{\dots}{ Extra arguments to be passed to the function (not implemented yet).
%%     ~~Describe \code{\dots} here~~
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
The returned value is a \code{S4 object} of class \code{\link[=predrgasp-class]{predrgasp-class}} with 
%\section{Slots}{
	 \describe{
	  \item{\code{call}:}{ \code{call} to the function.}
    \item{\code{mean}:}{ predictive mean for the testing inputs.}
    \item{\code{lower95}:}{lower bound of the 95\% posterior credible interval.}
    \item{\code{upper95}:}{upper bound of the 95\% posterior credible interval.}
    \item{\code{sd}:}{standard deviation of each \code{testing_input}.}
  }
%}
}

\author{
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
  # m1<- rgasp(design = input, response = output, lower_bound=FALSE, multiple_starts=FALSE)
  
  # number of points to be predicted 
  num_testing_input <- 5000    
  # generate points to be predicted
  testing_input <- matrix(runif(num_testing_input*dim_inputs),num_testing_input,dim_inputs)
  # Perform prediction
  m1.predict<-predict(m1, testing_input, outasS3 = FALSE)
  # Notice the call slot of the object
  print(m1.predict@call)
  
  # To convert the prediction to a S3 object 
  m1.predict.aslist <- as.S3prediction(m1.predict)
  # To recover back the prediction as a predrgasp-class object
  m1.predict.aspredgasp <- as.S4prediction.predict(m1.predict.aslist)
  str(m1.predict.aslist)
  # Notice that in this case the @call slot is different than the initial
  print(m1.predict.aspredgasp@call)
  
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

\keyword{internal}