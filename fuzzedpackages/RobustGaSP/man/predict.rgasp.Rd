\name{predict.rgasp}
\alias{predict}
\alias{predict.rgasp}
\alias{predict.rgasp-class}
\alias{predict,rgasp-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Prediction for Robust GaSP model
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Function to make prediction on the robust GaSP model after the robust GaSP model has been constructed.
}
\usage{
\S4method{predict}{rgasp}(object,testing_input,testing_trend= matrix(1,dim(testing_input)[1],1),
r0=NA,interval_data=T,
outasS3 = T,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ an object of  class \code{rgasp}.}
  \item{testing_input}{a matrix containing the inputs where the \code{rgasp} is to perform prediction.}
  \item{testing_trend}{a matrix of mean/trend for prediction.}
  \item{r0}{
  the distance between input and testing input. If the value is \code{NA}, it will be calculated later. It can also be specified by the user. If specified by user, it is either a \code{matrix} or \code{list}. The default value is \code{NA}.
  }
\item{interval_data}{
   a boolean value. If \code{T}, the interval of the data will be calculated. Otherwise, the interval of the mean of the data will be calculted.
}
  \item{outasS3}{a boolean parameter indicating whether the output of the function should be as an \code{S3 object}.}
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
If the parameter \code{outasS3=F}, then the returned value is a \code{S4 object} of class \code{\link[=predrgasp-class]{predrgasp-class}} with 
%\section{Slots}{
	 \describe{
	  \item{\code{call}:}{ \code{call} to \code{predict.rgasp} function where the returned object has been created.}
    \item{\code{mean}:}{ predictive mean for the testing inputs.}
    \item{\code{lower95}:}{lower bound of the 95\% posterior credible interval.}
    \item{\code{upper95}:}{upper bound of the 95\% posterior credible interval.}
    \item{\code{sd}:}{standard deviation of each \code{testing_input}.}
  }
%}
If the parameter \code{outasS3=T}, then the returned value is a \code{list} with 
\item{mean }{ predictive mean for the testing inputs.}
\item{lower95 }{lower bound of the 95\% posterior credible interval.}
\item{upper95 }{upper bound of the 95\% posterior credible interval.}
\item{sd }{standard deviation of each \code{testing_input}.}
}
\references{
M. Gu. (2016). Robust Uncertainty Quantification and Scalable Computation for Computer Models with Massive Output. Ph.D. thesis. Duke University.

M. Gu. and J.O. Berger (2016). Parallel partial Gaussian process emulation for computer models with massive output. \emph{Annals of Applied Statistics}, 10(3), 1317-1347.

M. Gu, X. Wang and J.O. Berger (2018), \emph{Robust Gaussian Stochastic Process Emulation}, \emph{Annals of Statistics}, 46(6A), 3038-3066.

M. Gu (2018), \emph{Jointly Robust Prior for Gaussian Stochastic Process in Emulation, Calibration and Variable Selection}, arXiv:1804.09329.

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
  # m1<- rgasp(design = input, response = output, lower_bound=FALS)
  
  # number of points to be predicted 
  num_testing_input <- 5000    
  # generate points to be predicted
  testing_input <- matrix(runif(num_testing_input*dim_inputs),num_testing_input,dim_inputs)
  # Perform prediction
  m1.predict<-predict(m1, testing_input)
  # Predictive mean
  m1.predict$mean  
  
  # The following tests how good the prediction is 
  testing_output <- matrix(0,num_testing_input,1)
  for(i in 1:num_testing_input){
    testing_output[i]<-dettepepel.3.data(testing_input[i,])
  }
  
  # compute the MSE, average coverage and average length
  # out of sample MSE
  MSE_emulator <- sum((m1.predict$mean-testing_output)^2)/(num_testing_input)  
  
  # proportion covered by 95\% posterior predictive credible interval
  prop_emulator <- length(which((m1.predict$lower95<=testing_output)
                   &(m1.predict$upper95>=testing_output)))/num_testing_input
  
  # average length of  posterior predictive credible interval
  length_emulator <- sum(m1.predict$upper95-m1.predict$lower95)/num_testing_input
  
  # output of prediction
  MSE_emulator
  prop_emulator
  length_emulator  
  # normalized RMSE
  sqrt(MSE_emulator/mean((testing_output-mean(output))^2 ))


  #-----------------------------------
  # a 2 dimensional example with trend
  #-----------------------------------
  # dimensional of the inputs
  dim_inputs <- 2    
  # number of the inputs
  num_obs <- 20       
  
  # uniform samples of design
  input <-matrix(runif(num_obs*dim_inputs), num_obs,dim_inputs) 
  # Following codes use maximin Latin Hypercube Design, which is typically better than uniform
  # library(lhs)
  # input <- maximinLHS(n=num_obs, k=dim_inputs)  ##maximin lhd sample
  
  # outputs from the 2 dim Brainin function
  
  output <- matrix(0,num_obs,1)
  for(i in 1:num_obs){
    output[i]<-limetal.2.data (input[i,])
  }
  
  #mean basis (trend)
  X<-cbind(rep(1,num_obs), input )
  
  
  # use constant mean basis with trend, with no constraint on optimization
  m2<- rgasp(design = input, response = output,trend =X,  lower_bound=FALSE)
  
  
  # number of points to be predicted 
  num_testing_input <- 5000    
  # generate points to be predicted
  testing_input <- matrix(runif(num_testing_input*dim_inputs),num_testing_input,dim_inputs)
  
  # trend of testing
  testing_X<-cbind(rep(1,num_testing_input), testing_input )
  
  
  # Perform prediction
  m2.predict<-predict(m2, testing_input,testing_trend=testing_X)
  # Predictive mean
  #m2.predict$mean  
  
  # The following tests how good the prediction is 
  testing_output <- matrix(0,num_testing_input,1)
  for(i in 1:num_testing_input){
    testing_output[i]<-limetal.2.data(testing_input[i,])
  }
  
  # compute the MSE, average coverage and average length
  # out of sample MSE
  MSE_emulator <- sum((m2.predict$mean-testing_output)^2)/(num_testing_input)  
  
  # proportion covered by 95\% posterior predictive credible interval
  prop_emulator <- length(which((m2.predict$lower95<=testing_output)
                   &(m2.predict$upper95>=testing_output)))/num_testing_input
  
  # average length of  posterior predictive credible interval
  length_emulator <- sum(m2.predict$upper95-m2.predict$lower95)/num_testing_input
  
  # output of prediction
  MSE_emulator
  prop_emulator
  length_emulator  
  # normalized RMSE
  sqrt(MSE_emulator/mean((testing_output-mean(output))^2 ))


    ###here try the isotropic kernel (a function of Euclidean distance)
  m2_isotropic<- rgasp(design = input, response = output,trend =X,  
             lower_bound=FALSE,isotropic=TRUE)
  
  m2_isotropic.predict<-predict(m2_isotropic, testing_input,testing_trend=testing_X)
  
  # compute the MSE, average coverage and average length
  # out of sample MSE
  MSE_emulator_isotropic <- sum((m2_isotropic.predict$mean-testing_output)^2)/(num_testing_input)
  
  # proportion covered by 95\% posterior predictive credible interval
  prop_emulator_isotropic <- length(which((m2_isotropic.predict$lower95<=testing_output)
                                &(m2_isotropic.predict$upper95>=testing_output)))/num_testing_input
  
  # average length of  posterior predictive credible interval
  length_emulator_isotropic <- sum(m2_isotropic.predict$upper95-
  m2_isotropic.predict$lower95)/num_testing_input
  
  MSE_emulator_isotropic
  prop_emulator_isotropic
  length_emulator_isotropic
  ##the result of isotropic kernel is not as good as the product kernel for this example


  #--------------------------------------------------------------------------------------
  # an 8 dimensional example using only a subset inputs and a noise with unknown variance
  #--------------------------------------------------------------------------------------
  set.seed(1)
  # dimensional of the inputs
  dim_inputs <- 8    
  # number of the inputs
  num_obs <- 50       
  
  # uniform samples of design
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
  m3<- rgasp(design = input[,c(1,4,6,7,8)], response = output,
             nugget.est=TRUE, lower_bound=FALSE)
  
  
  # number of points to be predicted 
  num_testing_input <- 5000    
  # generate points to be predicted
  testing_input <- matrix(runif(num_testing_input*dim_inputs),num_testing_input,dim_inputs)
  
  # resale the points to the region to be predict
  testing_input[,1]<-0.05+(0.15-0.05)*testing_input[,1];
  testing_input[,2]<-100+(50000-100)*testing_input[,2];
  testing_input[,3]<-63070+(115600-63070)*testing_input[,3];
  testing_input[,4]<-990+(1110-990)*testing_input[,4];
  testing_input[,5]<-63.1+(116-63.1)*testing_input[,5];
  testing_input[,6]<-700+(820-700)*testing_input[,6];
  testing_input[,7]<-1120+(1680-1120)*testing_input[,7];
  testing_input[,8]<-9855+(12045-9855)*testing_input[,8];
  
  
  # Perform prediction
  m3.predict<-predict(m3, testing_input[,c(1,4,6,7,8)])
  # Predictive mean
  #m3.predict$mean  
  
  # The following tests how good the prediction is 
  testing_output <- matrix(0,num_testing_input,1)
  for(i in 1:num_testing_input){
    testing_output[i]<-borehole(testing_input[i,])
  }
  
  # compute the MSE, average coverage and average length
  # out of sample MSE
  MSE_emulator <- sum((m3.predict$mean-testing_output)^2)/(num_testing_input)  
  
  # proportion covered by 95\% posterior predictive credible interval
  prop_emulator <- length(which((m3.predict$lower95<=testing_output)
                   &(m3.predict$upper95>=testing_output)))/num_testing_input
  
  # average length of  posterior predictive credible interval
  length_emulator <- sum(m3.predict$upper95-m3.predict$lower95)/num_testing_input
  
  # output of sample prediction
  MSE_emulator
  prop_emulator
  length_emulator  
  # normalized RMSE
  sqrt(MSE_emulator/mean((testing_output-mean(output))^2 ))

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
