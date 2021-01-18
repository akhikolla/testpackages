\name{RobustGaSP-package}
\alias{RobustGaSP-package}
\alias{RobustGaSP}
\docType{package}
\title{
\packageTitle{RobustGaSP}
}
\description{
\packageDescription{RobustGaSP}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{RobustGaSP}
\packageIndices{RobustGaSP}
%%~~ An overview of how to use the package, including the most important functions ~~
}
\author{
\packageAuthor{RobustGaSP}

Maintainer: \packageMaintainer{RobustGaSP}
}

\references{

J.O. Berger, V. De Oliveira and B. Sanso (2001), \emph{Objective Bayesian analysis of spatially correlated data}, \emph{Journal of the American Statistical Association}, {96}, 1361-1374.

M. Gu. and J.O. Berger (2016). Parallel partial Gaussian process emulation for computer models with massive output. \emph{Annals of Applied Statistics}, 10(3), 1317-1347.

M. Gu. (2016). Robust uncertainty quantification and scalable computation for computer models with massive output. Ph.D. thesis. Duke University.

M. Gu, X. Wang and J.O. Berger (2018), \emph{Robust Gaussian stochastic process emulation}, \emph{Annals of Statistics}, 46(6A), 3038-3066.

M. Gu (2018), \emph{Jointly robust prior for Gaussian stochastic process in emulation, calibration and variable selection}, arXiv:1804.09329.

R. Paulo (2005), \emph{Default priors for Gaussian processes}, \emph{Annals of statistics}, 33(2), 556-582.

J. Sacks, W.J. Welch, T.J. Mitchell, and H.P. Wynn (1989), \emph{Design and analysis of computer experiments}, \emph{Statistical Science}, \bold{4}, 409-435.



}

%%~~ Optionally other standard keywords, one per line, from file KEYWORDS in the R documentation directory ~~
\keyword{package}
\keyword{computer model}
\keyword{emulation}
\keyword{simulation}
%\seealso{
%%~~ Optional links to other man pages, e.g. ~~
%%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
%\code{\link[RobustGaSP]{RobustGaSP}}
%%\code{\link[RobustGaSP]{RobustGaSP-class}}
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
  
  ####
  # outputs from the 3 dim dettepepel.3.data function
  
  output = matrix(0,num_obs,1)
  for(i in 1:num_obs){
    output[i]<-dettepepel.3.data(input[i,])
  }
  
  # use constant mean basis, with no constraint on optimization
  m1<- rgasp(design = input, response = output, lower_bound=FALSE)
  
  # the following use constraints on optimization
  # m1<- rgasp(design = input, response = output, lower_bound=TRUE)
  
  # the following use a single start on optimization
  # m1<- rgasp(design = input, response = output, lower_bound=FALSE)
  
  # number of points to be predicted 
  num_testing_input <- 5000    
  # generate points to be predicted
  testing_input <- matrix(runif(num_testing_input*dim_inputs),num_testing_input,dim_inputs)
  # Perform prediction
  m1.predict<-predict(m1, testing_input, outasS3 = FALSE)
  # Predictive mean
  m1.predict@mean  
  
  # The following tests how good the prediction is 
  testing_output <- matrix(0,num_testing_input,1)
  for(i in 1:num_testing_input){
    testing_output[i]<-dettepepel.3.data(testing_input[i,])
  }
  
  # compute the MSE, average coverage and average length
  # out of sample MSE
  MSE_emulator <- sum((m1.predict@mean-testing_output)^2)/(num_testing_input)  
  
  # proportion covered by 95\% posterior predictive credible interval
  prop_emulator <- length(which((m1.predict@lower95<=testing_output)
                   &(m1.predict@upper95>=testing_output)))/num_testing_input
  
  # average length of  posterior predictive credible interval
  length_emulator <- sum(m1.predict@upper95-m1.predict@lower95)/num_testing_input
  
  # output of prediction
  MSE_emulator
  prop_emulator
  length_emulator  
  # normalized RMSE
  sqrt(MSE_emulator/mean((testing_output-mean(output))^2 ))
}
