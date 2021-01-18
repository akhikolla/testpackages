\name{predict.ppgasp}
\alias{predict.ppgasp}
\alias{predict.ppgasp-class}
\alias{predict,ppgasp-method}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Prediction for PP GaSP model
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Function to make prediction on the PP GaSP model after the PP GaSP model has been constructed.
}
\usage{
\S4method{predict}{ppgasp}(object, testing_input, 
testing_trend= matrix(1,dim(testing_input)[1],1),r0=NA, 
interval_data=T,
outasS3 = T, ...)
%- maybe also 'usage' for other objects documented here.
}
\arguments{
  \item{object}{ an object of  class \code{ppgasp}.}
  \item{testing_input}{a matrix containing the inputs where the \code{rgasp} is to perform prediction.}
  \item{testing_trend}{a matrix of mean/trend for prediction.}
  \item{r0}{
      the distance between input and testing input. If the value 
      is \code{NA}, it will be calculated later. It can also be 
      specified by the user. If specified by user, it is either a 
      \code{matrix} or \code{list}. The default value is \code{NA}.
    }
\item{interval_data}{
   a boolean value. If \code{T}, the interval of the data will be   
   calculated. Otherwise, the interval of the mean of the data will 
   be calculted.
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
If the parameter \code{outasS3=F}, then the returned value is a \code{S4 object} of class \code{\link[=predppgasp-class]{predppgasp-class}} with 

%	 \describe{
	  \item{\code{call}:}{ \code{call} to \code{predict.ppgasp} function where the returned object has been created.}
    \item{\code{mean}:}{ predictive mean for the testing inputs.}
    \item{\code{lower95}:}{lower bound of the 95\% posterior credible interval.}
    \item{\code{upper95}:}{upper bound of the 95\% posterior credible interval.}
    \item{\code{sd}:}{standard deviation of each \code{testing_input}.}

If the parameter \code{outasS3=T}, then the returned value is a \code{list} with 
\item{mean }{ predictive mean for the testing inputs.}
\item{lower95 }{lower bound of the 95\% posterior credible interval.}
\item{upper95 }{upper bound of the 95\% posterior credible interval.}
\item{sd }{standard deviation of each \code{testing_input}.}
}

\references{

M. Gu. and J.O. Berger (2016). Parallel partial Gaussian process emulation for computer models with massive output. \emph{Annals of Applied Statistics}, 10(3), 1317-1347.

M. Gu. (2016). Robust Uncertainty Quantification and Scalable Computation for Computer Models with Massive Output. Ph.D. thesis. Duke University.

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
  #----------------------------------
  # an example of environmental model
  #----------------------------------
  
  set.seed(1)
  #Here the sample size is very small. Consider to use more observations 
  n=80
  p=4
  ##using the latin hypercube will be better
  #library(lhs)
  #input_samples=maximinLHS(n,p)
  input_samples=matrix(runif(n*p),n,p)
  input=matrix(0,n,p)
  input[,1]=7+input_samples[,1]*6
  input[,2]=0.02+input_samples[,2]*1
  input[,3]=0.01+input_samples[,3]*2.99
  input[,4]=30.01+input_samples[,4]*0.285
  
  k=300
  output=matrix(0,n,k)
  ##environ.4.data is an environmental model on a spatial-time vector
  ##? environ.4.data
  for(i in 1:n){
    output[i,]=environ.4.data(input[i,],s=seq(0.15,3,0.15),t=seq(4,60,4)  )
  }
  
  ##samples some test inputs
  n_star=1000
  sample_unif=matrix(runif(n_star*p),n_star,p)
  
  testing_input=matrix(0,n_star,p)
  testing_input[,1]=7+sample_unif[,1]*6
  testing_input[,2]=0.02+sample_unif[,2]*1
  testing_input[,3]=0.01+sample_unif[,3]*2.99
  testing_input[,4]=30.01+sample_unif[,4]*0.285
  
  
  testing_output=matrix(0,n_star,k)
  
  for(i in 1:n_star){
    testing_output[i,]=environ.4.data(testing_input[i,],s=seq(0.15,3,0.15
    ),t=seq(4,60,4) )
  }
  
  ##we do a transformation of the output 
  ##one can change the number of initial values to test
  log_output_1=log(output+1)
  #since we have lots of output, we use 'nelder-mead' for optimization
  m.ppgasp=ppgasp(design=input,response=log_output_1,kernel_type
                  ='pow_exp',num_initial_values=2,optimization='nelder-mead')
  
  m_pred.ppgasp=predict(m.ppgasp,testing_input)
  ##we transform back for the prediction
  m_pred_ppgasp_median=exp(m_pred.ppgasp$mean)-1
  ##mean squared error
  mean( (m_pred_ppgasp_median-testing_output)^2)
  ##variance of the testing outputs
  var(as.numeric(testing_output))
  
  ##makes plots for the testing 
  par(mfrow=c(1,2))
  t=seq(3,60,3)
  s=seq(0.15,3,0.15)
  testing_plot_1=matrix(testing_output[1,],  length(t), length(s) )
  
  max_testing_plot_1=max(testing_plot_1)
  min_testing_plot_1=min(testing_plot_1)
  
  image(x=t,y=s,testing_plot_1,  col = hcl.colors(100, "terrain"),main='test outputs')
  contour(x=t,y=s,testing_plot_1, levels = seq(min_testing_plot_1, max_testing_plot_1,
                                               by = (max_testing_plot_1-min_testing_plot_1)/5),
          add = TRUE, col = "brown")
  
  ppgasp_plot_1=matrix(m_pred_ppgasp_median[1,],  length(t), length(s) )
  max_ppgasp_plot_1=max(ppgasp_plot_1)
  min_ppgasp_plot_1=min(ppgasp_plot_1)
  
  image(x=t,y=s,ppgasp_plot_1,  col = hcl.colors(100, "terrain"),main='prediction')
  contour(x=t,y=s,ppgasp_plot_1, levels = seq(min_testing_plot_1, max_ppgasp_plot_1,
                                              by = (max_ppgasp_plot_1-min_ppgasp_plot_1)/5),
          add = TRUE, col = "brown")
  dev.off()

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
