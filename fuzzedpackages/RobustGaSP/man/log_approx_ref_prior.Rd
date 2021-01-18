\name{log_approx_ref_prior}
\alias{log_approx_ref_prior}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
The natural logarithm of the jointly robust prior (up to a normalizing constant)
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
A function to the natural logarithm of the jointly robust prior (up to a normalizing constant).}
\usage{
log_approx_ref_prior(param, nugget, nugget_est, CL, a, b)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{param}{
%%     ~~Describe \code{param} here~~
A vector of natural logarithm of inverse-range parameters and natural logarithm of the nugget-variance ratio parameter.
}
  \item{nugget}{
%%     ~~Describe \code{nugget} here~~
The nugget-variance ratio parameter if this parameter is fixed.
}
  \item{nugget_est}{
%%     ~~Describe \code{nugget.est} here~~
Boolean value of whether the nugget is estimated or fixed.
}
  \item{CL}{
%%     ~~Describe \code{CL} here~~
Prior parameter in the jointly robust prior.
}
  \item{a}{
%%     ~~Describe \code{a} here~~
Prior parameter in the jointly robust prior.
}
  \item{b}{
%%     ~~Describe \code{b} here~~
Prior parameter in the jointly robust prior.
}
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
The numerical value of the derivative of the approximate reference prior with regard to inverse-range parameters and the nugget-variance ratio parameter. When the nugget is fixed, the derivative is on inverse-range parameters.
}
\references{
M. Gu. (2016). Robust uncertainty quantification and scalable computation for computer models with massive output. Ph.D. thesis. Duke University.

M. Gu (2018), \emph{Jointly robust prior for Gaussian stochastic process in emulation, calibration and variable selection}, arXiv:1804.09329.

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

\seealso{\code{\link{rgasp}} }

\examples{
# inputs
x<-runif(10);
n<-length(x);

# default prior parameters
a<-0.2
b<-n^{-1}*(a+1)
R0<-as.matrix(abs(outer(x,x, "-")))
CL<- mean(R0[which(R0>0)])

# compute the density of log reference prior up to a normalizing constant
param <- seq(-10,10,0.01)
prior <- rep(0,length(param))
for(i in 1:length(param)){
  prior[i] <- exp(log_approx_ref_prior(param[i],nugget=0,nugget_est=FALSE,CL,a,b) )
}
# plot
plot(param,prior,type='l',
                xlab='Logarithm of inverse range parameters',
                ylab='Prior density up to a normalizing constant')
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")

\keyword{internal}