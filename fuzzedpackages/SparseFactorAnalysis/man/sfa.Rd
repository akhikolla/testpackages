\name{sfa}
\alias{sfa}
\title{Sparse factor analysis for mixed binary and count data.}
\description{Scaling mixed binary and count data while estimating the underlying latent dimensionality.}
\usage{
sfa(M, missing.mat=NULL, gibbs=100, burnin=100, max.optim=50, 
  thin=1, save.curr="UDV_curr", save.each=FALSE, thin.save=25, 
  maxdim=NULL)
}
\arguments{
  \item{M}{Matrix to be scaled.}
\item{missing.mat}{Matrix indicating missing data.  Should be the same size as \env{M}, with a 1 denoting a missing observation and a 0 otherwise.  Defaults to all zeroes.}
\item{gibbs}{Number of posterior samples to draw}
\item{burnin}{
Number of burnin samples.}
\item{max.optim}{Number of iterations to fit the cutpoints using \env{optim}.  This is generally faster than the Hamiltonian Monte Carlo estimates, and is useful for the first part of the burnin phase.}
\item{thin}{
Extent of thinning of the MCMC chain.  Only every \env{thin} draw is saved to the output.
}
\item{save.curr}{Name of file in which to save object.
}
\item{save.each}{Whether to save with a new name at each thinned draw.
}
\item{thin.save}{How many thinned draws to wait between saving output.}
\item{maxdim}{Number of latent dimensions to fit.  Should be greater than the number of estimated dimensions.}
}
\details{ The function \env{sfa} is the main function in the package, \env{SparseFactorAnalysis}.  It takes in a matrix which in \emph{rows} has the same data type--either binary or count.  For example, every row may consist of roll call votes or word counts, and the columns may correspond with legislators.  The method combines the two data types, scales both, and selects the underlying latent dimensionality.
	}

\value{
  \item{dim.sparse}{Output for sparse estimates of dimensionality.}
  \item{dim.mean}{Non-sparse estimates of posterior mean of dimensionality.}
  \item{rowdim1}{Posterior samples of first dimension of spatial locations for each observation i.}
  \item{rowdim2}{Posterior samples of second dimension of spatial locations for row unit of observation.}
  \item{coldim1}{Posterior samples of first dimension of spatial locations for column unit of observation.}
  \item{coldim2}{Posterior samples of second dimension of spatial locations for column unit of observation.}
  \item{lambda.lasso}{Posterior samples for tuning parameter used for dimension selection.}
  \item{Z}{Posterior mean of fitted values, on a z-scale.}
  \item{rowdims.all}{Posterior mean of all row spatial locations.}
  \item{coldims.all}{Posterior mean of all column spatial locations.}
 }
\author{
Marc Ratkovic and Yuki Shiraito
}
\references{
In Song Kim, John Londregan, and Marc Ratkovic.  2015.  "Voting, Speechmaking, and the Dimensions of Conflict in the US Senate."  Working paper.
}
\seealso{\link{plot.sfa}, \link{summary.sfa}}
\examples{

\dontrun{
##Sample size and dimensions.
 set.seed(1)
 n.sim<-50
 k.sim<-500
 
##True vector of dimension weights.
 d.sim<-rep(0,n.sim)
 d.sim[1:3]<-c(2, 1.5, 1)*3

##Formulate true latent dimensions.
 U.sim<-matrix(rnorm(n.sim^2,sd=.5), nr=n.sim, nc=n.sim)
 V.sim<-matrix(rnorm(n.sim*k.sim,sd=.5), nr=k.sim, nc=n.sim)
 Theta.sim<-U.sim\%*\%diag(d.sim)\%*\%t(V.sim)

##Generate binary outcome and count data.
 probs.sim<-pnorm((-1+Theta.sim+rep(1,n.sim)\%*\%t(rnorm(k.sim,sd=.5)) + 
   rnorm(n.sim,sd=.5)\%*\%t(rep(1,k.sim))   ))
 votes.mat<- 
    apply(probs.sim[1:25,],c(1,2),FUN=function(x) rbinom(1,1,x))
 count.mat<- 
    apply(probs.sim[26:50, ],c(1,2),FUN=function(x) rpois(1,20*x))
 M<-rbind(votes.mat,count.mat)
 
## Run sfa
 sparse1<-sfa(M, maxdim=10)
 
##Analyze results.
 summary(sparse1)
 plot(sparse1,type="dim")
 plot(sparse1,type="scatter")

##Compare to true data generating process

plot(sparse1$Z,Theta.sim)
abline(c(0,1))

}

}
