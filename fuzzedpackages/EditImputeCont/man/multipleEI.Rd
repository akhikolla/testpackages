\name{multipleEI}
\alias{multipleEI}
\title{Generate Multiple Edit-imputed Datasets}
\description{
multipleEI returns \eqn{m} multiple edit-imputed datasets.
}
\usage{
multipleEI(model.obj, n.burnin, m, int.btw.EI, show.iter=TRUE)
}
\arguments{
	\item{model.obj}{\code{\link{Rcpp_bei}} model-object generated from \code{\link{createModel}}.}
	
	\item{n.burnin}{number of burn in iterations.}
	
	\item{m}{number of multiple edit-imputed datasets.}
		
	\item{int.btw.EI}{interval (number of iterations) between EI datasets.}
			
	\item{show.iter}{\code{logical} specifying if the iteration number of burning-in is displayed.}
}
\details{
The total number of MCMC iterations is (\code{n.burnin} + \code{m} * \code{int.btw.EI}). Please see the example in the demo folder for more detailed explanation. 
}
\value{
\code{array} of (\eqn{m}, \eqn{n}, \eqn{p}) dimension where \eqn{m} is the number of edit-imputed data sets, \eqn{n} is the number of records and \eqn{p} is the number of variables.
}
\seealso{\code{\link{createModel}}} 
\examples{
data(SimpleEx)
data1 = readData(SimpleEx$D.obs, SimpleEx$Ratio.edit, NULL, SimpleEx$Balance.edit)
# model1 = createModel(data1, 15)

## get 3 edit-imputed data from MCMC by storing every 100 iterations after 50 burn-in

# result1 = multipleEI(model1, n.burnin=50, m=3, int.btw.EI=100)

# dim(result1) 
## [1]   3 1000   4  
## m=3 Edit-imputed datasets of n=1000 records with p=4 variables
}
% \keyword{internal}