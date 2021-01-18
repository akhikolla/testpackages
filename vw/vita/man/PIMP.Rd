\name{PIMP}
\alias{PIMP}
\alias{PIMP.default}
\alias{print.PIMP}
\title{
PIMP-algorithm for the permutation variable importance measure
}
\description{
\code{PIMP} implements the test approach of Altmann et al. (2010) for the permutation variable importance measure \code{VarImp}
in a random forest for classification and regression.
}
\usage{
\method{PIMP}{default}(X, y, rForest, S = 100, parallel = FALSE, ncores=0, seed = 123, ...)
\method{print}{PIMP}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{a data frame or a matrix of predictors}
  \item{y}{a response vector. If a factor, classification is assumed,
              otherwise regression is assumed.}
  \item{rForest}{an object of class \code{\link[randomForest]{randomForest}}, \code{importance} must
                be set to True.}
  \item{S}{The number of permutations for the response vector \code{y}. Default is \code{S=100}. }
  \item{parallel}{Should the PIMP-algorithm run parallel?  Default is \code{parallel=FALSE} and the number of cores is
                 set to one. The parallelized version of the PIMP-algorithm are based on
                 \code{\link[parallel]{mclapply}} and so is not available on Windows.}
  \item{ncores}{The number of cores to use, i.e. at most how many child processes will be run
                simultaneously. Must be at least one, and parallelization requires at least two cores.
                If \code{ncores=0}, then the half of CPU cores on the current host are used.}
   \item{seed}{a single integer value to specify seeds. The "combined multiple-recursive generator"
              from L'Ecuyer (1999) is set as random number generator for the parallelized version of
              the PIMP-algorithm.  Default is \code{ seed = 123}.}
  \item{\dots}{optional parameters for \code{\link[randomForest]{randomForest} } }
  \item{x}{for the print method, an \code{PIMP} object }
}
\details{
The PIMP-algorithm by Altmann et al. (2010) permutes \eqn{S} times the response variable \eqn{y}.
For each permutation of the response vector \eqn{y^{*s}}, a new forest is grown and the permutation
variable importance measure (\eqn{VarImp^{*s}}) for all predictor variables \eqn{X} is computed.
The vector \code{perVarImp} of \eqn{S} \eqn{VarImp} measures for every predictor variables are used
to approximate the null importance distributions (\code{\link[vita]{PimpTest}}).
}
\value{
\item{VarImp }{ the \emph{original permutation variable importance} measures of the random forest. }
\item{PerVarImp }{a matrix, where each row is a vector containing the \code{S} permuted VarImp
                  measures for each predictor variables. }
\item{type}{one of regression, classification}
}
\references{
 Breiman L. (2001), \emph{Random Forests}, Machine Learning 45(1),5-32, <doi:10.1023/A:1010933404324>

Altmann A.,Tolosi L., Sander O. and  Lengauer T. (2010),\emph{Permutation importance: a corrected feature importance measure}, Bioinformatics Volume 26 (10), 1340-1347, <doi:10.1093/bioinformatics/btq134>
}

\seealso{
 \code{\link[vita]{PimpTest}}, \code{\link[randomForest]{importance}}, \code{\link[randomForest]{randomForest}}, \code{\link[parallel]{mclapply}}
}
\examples{

###############################
#      Regression            #
##############################
##############################
## Simulating data
X = replicate(12,rnorm(100))
X = data.frame(X) #"X" can also be a matrix
y = with(X,2*X1 + 1*X2 + 2*X3 + 1*X4 - 2*X5 - 1*X6 - 1*X7 + 2*X8 )

##############################
## Regression with Random Forest:
library("randomForest")
reg.rf = randomForest(X,y,mtry = 3,ntree=500,importance=TRUE)
##############################
## PIMP-Permutation variable importance measure
# the parallelized version of the PIMP-algorithm
system.time(pimp.varImp.reg<-PIMP(X,y,reg.rf,S=10, parallel=TRUE, ncores=2))
# the non parallelized version of the PIMP-algorithm
system.time(pimp.varImp.reg<-PIMP(X,y,reg.rf,S=10, parallel=FALSE))

##############################
#      Classification        #
##############################
## Simulating data
X = replicate(12,rnorm(100))
X= data.frame( X) #"X" can also be a matrix
z  = with(X,2*X1 + 3*X2 + 2*X3 + 1*X4 -
            2*X5 - 2*X6 - 2*X7 + 1*X8 )
pr = 1/(1+exp(-z))         # pass through an inv-logit function
y = as.factor(rbinom(100,1,pr))

##############################
## Classification with Random Forest:
cl.rf = randomForest(X,y,mtry = 3,ntree = 500, importance = TRUE)
##############################
## PIMP-Permutation variable importance measure
# the parallelized version of the PIMP-algorithm
system.time(pimp.varImp.cl<-PIMP(X,y,cl.rf,S=10, parallel=TRUE, ncores=2))
# the non parallelized version of the PIMP-algorithm
system.time(pimp.varImp.cl<-PIMP(X,y,cl.rf,S=10, parallel=FALSE))

}


