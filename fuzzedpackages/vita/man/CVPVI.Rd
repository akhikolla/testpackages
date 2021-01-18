\name{CVPVI}
\alias{CVPVI}
\alias{CVPVI.default}
\alias{print.CVPVI}
\title{
Cross-validated permutation variable importance measure
}
\description{
Compute cross-validated permutation variable importance measure from a random forest for classification and regression.
}
\usage{

\method{CVPVI}{default}(X, y, k = 2, mtry= if (!is.null(y) && !is.factor(y))
                        max(floor(ncol(X)/3), 1) else floor(sqrt(ncol(X))),
    ntree = 500, nPerm = 1, parallel = FALSE, ncores = 0, seed = 123, ...)
\method{print}{CVPVI}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{a data frame or a matrix of predictors. }
  \item{y}{a response vector.}
  \item{k}{an integer for the number of folds. Default is \code{k = 2}}
  \item{mtry}{Number of variables randomly sampled as candidates at each split for the l-th forest. Note that
              the default values are different for classification (\code{mtry=sqrt(p)} where \code{p} is number
              of variables in \code{x}) and regression (\code{mtry=p/3}).}
  \item{ntree}{Number of trees to grow for the l-th forest. Default is \code{ntree=500}.}
  \item{nPerm}{Number of times the l-th data set are permuted per tree for assessing variable fold-specific
                permutation variable importance. Default is \code{nPerm=1}.}
  \item{parallel}{Should the CVPVI implementation run parallel? Default is \code{parallel=FALSE} and the number of cores is
                 set to one. The parallelized version of the CVPVI implementation are based on
                 \code{\link[parallel]{mclapply}} and so are not available on Windows.}
  \item{ncores}{The number of cores to use, i.e. at most how many child processes will be run
                simultaneously. Must be at least one, and parallelization requires at least two cores.
                If \code{ncores=0}, then the half of CPU cores on the current host are used.}
  \item{seed}{a single integer value to specify seeds. The "combined multiple-recursive generator"
              from L'Ecuyer (1999) is set as random number generator for the parallelized version of
              the CVPVI implementation. Default is \code{ seed = 123}.}
  \item{\dots}{optional parameters for \code{\link[randomForest]{randomForest} } }
  \item{x}{for the print method, an \code{CVPVI} object }
}
\details{
This method randomly splits the dataset into k sets of equal size. The method constructs k random forests, where the l-th forest is constructed based on observations that are not part of the l-th set. For each forest the fold-specific permutation variable importance measure is computed using all observations in the l-th data set: For each tree, the prediction error on the l-th data set is recorded. Then the same is done after permuting the values of each predictor variable.
The differences between the two prediction errors are then averaged over all trees. The cross-validated permutation variable importance is the average of all k-fold-specific permutation variable importances. For classification the mean decrease in accuracy over all classes is used and for regression the mean decrease in MSE.
}
\value{
\item{fold_varim }{a p by k matrix of fold-specific permutation variable importances. For classification the mean decrease in accuracy over all classes. For regression mean decrease in MSE. }
\item{cv_varim }{cross-validated permutation variable importances. For classification the mean decrease in accuracy over all classes. For regression mean decrease in MSE. }
\item{type}{one of regression, classification}
}
\references{
Janitza S, Celik E, Boulesteix A-L, (2015),\emph{ A computationally fast variable importance test for random forest for high dimensional data},Technical Report 185, University of Munich, <\url{http://nbn-resolving.de/urn/resolver.pl?urn=nbn:de:bvb:19-epub-25587-4}>
}


\seealso{
\code{\link[vita]{VarImpCVl}}, \code{\link[randomForest]{importance}}, \code{\link[randomForest]{randomForest}}, \code{\link[parallel]{mclapply}}
}
\examples{
##############################
#      Classification        #
##############################
## Simulating data
X = replicate(10,rnorm(100))
X= data.frame( X) #"X" can also be a matrix
z  = with(X,5*X1 + 3*X2 + 2*X3 + 1*X4 -
            5*X5 - 9*X6 - 2*X7 + 1*X8 )
pr = 1/(1+exp(-z))         # pass through an inv-logit function
y = as.factor(rbinom(100,1,pr))
##################################################################
# cross-validated permutation variable importance
cv_vi = CVPVI(X,y,k = 2,mtry = 3,ntree = 1000,ncores = 4)
print(cv_vi)

##################################################################
#compare them with the original permutation variable importance
library("randomForest")
cl.rf = randomForest(X,y,mtry = 3,ntree = 1000, importance = TRUE)

round(cbind(importance(cl.rf, type=1, scale=FALSE),cv_vi$cv_varim),digits=5)


###############################
#      Regression            #
##############################

##################################################################
## Simulating data:
X = replicate(10,rnorm(100))
X = data.frame( X) #"X" can also be a matrix
y = with(X,2*X1 + 2*X2 + 2*X3 + 1*X4 - 2*X5 - 2*X6 - 1*X7 + 2*X8 )

##################################################################
# cross-validated permutation variable importance
cv_vi = CVPVI(X,y,k = 3,mtry = 3,ntree = 1000,ncores = 2)
print(cv_vi)
##################################################################
#compare them with the original permutation variable importance
library("randomForest")
reg.rf = randomForest(X,y,mtry = 3,ntree = 1000, importance = TRUE)

round(cbind(importance(reg.rf, type=1, scale=FALSE),cv_vi$cv_varim),digits=5)

}

