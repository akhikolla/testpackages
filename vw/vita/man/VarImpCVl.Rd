\name{VarImpCVl}
\alias{VarImpCVl}

\title{
Fold-specific permutation variable importance measure
}
\description{
Compute fold-specific permutation variable importance measure from a random forest for classification and regression.
}
\usage{
VarImpCVl(X_l, y_l, rForest, nPerm = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X_l}{ a data frame or a matrix of predictors from the l-th data set}
  \item{y_l}{ a response vector from the l-th data set. If a factor, classification is assumed,
              otherwise regression is assumed.}
  \item{rForest}{an object of class \code{\link[randomForest]{randomForest}}, keep.forest must be set to True.
                The l-th Forest based on observations that are not part of the l-th data set.}
  \item{nPerm}{ Number of permutations performed per tree for computing fold-specific
                permutation variable importance. Currently only implemented for regression. }
}
\details{
The fold-specific permutation variable importance measure is computed from permuting predictor values for the l-th data set:
For each tree, the prediction error on the l-th data set is recorded. Then the same is done
after permuting each predictor variable from the l-th data set. The difference between the two prediction errors are
then averaged over all trees.
}

\value{

\item{fold_importance }{ Fold-specific permutation variable importance measure. For classification the mean
                         decrease in accuracy over all classes is used, for regression the mean decrease in MSE. }
\item{type}{one of regression, classification}
}

\references{
Janitza S, Celik E, Boulesteix A-L, (2015),\emph{ A computationally fast variable
importance test for random forest for high dimensional data},Technical Report 185,
University of Munich, <\url{http://nbn-resolving.de/urn/resolver.pl?urn=nbn:de:bvb:19-epub-25587-4}>

}


\seealso{
\code{\link[randomForest]{importance}}, \code{\link[randomForest]{randomForest}}
}
\examples{
##############################
#      Classification        #
##############################
## Simulating data
X = replicate(8,rnorm(100))
X= data.frame( X) #"X" can also be a matrix
z  = with(X,5*X1 + 3*X2 + 2*X3 + 1*X4 -
          5*X5 - 9*X6 - 2*X7 + 1*X8 )
pr = 1/(1+exp(-z))         # pass through an inv-logit function
y = as.factor(rbinom(100,1,pr))

##################################################################
## Split indexes 2- folds
k = 2
cuts = round(length(y)/k)
from = (0:(k-1)*cuts)+1
to = (1:k*cuts)
rs = sample(1:length(y))
l = 1
##################################################################
## Compute fold-specific permutation variable importance
library("randomForest")

lth = rs[from[l]:to[l]]
# without the l-th data set
Xl = X[-lth,]
yl = y[-lth]
cl.rf_l = randomForest(Xl,yl,keep.forest = TRUE)
# the l-th data set
X_l = X[lth,]
y_l = y[lth]
# Compute l-th fold-specific variable importance
cvl_varim=VarImpCVl(X_l,y_l,cl.rf_l)

##############################
#      Regression            #
##############################
##################################################################
## Simulating data:
X = replicate(15,rnorm(120))
X = data.frame( X) #"X" can also be a matrix
y = with(X,2*X1 + 2*X2 + 2*X3 + 1*X4 - 2*X5 - 2*X6 - 1*X7 + 2*X8 )
##################################################################
## Split indexes 2- folds
k = 2
cuts = round(length(y)/k)
from = (0:(k-1)*cuts)+1
to = (1:k*cuts)
rs = sample(1:length(y))
l = 1

##################################################################
## Compute fold-specific permutation variable importance
library("randomForest")

lth = rs[from[l]:to[l]]
# without the l-th data set
Xl = X[-lth,]
yl = y[-lth]
reg.rf_l = randomForest(Xl,yl,keep.forest = TRUE)
# the l-th data set
X_l = X[lth,]
y_l = y[lth]
# Compute l-th fold-specific variable importance
CVVI_l = VarImpCVl(X_l,y_l,reg.rf_l)

}

