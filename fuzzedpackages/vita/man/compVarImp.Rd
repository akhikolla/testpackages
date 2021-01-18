\name{compVarImp}
\alias{compVarImp}

\title{
Compute permutation variable importance measure
}
\description{
Compute permutation variable importance measure from a random forest for classification and regression.
}
\usage{
compVarImp(X, y,rForest,nPerm=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{a data frame or a matrix of predictors.}
  \item{y}{a response vector. If a factor, classification is assumed, otherwise regression is assumed.}
  \item{rForest}{an object of class \code{\link[randomForest]{randomForest}}, keep.forest,keep.inbag must
                be set to True.}
  \item{nPerm}{Number of times the OOB data are permuted per tree for assessing variable importance.
              Number larger than 1 gives slightly more stable estimate, but not very effective.
              Currently only implemented for regression. }
}
\details{
The permutation variable importance measure is computed from permuting OOB data: For each tree,
the prediction error on the out-of-bag observations is recorded. Then the same is done
after permuting a predictor variable. The differences between the two error rates are then averaged over all
trees.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use

\item{importance }{The permutation variable importance measure. A matrix with nclass + 1
                  (for classification) or one (for regression) columns. For classification, the first
                  nclass columns are the class-specific measures computed as mean decrease in accuracy.
                  The nclass + 1st column is the mean decrease in accuracy over all classes. For regression
                  the mean decrease in MSE is given. }
\item{importanceSD }{The "standard errors" of the permutation-based importance measure. For classification, a p
                     by nclass + 1 matrix corresponding to the first nclass + 1 columns of the importance matrix.
                    For regression a vector of length p.}
\item{type}{one of regression, classification}

%% ...
}
\references{
 Breiman L. (2001), \emph{Random Forests}, Machine Learning 45(1),5-32, <doi:10.1023/A:101093340432>
}

\seealso{
\code{\link[randomForest]{importance}}, \code{\link[randomForest]{randomForest}},\code{\link[vita]{CVPVI}}
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
##############################
## Classification with Random Forest:
library("randomForest")
cl.rf= randomForest(X,y,mtry = 3,ntree=100,
                    importance=TRUE,keep.inbag = TRUE)

##############################
## Permutation variable importance measure
vari= compVarImp(X,y,cl.rf)

##############################
#compare them with the original results
cbind(cl.rf$importance[,1:3],vari$importance)
cbind(cl.rf$importance[,3],vari$importance[,3])
cbind(cl.rf$importanceSD,vari$importanceSD)
cbind(cl.rf$importanceSD[,3],vari$importanceSD[,3])
cbind(cl.rf$type,vari$type)


###############################
#      Regression             #
###############################
## Simulating data
X = replicate(8,rnorm(100))
X= data.frame( X) #"X" can also be a matrix
y= with(X,5*X1 + 3*X2 + 2*X3 + 1*X4 -
          5*X5 - 9*X6 - 2*X7 + 1*X8 )
##############################
## Regression with Random Forest:
library("randomForest")
reg.rf= randomForest(X,y,mtry = 3,ntree=100,
                     importance=TRUE,keep.inbag = TRUE)

##############################
## Permutation variable importance measure
vari= compVarImp(X,y,reg.rf)

##############################
#compare them with the original results
cbind(importance(reg.rf, type=1, scale=FALSE),vari$importance)
cbind(reg.rf$importanceSD,vari$importanceSD)
cbind(reg.rf$type,vari$type)

}
