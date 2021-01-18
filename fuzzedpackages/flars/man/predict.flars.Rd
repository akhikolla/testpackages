\name{predict.flars}
\alias{predict.flars}
\title{
Prediction for functional least angle regression.
}
\description{
This is the function that carries out the prediction of the new observations. 
}
\usage{
\method{predict}{flars}(object,newdata,...)
}
\arguments{
  \item{object}{
This must be a \code{flars} object from the function \code{\link{flars}}.
}
  \item{newdata}{
A list of new observations. The format of this set of data must be the same as the training data, including the order of the variables.
}
  \item{...}{
Other arguments to input.
}
}

\value{
A matrix of predictions. Since the input \code{flars} object may have more than one estimated coefficients, the number of predictions may be more than one set. Each column of the outcome is corresponding to one set of coefficients. 
}

\examples{
library(flars)
library(fda)
## Generate some data.
dataL=data_generation(seed = 1,uncorr = TRUE,nVar = 8,nsamples = 120,
      var_type = 'm',cor_type = 3)

## Split the training data and the testing data
nTrain=80
nsamples=120

TrainIdx=seq(nTrain)
TestIdx=seq(nsamples)[-TrainIdx]
fsmTrain=lapply(dataL$x,function(fsmI) fsmI[TrainIdx,,drop=FALSE])
fsmTest=lapply(dataL$x,function(fsmI) fsmI[TestIdx,,drop=FALSE])
yTrain=dataL$y[TrainIdx]
yTest=dataL$y[TestIdx]
  
## Do the variable selection
out=flars(fsmTrain,yTrain,method='basis',max_selection=9,
    normalize='norm',lasso=FALSE)
	
## Do the prediction
pred=predict(out,newdata = fsmTest)

# apply(pred,2,flars:::rmse,yTest)
}
