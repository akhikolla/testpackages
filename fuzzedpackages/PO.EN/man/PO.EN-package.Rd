\name{PO.EN-package}
\alias{PO.EN-package}
\alias{PO.EN-package}
\docType{package}
\title{
  An Elastic-Net Regularized  presence-only Model
}
\description{
 This package fits a  presence-only model with elastic-net penalty using coordinate descent. This package also provides a feature of tuning the prevalence parameter through a two-dimensional cross-validation. The package can be used in genetics study mainly for predicting regulatory effects of genetic variants given a large number of epigenetic features.
}
\details{
 Accept typical presence-only response vector \code{y}, a vector consisted of presence and background observations, and design matrix \code{x}. Three main functions:
\tabular{ll}{
\tab \cr
\code{cv.PO.EN} \tab The cross-validation tuning function\cr
  \tab \cr
\code{PO.EN} \tab The main model-fitting function\cr
  \tab \cr
\code{PO.EN.predict} \tab The predicting function\cr
}
}
\author{
Zikun Yang, Chen Wang, Iuliana Ionita-Laza

Maintainer: Zikun Yang <yangzikun1125@gmail.com>
}
\references{
  Zikun Yang, Chen Wang, Iuliana Ionita-Laza. A robust  presence-only model to predict regulatory effects of genetic variants at single nucleotide resolution by integrating epigenetic information and massively parallel reporter assays. 2020
}
\keyword{presence-only elastic-net }
\examples{
data(example.data) # example training dataset, including training dataset and testing dataset
train_data<-example.data$train.data
y_train=train_data$response;x_train=train_data[,-1]  # response and design matrix of training data
test_data<-example.data$test.data
y_test=test_data$response;x_test=test_data[,-1]  # response and design matrix of testing data
\donttest{PO.EN.cv<-cv.PO.EN(x_train,y_train,input.pi=seq(0.01,0.4,length.out=4))}
\dontshow{PO.EN.cv<-list()
          PO.EN.cv$lambda.min<-0.1
          PO.EN.cv$pi<-0.2
      }
PO.EN.beta<-PO.EN(x_train,y_train,lambda=PO.EN.cv$lambda.min,
            true.prob=PO.EN.cv$pi,beta_start=rep(0,ncol(x_train)+1))
predictions<-PO.EN.predict(x_test,PO.EN.beta)
pROC::roc(y_test~predictions)
}
