\name{flars_TrainTest}
\alias{flars_TrainTest}
\title{
Internal function for doing simulation using functional lars.
}
\description{
This is a function built for doing data generation and variable selection using functional lars with different settings and data with different correlation structures.  
}
\usage{
flars_TrainTest(seed=1,nsamples=120,nTrain=80,var_type=c('f','m'),
                VarThreshold0=0.1,SignThreshold0=0.8,cor_type=1:5,
                lasso=TRUE, check = 1,uncorr=T,nVar=8,Discrete_Norm_ID=1:12,
                NoRaw_max=12,raw_max=9,hyper=NULL,RealX=NULL,RealY=NULL,
                dataL=NULL,nCor=0,control=list())
}
\arguments{
  \item{seed}{
Set the seed for random numbers.
}
  \item{nsamples}{
Sample size of the data to generate.
}
  \item{nTrain}{
Sample size of the training data.
}
  \item{var_type}{
Two choices of the variable types. See details for more information.
}
\item{cor_type}{
Correlation structures. See details for more information.
}
\item{VarThreshold0}{
Threshold for removing variables based on variation explained. See \code{\link{flars}} for more details.
}
\item{SignThreshold0}{
Same as \code{VarThreshold0}
}
  \item{lasso}{
Use lasso modification or not. In other words, can variables selected in the former iterations be removed in the later iterations. 
}
  \item{check}{
Type of lasso check. 1 means variance check, 2 means sign check. \code{check=1} is much better than the other one. 
}
  \item{uncorr}{
If the variables are uncorrelated or not. See details for more information.
}
  \item{nVar}{
Number of variables to generate.
}
  \item{Discrete_Norm_ID}{
Which discrete method and which norm to use. 1 to 12.
}
  \item{NoRaw_max}{
Number of variables to select when not using RDP discretising method.
}
  \item{raw_max}{
Number of variables to select when using RDP discretising method.
}
  \item{hyper}{
Hyper parameters used in the Gaussian process. GP is used for building the covariance structure of the functional variables.
}
  \item{RealX}{
Real data input X.
}
  \item{RealY}{
Real data input Y.
}
  \item{dataL}{
Real input data list rather than generate in the function. It should has the same structure as that generated.
}
  \item{nCor}{
Number of cores to use.
}
  \item{control}{
List of control items. See \code{\link{fccaGen}} for more details.
}
}

\value{
A list of results using different normalization methods and different representation methods for the functional coefficients. 
}