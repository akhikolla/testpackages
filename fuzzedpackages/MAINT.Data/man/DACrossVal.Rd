\name{DACrossVal}

\title{Cross Validation for Discriminant Analysis Classification Rules}

\alias{DACrossVal}

\description{
\sQuote{DACrossVal} evaluates the performance of a Discriminant Analysis training sample algorithm by k-fold 
Cross-Validation.  
}

\usage{DACrossVal(data, grouping, TrainAlg, EvalAlg=EvalClrule, 
Strfolds=TRUE, kfold=10, CVrep=20, prior="proportions", loo=FALSE, \dots)}

\arguments{

  \item{data}{Matrix, data frame or Interval Data object of observations.}

  \item{grouping}{Factor specifying the class for each observation.}

  \item{TrainAlg}{A function with the training algorithm. It should return an object that can be used as input to the argument of \sQuote{EValAlg}.}

  \item{EvalAlg}{A function with the evaluation algorithm. By default set to \sQuote{EvalClrule} which returns a list with components \dQuote{err} (estimates of error rates by class) and \dQuote{Nk} (number of out-sample observations by class). This default can be used for all \sQuote{TrainAlg} 
arguments that return an object with a predict method returning a list with a \sQuote{class} component (a factor) containing the classification results.}

  \item{Strfolds}{Boolean flag indicating if the folds should be stratified according to the original class
proportions (default), or randomly generated from the whole training sample, ignoring class membership.}

  \item{kfold}{Number of training sample folds to be created in each replication.}

  \item{CVrep}{Number of replications to be performed.}

  \item{prior}{The prior probabilities of class membership.  If unspecified, the class proportions for the training set are used.  If present, the probabilities should be specified in the order of the factor levels.}

  \item{loo}{A boolean flag indicating if a leave-one-out strategy should be employed. When set to \dQuote{TRUE} overrides the kfold and CVrep arguments.}

  \item{\dots}{Further arguments to be passed to \sQuote{TrainAlg} and \sQuote{EvalAlg}.} 
}

\value{A three dimensional array with the number of tested observations, and estimated classification errors for each combination of fold and replication tried. The array dimensions are defined as follows:\cr
The first dimension runs through the different fold-replication combinations.\cr
The second dimension represents the classes. \cr
The third dimension has two named levels representing respectively the number of observations tested (\dQuote{Nk}), and the estimated classification errors (\dQuote{Clerr}).}

\author{A. Pedro Duarte Silva}

\seealso{\code{\link{lda}}, \code{\link{qda}}, \code{\link{IData}}}

\examples{

\dontrun{

# Compare performance of linear and quadratic discriminant analysis with 
#  Covariance cases C1 and c4 on the ChinaT data set by 5-fold cross-validation 
#  replicated twice

# Create an Interval-Data object containing the intervals for 899 observations 
# on the temperatures by quarter in 60 Chinese meteorological stations.

ChinaT <- IData(ChinaTemp[1:8])

# Classical (configuration 1) Linear Discriminant Analysis 

CVldaC1 <- DACrossVal(ChinaT,ChinaTemp$GeoReg,TrainAlg=lda,CovCase=1,kfold=5,CVrep=2)
summary(CVldaC1[,,"Clerr"])

# Linear Discriminant Analysis with covariance case 3

CVldaC4 <- DACrossVal(ChinaT,ChinaTemp$GeoReg,TrainAlg=lda,CovCase=3,kfold=5,CVrep=2)
summary(CVldaC4[,,"Clerr"])

# Classical (configuration 1) Quadratic Discriminant Analysis 

CVqdaC1 <- DACrossVal(ChinaT,ChinaTemp$GeoReg,TrainAlg=qda,CovCase=1,kfold=5,CVrep=2)
summary(CVqdaC1[,,"Clerr"])

# Quadratic Discriminant Analysis with covariance case 3

CVqdaC4 <- DACrossVal(ChinaT,ChinaTemp$GeoReg,TrainAlg=qda,CovCase=3,kfold=5,CVrep=2)
summary(CVqdaC4[,,"Clerr"])

}

}

\keyword{discriminant analysis}
\keyword{cross-validation}

