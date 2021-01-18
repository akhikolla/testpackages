\encoding{UTF-8}
\name{pauclog}
\alias{pauclog}

\title{
 Calculates the p-values
}
\description{
  This auxiliary function calculates the logarithm of p-values of the statistical significance test of the difference of samples from two classes using AUC values (for each input feature). It takes as an input the results of the AUC value calculation by the function \code{\link{compute.aucs}}.
  It can be reasonably used only for two-class problem.
  The results is in the form of \dQuote{numeric vector} with the logarithms of the p-values for each features.
}
\usage{
pauclog(auc,n=100,n.plus=0.5,labels=numeric(),pos=numeric())
}
\arguments{
  \item{auc}{a numeric vector of AUC values.}
  \item{n}{ the whole number of observations for the test.}
  \item{n.plus}{the number of cases in the sample with the positive class.}
  \item{labels}{the factor with the class labels.}
  \item{pos}{the numeric vector with the level of the positive class.}
}
\details{
  This auxiliary function's main job is to calculate the logarithm of p-values of the statistical significance test of two samples, defined by negative and positive class labels, i.e. two-class problem. See the   \dQuote{Value} section to this page for more details.

}
\value{

  A returned data consists is the following:
  \item{pauclog}{a numeric vector with the logarithm of p-value for each feature}
  }

\references{
   David J. Hand and Robert J. Till (2001). A Simple Generalisation of the Area Under the ROC Curve for Multiple Class Classification Problems. Machine Learning 45(2), p. 171–186.
}

\seealso{
\code{\link{compute.aucs}}, \code{\link{pauc}}
}

\examples{
# example
data(datasetF6)

# class label must be factor
datasetF6[,ncol(datasetF6)]<-as.factor(datasetF6[,ncol(datasetF6)])

auc.val=compute.aucs(dattable=datasetF6)
vauc<-auc.val[,"AUC"]
val=levels(datasetF6[,ncol(datasetF6)])

if(length(val)==2)
{
	 pos=auc.val[,"Positive class"]
	 paucv<-pauclog(auc=vauc,labels=datasetF6[,ncol(datasetF6)],pos=pos)
}else{
	 num.size=100
	 num.prop=0.5
	 paucv<-pauclog(auc=vauc,n=num.size,n.plus=num.prop)
}
}

\keyword{feature selection}
\keyword{classification}
\keyword{AUC values}
