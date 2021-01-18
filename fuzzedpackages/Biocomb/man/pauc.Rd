\encoding{UTF-8}
\name{pauc}
\alias{pauc}

\title{
 Calculates the p-values
}
\description{
  This auxiliary function calculates the p-value of the statistical significance test of the difference of samples from two classes using AUC values (for each input feature). It takes as an input the results of the AUC value calculation using function \code{\link{compute.aucs}}.
  It can be reasonably used only for two-class problem.
  The results is in the form of \dQuote{numeric vector} with p-values for each features.
}
\usage{
pauc(auc,n=100,n.plus=0.5,labels=numeric(),pos=numeric())
}
\arguments{
  \item{auc}{a numeric vector of AUC values.}
  \item{n}{ the whole number of observations for the test.}
  \item{n.plus}{the number of cases in the sample with the positive class.}
  \item{labels}{the factor with the class labels.}
  \item{pos}{the numeric vector with the level of the positive class.}
}
\details{
  This auxiliary function's main job is to calculate the p-values of the statistical significance test of two samples, defined by negative and positive class labels, i.e. two-class problem. See the   \dQuote{Value} section to this page for more details.

}
\value{

  A returned data consists is the following:
  \item{pauc}{a numeric vector with the p-values for each feature}
  }

\references{
   David J. Hand and Robert J. Till (2001). A Simple Generalisation of the Area Under the ROC Curve for Multiple Class Classification Problems. Machine Learning 45(2), p. 171–186.
}

\seealso{
\code{\link{compute.aucs}}, \code{\link{pauclog}}
}

\examples{
# example
data(data_test)

# class label must be factor
data_test[,ncol(data_test)]<-as.factor(data_test[,ncol(data_test)])

auc.val=compute.aucs(dattable=data_test)
vauc<-auc.val[,"AUC"]
val=levels(data_test[,ncol(data_test)])

if(length(val)==2)
{
	 pos=auc.val[,"Positive class"]
	 paucv<-pauc(auc=vauc,labels=data_test[,ncol(data_test)],pos=pos)
}else{
	 num.size=100
	 num.prop=0.5
	 paucv<-pauc(auc=vauc,n=num.size,n.plus=num.prop)
}
}

\keyword{feature selection}
\keyword{classification}
\keyword{AUC values}
