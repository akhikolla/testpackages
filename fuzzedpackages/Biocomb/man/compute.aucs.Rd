\encoding{UTF-8}
\name{compute.aucs}
\alias{compute.aucs}

\title{
 Ranks the features
}
\description{
  This function calculates the features weights using the AUC (Area Under the ROC Curve) values.
  It can handle only numerical values. This function performs two-class or multiclass AUC. A multiclass AUC is a mean of AUCs for all combinations of the two class labels.  This function measures the worth of a feature
by computing the AUC values with respect to the class.The results is in the form of \dQuote{data.frame}. In the case of two-class problem it consists of the three fields: features (Biomarker) names, AUC values and level of the positive class. In the case of more than two classes it consists of two fields: features (Biomarker) names, AUC values. This function is used internally to perform the classification with feature selection using the function \dQuote{classifier.loop} with argument \dQuote{auc} for feature selection.
}
\usage{
compute.aucs(dattable)
}
\arguments{
  \item{dattable}{a dataset, a matrix of feature values for several cases, the last column is for the class labels. Class labels could be numerical or character values. The maximal number of classes is ten.}
}
\details{
  This function's main job is to calculate the weights of the features according to AUC values. See the   \dQuote{Value} section to this page for more details.

  Data can be provided in matrix form, where the rows correspond to cases with feature values and class label. The columns contain the values of individual features and the last column must contain class labels. The maximal number of class labels equals 10.
  The class label features must be defined as factors.
}
\value{
  The data can be provided with reasonable number of missing values that must be at first preprocessed with one of the imputing methods in the function  \code{\link{input_miss}}.
  A returned data.frame consists of th the following fields:
  \item{Biomarker}{a character vector of feature names}
  \item{AUC}{a numeric vector of AUC values for the features according to class}
  \item{Positive class}{a numeric vector of positive class levels for two-class problem}
  }

\references{
   David J. Hand and Robert J. Till (2001). A Simple Generalisation of the Area Under the ROC Curve for Multiple Class Classification Problems. Machine Learning 45(2), p. 171–186.
}

\seealso{
\code{\link{input_miss}}, \code{\link{select.process}}
}

\examples{
# example for dataset without missing values
data(data_test)

# class label must be factor
data_test[,ncol(data_test)]<-as.factor(data_test[,ncol(data_test)])

out=compute.aucs(dattable=data_test)

# example for dataset with missing values
data(leukemia_miss)
xdata=leukemia_miss

# class label must be factor
xdata[,ncol(xdata)]<-as.factor(xdata[,ncol(xdata)])

# the nominal features must be factors
attrs.nominal=101
xdata[,attrs.nominal]<-as.factor(xdata[,attrs.nominal])

delThre=0.2
out=input_miss(xdata,"mean.value",attrs.nominal,delThre)
if(out$flag.miss)
{
 xdata=out$data
}
xdata=xdata[,-attrs.nominal]
# the nominal features are not processed
out=compute.aucs(dattable=xdata)
}

\keyword{feature selection}
\keyword{classification}
\keyword{AUC values}
\keyword{missing values}
