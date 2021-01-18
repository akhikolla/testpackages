\encoding{UTF-8}
\name{compute.auc.permutation}
\alias{compute.auc.permutation}

\title{
 Calculates the p-values
}
\description{
  This auxiliary function calculates the p-value of the significance of the AUC values using the permutation test (for each input feature). It takes as an input the results of the AUC value calculation using function \code{\link{compute.aucs}}.

  The results is in the form of \dQuote{numeric vector} with p-values for each AUC value.
}
\usage{
compute.auc.permutation(aucs,dattable,repetitions=1000)
}
\arguments{
  \item{aucs}{a numeric vector of AUC values.}
  \item{dattable}{a dataset, a matrix of feature values for several cases, the last column is for the class labels. Class labels could be numerical or character values.}
  \item{repetitions}{the number of permutations of feature values.}
}
\details{
  This auxiliary function's main job is to calculate the p-values of the statistical significance test of the AUC values for each input feature). See the   \dQuote{Value} section to this page for more details.

  Data can be provided in matrix form, where the rows correspond to cases with feature values and class label. The columns contain the values of individual features and the last column must contain class labels. The maximal number of class labels equals 10.
}
\value{
The data can be provided with reasonable number of missing values that must be at first preprocessed with one of the imputing methods in the function  \code{\link{input_miss}}.
  A returned data is the following:
  \item{p.values}{a numeric vector with the p-values for each feature AUC value}
  }

\references{
   David J. Hand and Robert J. Till (2001). A Simple Generalisation of the Area Under the ROC Curve for Multiple Class Classification Problems. Machine Learning 45(2), p. 171–186.
}

\seealso{
\code{\link{compute.aucs}}, \code{\link{pauclog}}, \code{\link{pauc}}, \code{\link{compute.auc.random}}
}

\examples{
# example
data(data_test)

# class label must be factor
data_test[,ncol(data_test)]<-as.factor(data_test[,ncol(data_test)])

auc.val=compute.aucs(dattable=data_test)
vauc<-auc.val[,"AUC"]
rep.num<-20

p.values=compute.auc.permutation(aucs=vauc,dattable=data_test,rep.num)
}

\keyword{feature selection}
\keyword{classification}
\keyword{AUC values}
