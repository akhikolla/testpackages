\encoding{UTF-8}
\name{compute.auc.random}
\alias{compute.auc.random}

\title{
 Calculates the p-values
}
\description{
  This auxiliary function calculates the p-value of the significance of the AUC values using the comparison with random sample generation (for each input feature). It takes as an input the results of the AUC value calculation using function \code{\link{compute.aucs}}.

  The results is in the form of \dQuote{numeric vector} with p-values for each AUC value.
}
\usage{
compute.auc.random(aucs,dattable,repetitions=10000,correction="none")
}
\arguments{
  \item{aucs}{a numeric vector of AUC values.}
  \item{dattable}{a dataset, a matrix of feature values for several cases, the last column is for the class labels. Class labels could be numerical or character values.}
  \item{repetitions}{the number of repetitions of random sample' generation.}
  \item{correction}{the method of p-value correction for multiple testing, including Bonferroni-Holm, Bonferroni corrections or without correction.}
}
\details{
  This auxiliary function's main job is to calculate the p-values of the statistical significance test of the AUC values for each input feature for two-class problem.

  Data can be provided in matrix form, where the rows correspond to cases with feature values and class label. The columns contain the values of individual features and the last column must contain class labels (with two class labels).

  The correction methods include the Bonferroni correction ("bonferroni") in which the p-values are multiplied by the number of comparisons and the less conservative corrections by Bonferroni-Holm method ("bonferroniholm"). A pass-through option ("none") is also included.

The correction methods are designed to give strong control of the family-wise error rate.
See the   \dQuote{Value} section to this page for more details.
}
\value{
The data can be provided with reasonable number of missing values that must be at first preprocessed with one of the imputing methods in the function  \code{\link{input_miss}}.
  A returned data is the following:
  \item{pvalues.raw}{a numeric vector with the corrected p-values for each feature AUC value}
  }

\references{
   Benjamini, Y., and Yekutieli, D. (2001). The control of the false discovery rate in multiple testing under dependency. Annals of Statistics 29, 1165–1188.
}

\seealso{
\code{\link{compute.aucs}}, \code{\link{pauclog}}, \code{\link{pauc}}, \code{\link{compute.auc.permutation}}
}

\examples{
# example
data(datasetF6)

# class label must be factor
datasetF6[,ncol(datasetF6)]<-as.factor(datasetF6[,ncol(datasetF6)])

auc.val=compute.aucs(dattable=datasetF6)
vauc<-auc.val[,"AUC"]

cors<-"none"
rep.num<-100

pvalues.raw<-compute.auc.random(aucs=vauc,dattable=datasetF6,
 repetitions=rep.num,correction=cors)
}

\keyword{feature selection}
\keyword{classification}
\keyword{AUC values}
