\encoding{UTF-8}
\name{chi2.algorithm}
\alias{chi2.algorithm}

\title{
 Select the subset of features
}
\description{
  This function selects the subset of features on the basis of the Chi2 discretization algorithm. The algorithm provides the way to select numerical features while discretizing them. It is based on the \eqn{\chi^2} statistic, and consists of two phases of discretization. According to the value of \eqn{\chi^2} statistic for each pair of adjacent intervals the merging of the intervals continues until an inconsistency rate is exceeded. Chi2 algorithms automatically determines a proper \eqn{\chi^2} threshold that keeps the fidelity of the original data. The nominal features must be determined as they didn't take part in the discretization process but in the process of inconsistency rate calculation. In the process of discretization the irrelevant features are removed. The results is in the form of \dQuote{list}, consisting of two fields: the processed dataset without irrelevant features and the names of the selected features. This function is used internally to perform the classification with feature selection using the function \dQuote{classifier.loop} with argument \dQuote{Chi2-algorithm} for feature selection.
}
\usage{
chi2.algorithm(matrix,attrs.nominal,threshold)
}
\arguments{
  \item{matrix}{a dataset, a matrix of feature values for several cases, the last column is for the class labels. Class labels could be numerical or character values. The maximal number of classes is ten.}
  \item{attrs.nominal}{a numerical vector, containing the column numbers of the nominal features, selected for the analysis.}
  \item{threshold}{a numeric threshold value for the inconsistency rate.}
}
\details{
  This function's main job is to select the subset of informative numerical features using the two phase process of feature values merging according to the \eqn{\chi^2} statistic for the pairs of adjacent intervals. The stopping criterion of merging is the inconsistency rate of the processed dataset. See the
  \dQuote{Value} section to this page for more details.

  Data can be provided in matrix form, where the rows correspond to cases with feature values and class label. The columns contain the values of individual features and the last column must contain class labels. The maximal number of class labels equals 10.
  The class label features and all the nominal features must be defined as factors.
}
\value{
  The data can be provided with reasonable number of missing values that must be at first preprocessed with one of the imputing methods in the function  \code{\link{input_miss}}.

  A returned data.frame consists of the the following fields:
  \item{data.out}{the processed dataset without irrelevant features (features which have been\cr merged into a single interval)}
  \item{subset}{a character vector of the selected feature names}
  }

\references{
   H. Liu and L. Yu. "Toward Integrating Feature Selection Algorithms for Classification and Clustering", IEEE Trans. on Knowledge and Data Engineering, pdf, 17(4), 491-502, 2005.
}

\seealso{
\code{\link{input_miss}}, \code{\link{select.process}}
}

\examples{
# example for dataset without missing values
#p1=Sys.time()
data(data_test)
# not all features to select
xdata=data_test[,c(1:6,ncol(data_test))]
# class label must be factor
xdata[,ncol(xdata)]<-as.factor(xdata[,ncol(xdata)])
attrs.nominal=numeric()
threshold=0
out=chi2.algorithm(matrix=xdata,attrs.nominal=attrs.nominal,threshold=threshold)
#Sys.time()-p1
}

\keyword{feature selection}
\keyword{classification}
\keyword{chi-squared}
\keyword{missing values}
