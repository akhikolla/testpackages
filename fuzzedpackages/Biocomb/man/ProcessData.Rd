\encoding{UTF-8}
\name{ProcessData}
\alias{ProcessData}

\title{
 Select the subset of features
}
\description{
  The auxiliary function performs the discretization of the numerical features and is called from the several functions for feature selection. The discretization options include minimal description length (MDL), equal frequency and equal interval width methods.
  The results is in the form of \dQuote{list}, consisting of two fields: the processed dataset and the column numbers of the features. When the value of the input parameter \dQuote{flag}=TRUE the second field will include the column numbers of the features, which have more than single interval after discretization.
}
\usage{
ProcessData(matrix,disc.method,attrs.nominal,flag=FALSE)
}
\arguments{
  \item{matrix}{a dataset, a matrix of feature values for several cases, the last column is for the class labels. Class labels could be numerical or character values. The maximal number of classes is ten.}
  \item{disc.method}{a method used for feature discretization.The discretization options include minimal description length (MDL), equal frequency and equal interval width methods.}
  \item{attrs.nominal}{a numerical vector, containing the column numbers of the nominal features, selected for the analysis.}
  \item{flag}{a binary logical value. If equals TRUE the output list will contain the processed dataset with the features, having more than one interval after discretization together with their names. In the case of FALSE value the processed dataset with all the features will be returned.}
}
\details{
  This auxiliary function's main job is to descritize the numerical features using the one of the discretization methods. See the
  \dQuote{Value} section to this page for more details.


  Data can be provided in matrix form, where the rows correspond to cases with feature values and class label. The columns contain the values of individual features and the last column must contain class labels. The maximal number of class labels equals 10.
  The class label features and all the nominal features must be defined as factors.
}
\value{
  The data can be provided with reasonable number of missing values that must be at first preprocessed with one of the imputing methods in the function  \code{\link{input_miss}}.

  A returned list consists of the the following fields:
  \item{m3}{a processed dataset}
  \item{sel.feature}{a numeric vector with the column numbers of the features, having more than one interval value (when \dQuote{flag}=TRUE). If \dQuote{flag}=FALSE it return all the column numbers of the dataset.}
  }

\references{
   H. Liu, F. Hussain, C. L. Tan, and M. Dash, "Discretization: An enabling technique," Data Mining and Knowledge Discovery, Vol. 6, No. 4, 2002, pp. 393-423.
}

\seealso{
\code{\link{select.inf.gain}}, \code{\link{select.inf.symm}}, \code{\link{select.inf.chi2}}, \cr \code{\link{select.fast.filter}}, \code{\link{select.process}}
}

\examples{
# example for dataset without missing values
data(data_test)

# class label must be factor
data_test[,ncol(data_test)]<-as.factor(data_test[,ncol(data_test)])

disc<-"MDL"
attrs.nominal=numeric()
flag=FALSE
out=ProcessData(matrix=data_test,disc.method=disc,
attrs.nominal=attrs.nominal,flag=flag)
}

\keyword{feature selection}
\keyword{classification}
\keyword{discretization}
