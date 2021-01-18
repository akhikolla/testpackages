\encoding{UTF-8}
\name{select.cfs}
\alias{select.cfs}

\title{
 Select the subset of features
}
\description{
  This function selects the subset of features using the best first search strategy on the basis of correlation measure (CFS). CFS evaluates a subset of features by considering the individual predictive ability of each feature along with the degree of redundancy between them. It can handle both numerical and nominal values. The results is in the form of \dQuote{data.frame}, consisting of the following fields: features (Biomarker) names and the positions of the features in the dataset.
  This function is used internally to perform the classification with feature selection using the function \dQuote{classifier.loop} with argument \dQuote{CFS} for feature selection. The variable \dQuote{Index} of the data.frame is passed to the classification function.
}
\usage{
select.cfs(matrix)
}
\arguments{
  \item{matrix}{a dataset, a matrix of feature values for several cases, the last column is for the class labels. Class labels could be numerical or character values. The maximal number of classes is ten.}
}
\details{
  This function's main job is to select the subset of informative features according to best first search strategy using the correlation measure (informative theoretic measure). The measure consideres the individual predictive ability of each feature along with the degree of redundancy between them. See the
  \dQuote{Value} section to this page for more details.


  Data can be provided in matrix form, where the rows correspond to cases with feature values and class label. The columns contain the values of individual features and the last column must contain class labels. The maximal number of class labels equals 10.
  The class label features and all the nominal features must be defined as factors.
}
\value{
  The data can be provided with reasonable number of missing values that must be at first preprocessed with one of the imputing methods in the function  \code{\link{input_miss}}.
  A returned list consists of the the following fields:
  \item{Biomarker}{a character vector of feature names}
  \item{Index}{a numerical vector of the positions of the features in the dataset}
  }

\references{
   Y. Wang, I.V. Tetko, M.A. Hall, E. Frank, A. Facius, K.F.X. Mayer, and H.W. Mewes, "Gene Selection from Microarray Data for Cancer Classification—A Machine Learning Approach," Computational Biology and Chemistry, vol. 29, no. 1, pp. 37-46, 2005.
}

\seealso{
\code{\link{input_miss}}, \code{\link{select.process}}
}

\examples{
# example for dataset without missing values
data(data_test)

# class label must be factor
data_test[,ncol(data_test)]<-as.factor(data_test[,ncol(data_test)])
out=select.cfs(matrix=data_test)
}

\keyword{feature selection}
\keyword{classification}
\keyword{correlation}
\keyword{missing values}
