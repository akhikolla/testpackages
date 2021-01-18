\encoding{UTF-8}
\name{select.forward.Corr}
\alias{select.forward.Corr}

\title{
 Select the subset of features
}
\description{
  This function selects the subset of features using the forward search strategy on the basis of correlation measure (CFS algorithm with forward search). CFS evaluates a subset of features by considering the individual predictive ability of each feature along with the degree of redundancy between them. It can handle both numerical and nominal values. At the beginning the discretization of the numerical features values is performed using the function \code{\link{ProcessData}}. At the first step of the method the one-feature subset is selected according to its informative score, which takes into account the average feature to class correlation and the average feature to feature correlation. In the following steps the subset is incrementally extended according to the forward search strategy until the stopping criterion is met.
  The result is in the form of character vector with the names of the selected features. This function is used internally to perform the classification with feature selection using the function \dQuote{classifier.loop}.
}
\usage{
select.forward.Corr(matrix,disc.method,attrs.nominal)
}
\arguments{
  \item{matrix}{a dataset, a matrix of feature values for several cases, the last column is for the class labels. Class labels could be numerical or character values. The maximal number of classes is ten.}
  \item{disc.method}{a method used for feature discretization.The discretization options include minimal description length (MDL), equal frequency and equal interval width methods.}
  \item{attrs.nominal}{a numerical vector, containing the column numbers of the nominal features, selected for the analysis.}

}
\details{
  This function's main job is to select the subset of informative features according to forward selection strategy using the correlation measure (informative theoretic measure). The measure consideres the individual predictive ability of each feature along with the degree of redundancy between them. See the
  \dQuote{Value} section to this page for more details.


  Data can be provided in matrix form, where the rows correspond to cases with feature values and class label. The columns contain the values of individual features and the last column must contain class labels. The maximal number of class labels equals 10.
  The class label features and all the nominal features must be defined as factors.
}
\value{
  The data can be provided with reasonable number of missing values that must be at first preprocessed with one of the imputing methods in the function  \code{\link{input_miss}}.

  A returned value is
  \item{subset}{a character vector of the names of selected features}
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
disc<-"MDL"
attrs.nominal=numeric()
out=select.forward.Corr(matrix=data_test,disc.method=disc,
attrs.nominal=attrs.nominal)
}

\keyword{feature selection}
\keyword{classification}
\keyword{correlation}
\keyword{missing values}
