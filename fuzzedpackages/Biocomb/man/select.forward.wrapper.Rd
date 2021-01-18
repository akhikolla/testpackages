\encoding{UTF-8}
\name{select.forward.wrapper}
\alias{select.forward.wrapper}

\title{
 Select the subset of features
}
\description{
  This function selects the subset of features using the wrapper method with decision tree algorithm and forward search strategy. It can handle both numerical and nominal values. The wrapper method makes use of the classification algorithm in order to estimate the quality measure of the feature subset. The method uses the built-in cross-validation procedure to estimate the accuracy of classification for the feature subset. At the first step of the method the one-feature subset is selected according to the quality measure. In the following steps the subset is incrementally extended according to the forward search strategy until the stopping criterion is met.
  The result is in the form of character vector with the names of the selected features. This function is used internally to perform the classification with feature selection using the function \dQuote{classifier.loop} with argument \dQuote{CorrSF} for feature selection.
}
\usage{
select.forward.wrapper(dattable)
}
\arguments{
  \item{dattable}{a dataset, a matrix of feature values for several cases, the last column is for the class labels. Class labels could be numerical or character values. The maximal number of classes is ten.}

}
\details{
  This function's main job is to select the subset of informative features according to forward selection strategy using the wrapper method. The decision tree is used as the classifier to estimate the quality of the feature subset. See the
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

out=select.forward.wrapper(dattable=data_test)
}

\keyword{feature selection}
\keyword{classification}
\keyword{information gain}
\keyword{missing values}
