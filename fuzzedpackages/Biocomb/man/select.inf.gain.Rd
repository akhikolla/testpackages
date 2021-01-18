\encoding{UTF-8}
\name{select.inf.gain}
\alias{select.inf.gain}

\title{
 Ranks the features
}
\description{
  This function calculates the features weights using the Information Gain criterion measure and performs the ranking of the features (in decreasing order of Information Gain criteria). It can handle both numerical and nominal values. At first it performs the discretization of the numerical features values, according to several optional discretization methods using the function \code{\link{ProcessData}}. This function measures the worth of a feature
by computing the Information Gain criterion measure with respect to the class.The results is in the form of \dQuote{data.frame}, consisting of the following fields: features (Biomarker) names, values of the Information Gain criterion measure and the positions of the features in the dataset. The features in the data.frame are sorted according to the Information Gain uncertainty criterion values. This function is used internally to perform the classification with feature selection using the function \dQuote{classifier.loop} with argument \dQuote{InformationGain} for feature selection. The variable \dQuote{NumberFeature} of the data.frame is passed to the classification function.
}
\usage{
select.inf.gain(matrix,disc.method,attrs.nominal)
}
\arguments{
  \item{matrix}{a dataset, a matrix of feature values for several cases, the last column is for the class labels. Class labels could be numerical or character values. The maximal number of classes is ten.}
  \item{disc.method}{a method used for feature discretization.The discretization options include minimal description length (MDL), equal frequency and equal interval width methods.}
  \item{attrs.nominal}{a numerical vector, containing the column numbers of the nominal features, selected for the analysis.}
}
\details{
  This function's main job is to rank the features according to Information Gain criterion. See the   \dQuote{Value} section to this page for more details. Before
  starting it calls the \code{\link{ProcessData}} function to make the discretization of numerical features.

  Data can be provided in matrix form, where the rows correspond to cases with feature values and class label. The columns contain the values of individual features and the last column must contain class labels. The maximal number of class labels equals 10.
  The class label features and all the nominal features must be defined as factors.
}
\value{
  The data can be provided with reasonable number of missing values that must be at first preprocessed with one of the imputing methods in the function  \code{\link{input_miss}}.
  A returned list consists of the the following fields:
  \item{Biomarker}{a character vector of feature names}
  \item{Information.Gain}{a numeric vector of Information gain values for the features}
  \item{NumberFeature}{a numerical vector of the positions of the features in the dataset}
  }

\references{
   Y. Wang, I.V. Tetko, M.A. Hall, E. Frank, A. Facius, K.F.X. Mayer, and H.W. Mewes, "Gene Selection from Microarray Data for Cancer Classification—A Machine Learning Approach," Computational Biology and Chemistry, vol. 29, no. 1, pp. 37-46, 2005.
}

\seealso{
 \code{\link{ProcessData}}, \code{\link{input_miss}}, \code{\link{select.process}}
}

\examples{
# example for dataset without missing values
data(data_test)

# class label must be factor
data_test[,ncol(data_test)]<-as.factor(data_test[,ncol(data_test)])
disc<-"equal interval width"
attrs.nominal=numeric()
out=select.inf.gain(data_test,disc.method=disc,attrs.nominal=attrs.nominal)

# example for dataset with missing values
data(leukemia_miss)
xdata=leukemia_miss

# class label must be factor
xdata[,ncol(xdata)]<-as.factor(xdata[,ncol(xdata)])

# nominal features must be factors
attrs.nominal=101
xdata[,attrs.nominal]<-as.factor(xdata[,attrs.nominal])

delThre=0.2
out=input_miss(xdata,"mean.value",attrs.nominal,delThre)
if(out$flag.miss)
{
 xdata=out$data
}
disc<-"equal interval width"
out=select.inf.gain(xdata,disc.method=disc,attrs.nominal=attrs.nominal)
}

\keyword{feature selection}
\keyword{classification}
\keyword{information gain}
\keyword{missing values}
