\encoding{UTF-8}
\name{select.fast.filter}
\alias{select.fast.filter}

\title{
 Select the subset of features
}
\description{
  This function selects the subset of features on the basis of the fast correlation-based filter method (FCBF). It can handle both numerical and nominal values. At first it performs the discretization of the numerical features values, according to several optional discretization methods using the function \code{\link{ProcessData}}. A fast filter can identify relevant features as well as redundancy among relevant features without pairwise correlation analysis. The overall complexity of FCBF is O(MN logN), where M - number of samples, N - number of features.The results is in the form of \dQuote{data.frame}, consisting of the features (Biomarker) names, values of the information gain and the positions of the features in the dataset. The information gain value is the correlation between the features and the class. This function is used internally to perform the classification with feature selection using the function \dQuote{classifier.loop} with argument \dQuote{FastFilter} for feature selection. The variable \dQuote{NumberFeature} of the data.frame is passed to the classification function.
}
\usage{
select.fast.filter(matrix,disc.method,threshold,attrs.nominal)
}
\arguments{
  \item{matrix}{a dataset, a matrix of feature values for several cases, the last column is for the class labels. Class labels could be numerical or character values. The maximal number of classes is ten.}
  \item{disc.method}{a method used for feature discretization.The discretization options include minimal description length (MDL), equal frequency and equal interval width methods.}
  \item{threshold}{a numeric threshold value for the correlation of feature with class to be included in the final subset.}
  \item{attrs.nominal}{a numerical vector, containing the column numbers of the nominal features, selected for the analysis.}
}
\details{
  This function's main job is to select the subset of informative features according to correlation between features and class, and between features themselves. See the
  \dQuote{Value} section to this page for more details. Before
  starting it calls the \code{\link{ProcessData}} function to make the discretization of numerical features.

  Data can be provided in matrix form, where the rows correspond to cases with feature values and class label. The columns contain the values of individual features and the last column must contain class labels. The maximal number of class labels equals 10.
  The class label features and all the nominal features must be defined as factors.
}
\value{
  The data can be provided with reasonable number of missing values that must be at first preprocessed with one of the imputing methods in the function  \code{\link{input_miss}}.

  A returned data.frame consists of the the following fields:
  \item{Biomarker}{a character vector of feature names}
  \item{Information.Gain}{a numeric vector of information gain values for the features according to class}
  \item{NumberFeature}{a numerical vector of the positions of the features in the dataset}
  }

\references{
   L. Yu and H. Liu. "Feature Selection for High-Dimensional Data: A Fast Correlation-Based Filter Solution". In Proceedings of The Twentieth International Conference on Machine Leaning (ICML-03), Washington, D.C. pp. 856-863. August 21-24, 2003.
}

\seealso{
 \code{\link{ProcessData}}, \code{\link{input_miss}}, \code{\link{select.process}}
}

\examples{
# example for dataset without missing values
data(data_test)

# class label must be factor
data_test[,ncol(data_test)]<-as.factor(data_test[,ncol(data_test)])
disc<-"MDL"
threshold=0.2
attrs.nominal=numeric()
out=select.fast.filter(data_test, disc.method=disc, threshold=threshold,
attrs.nominal=attrs.nominal)

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
disc<-"MDL"
threshold=0.2
out=select.fast.filter(xdata, disc.method=disc, threshold=threshold,
attrs.nominal=attrs.nominal)
}

\keyword{feature selection}
\keyword{classification}
\keyword{information gain}
\keyword{missing values}
