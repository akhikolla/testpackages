\encoding{UTF-8}
\name{select.process}
\alias{select.process}

\title{
 Feature ranking and feature selection
}
\description{
  The main function for the feature ranking or feature subset selection. It can handle both numerical and nominal values.
  It presents the infrastructure to perform the feature ranking or feature selection for the data set with two or more class labels. The function calls several feature ranking methods with different quality measures, including AUC values (functions \code{\link{compute.aucs}}), information gain (function \code{\link{select.inf.gain}}), symmetrical uncertainty (function \code{\link{select.inf.symm}}), chi-squared (\eqn{\chi^2}) statistic (function \code{\link{select.inf.chi2}}). It also calls the number of feature selection methods, including fast correlation-based filter method (FCBF) (function \code{\link{select.fast.filter}}), Chi2 discretization algorithm (function \code{\link{chi2.algorithm}}), CFS algorithm with forward search (function \code{\link{select.forward.Corr}}), wrapper method with decision tree algorithm and forward search strategy (function \code{\link{select.forward.wrapper}}).
  The results is in the form of \dQuote{numeric vector} with the column numbers of the selected features for features selection algorithms and ordered features' column numbers according to the criteria for feature ranking. The number of features can be limited to the \dQuote{max.no.features} , which is the function input parameter. The output of the function is used in function \dQuote{classifier.loop} in the process of classification.
}
\usage{
select.process(dattable,method="InformationGain",disc.method="MDL",
threshold=0.2,threshold.consis=0.05,attrs.nominal=numeric(),
max.no.features=10)
}
\arguments{
  \item{dattable}{a dataset, a matrix of feature values for several cases, the last column is for the class labels. Class labels could be numerical or character values. The maximal number of classes is ten.}
  \item{method}{a method of feature ranking or feature selection. There are 6 methods for feature ranking ("auc", "HUM", "Chi-square", "InformationGain", "symmetrical.uncertainty", "Relief") and 4 methods for feature selection ("FastFilter", "CFS", "CorrSF", "Chi2-algorithm")}
  \item{disc.method}{a method used for feature discretization. There are three options "MDL","equal interval width","equal frequency". The discretization options "MDL" assigned to the minimal description length (MDL) discretization algorithm, which is a supervised algorithm. The last two options refer to the unsupervized discretization algorithms.}
  \item{threshold}{a numeric threshold value for the correlation of feature with class to be included in the final subset. It is used by fast correlation-based filter method (FCBF)}
  \item{threshold.consis}{a numeric threshold value for the inconsistency rate. It is used by Chi2 discretization algorithm.}
  \item{attrs.nominal}{a numerical vector, containing the column numbers of the nominal features, selected for the analysis.}
  \item{max.no.features}{the maximal number of features to be selected or ranked.}
}
\details{
  This function's main job is to present the infrostructure to perform the feature ranking or feature selection for the data set with two or more class labels. See the
  \dQuote{Value} section to this page for more details.

  Data can be provided in matrix form, where the rows correspond to cases with feature values and class label. The columns contain the values of individual features and the last column must contain class labels. The maximal number of class labels equals 10.
  The class label features and all the nominal features must be defined as factors.
}
\value{
  The data can be provided with reasonable number of missing values that must be at first preprocessed with one of the imputing methods in the function  \code{\link{input_miss}}.

  A returned value is
  \item{sel.feat}{a vector of column numbers of the selected features for features selection algorithms and ordered features' column numbers according to the criteria for feature ranking}
  }

\references{
 H. Liu and L. Yu. "Toward Integrating Feature Selection Algorithms for Classification and Clustering", IEEE Trans. on Knowledge and Data Engineering, pdf, 17(4), 491-502, 2005.\cr
   L. Yu and H. Liu. "Feature Selection for High-Dimensional Data: A Fast Correlation-Based Filter Solution". In Proceedings of The Twentieth International Conference on Machine Leaning (ICML-03), Washington, D.C. pp. 856-863. August 21-24, 2003.
}

\seealso{
 \code{\link{select.inf.gain}}, \code{\link{select.inf.symm}}, \code{\link{select.inf.chi2}}, \cr \code{\link{select.fast.filter}}, \code{\link{chi2.algorithm}}, \code{\link{select.forward.Corr}}, \cr \code{\link{select.forward.wrapper}}, \code{\link{input_miss}}
}

\examples{
# example for dataset without missing values
data(data_test)

# class label must be factor
data_test[,ncol(data_test)]<-as.factor(data_test[,ncol(data_test)])

method="InformationGain"
disc<-"MDL"
thr=0.1
thr.cons=0.05
attrs.nominal=numeric()
max.f=15

out=select.process(data_test,method=method,disc.method=disc,
threshold=thr, threshold.consis=thr.cons,attrs.nominal=attrs.nominal,
max.no.features=max.f)


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

method="InformationGain"
disc<-"MDL"
thr=0.1
thr.cons=0.05
max.f=15

out=select.process(xdata,method=method,disc.method=disc,
threshold=thr, threshold.consis=thr.cons,attrs.nominal=attrs.nominal,
max.no.features=max.f)
}

\keyword{feature selection}
\keyword{classification}
\keyword{information gain}
\keyword{symmetrical uncertainty}
\keyword{chi-squared}
\keyword{missing values}
