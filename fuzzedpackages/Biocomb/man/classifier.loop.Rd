\encoding{UTF-8}
\name{classifier.loop}
\alias{classifier.loop}

\title{
 Classification and classifier validation
}
\description{
  The main function for the classification and classifier validation. It performs the classification using different classification algorithms (classifiers) with the embedded feature selection and using the different schemes for the performance validation.\cr
  It presents the infrostructure to perform classification of the data set with two or more class labels. The function calls several classification methods, including Nearest shrunken centroid ("nsc"), Naive Bayes classifier ("nbc"), Nearest Neighbour classifier ("nn"), Multinomial Logistic Regression ("mlr"), Support Vector Machines ("svm"), Linear Discriminant Analysis ("lda"), Random Forest ("rf").
  The function calls the \code{\link{select.process}} in order to perform feature selection for classification, which helps to improve the quality of the classifier.
  The classifier accuracy is estimated using the embedded validation procedures, including the Repeated random sub-sampling validation ("sub-sampling"), k-fold cross-validation ("fold-crossval") and Leave-one-out cross-validation ("leaveOneOut").

  The results is in the form of \dQuote{list} with the data.frame of classification results for each selected classifier \dQuote{predictions}, matrix with the statistics for the frequency each feature is selected \dQuote{no.selected}, vector or matrix with the number of true classifications for each selected classifier \dQuote{true.classified}.
}
\usage{
classifier.loop(dattable,classifiers=c("svm","lda","rf","nsc"),
 feature.selection=c("auc","InformationGain"),
 disc.method="MDL",threshold=0.3, threshold.consis=0,
 attrs.nominal=numeric(), no.feat=20,flag.feature=TRUE,
 method.cross=c("leaveOneOut","sub-sampling","fold-crossval"))

}
\arguments{
  \item{dattable}{a dataset, a matrix of feature values for several cases, the last column is for the class labels. Class labels could be numerical or character values. The maximal number of classes is ten.}
  \item{classifiers}{the names of the classifiers.}
  \item{feature.selection}{a method of feature ranking or feature selection.}
  \item{disc.method}{a method used for feature discretization. There are three options "MDL","equal interval width","equal frequency". The discretization options "MDL" assigned to the minimal description length (MDL) discretization algorithm, which is a supervised algorithm. The last two options refer to the unsupervized discretization algorithms.}
  \item{threshold}{a numeric threshold value for the correlation of feature with class to be included in the final subset. It is used by fast correlation-based filter method (FCBF)}
  \item{threshold.consis}{a numeric threshold value for the inconsistency rate. It is used by Chi2 discretization algorithm.}
  \item{attrs.nominal}{a numerical vector, containing the column numbers of the nominal features, selected for the analysis.}
  \item{no.feat}{the maximal number of features to be selected.}
  \item{flag.feature}{logical value; if TRUE the process of classifier construction and validation will be repeated for each subset of features, starting with one feature and upwards.}
  \item{method.cross}{a character value with the names of the model validation technique for assessing how the classification results will generalize to an independent data set. It includes Repeated random sub-sampling validation, k-fold cross-validation and Leave-one-out cross-validation.}
}
\details{
  This function's main job is to perform classification with feature selection and the estimation of classification results with the model validation techniques. See the
  \dQuote{Value} section to this page for more details.

  Data can be provided in matrix form, where the rows correspond to cases with feature values and class label. The columns contain the values of individual features and the last column must contain class labels. The maximal number of class labels equals 10.
  The class label features and all the nominal features must be defined as factors.
}
\value{
  The data can be provided with reasonable number of missing values that must be at first preprocessed with one of the imputing methods in the function  \code{\link{input_miss}}.

  A returned list consists of the the following fields:
  \item{predictions}{a data.frame of classification results for each selected classifier}
  \item{no.selected}{a matrix with the statistics for each feature selection frequency}
  \item{true.classified}{a vector or matrix with the number of true classifications for each selected classifier}
  }

\references{
 S. Dudoit, J. Fridlyand, and T. P. Speed. Comparison of discrimination methods for the
classification of tumors using gene expression data. Journal of the American Statistical Association, 97(457):77–87, 2002.
}

\seealso{
 \code{\link{select.process}}, \code{\link{input_miss}}
}

\examples{
# example for dataset without missing values
data(leukemia72_2)

# class label must be factor
leukemia72_2[,ncol(leukemia72_2)]<-as.factor(leukemia72_2[,ncol(leukemia72_2)])

class.method="svm"
method="InformationGain"
disc<-"MDL"
cross.method<-"fold-crossval"

thr=0.1
thr.cons=0.05
attrs.nominal=numeric()
max.f=10

out=classifier.loop(leukemia72_2,classifiers=class.method,
feature.selection=method,disc.method=disc,
threshold=thr, threshold.consis=thr.cons,attrs.nominal=attrs.nominal,
no.feat=max.f,flag.feature=FALSE,method.cross=cross.method)


# example for dataset with missing values
\dontrun{
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

class.method="svm"
method="InformationGain"
disc<-"MDL"
cross.method<-"fold-crossval"

thr=0.1
thr.cons=0.05
max.f=10

out=classifier.loop(xdata,classifiers=class.method,
feature.selection=method,disc.method=disc,
threshold=thr, threshold.consis=thr.cons,attrs.nominal=attrs.nominal,
no.feat=max.f,flag.feature=FALSE,method.cross=cross.method)
}
}

\keyword{feature selection}
\keyword{classification}
\keyword{missing values}
