\encoding{UTF-8}
\name{plotClass.result}
\alias{plotClass.result}

\title{
 Plots the results of classifier validation schemes
}
\description{
  This function plots the barplots and boxplots, which help in estimation of the results of classifiers' validation, performed by different validation models. It must be called after the performing the classification validation with function \code{\link{classifier.loop}}.
}
\usage{
plotClass.result(true.classified, cross.method, class.method,
flag.feature, feat.num)
}
\arguments{
  \item{true.classified}{a vector or matrix of classification results for one or several classifiers and one or several feature sets. The matrix is the output value of the function \code{\link{classifier.loop}}.}
  \item{cross.method}{a character value with the names of the model validation technique for assessing how the classification results will generalize to an independent data set. It includes Repeated random sub-sampling validation, k-fold cross-validation and Leave-one-out cross-validation.}
  \item{class.method}{the names of the classifiers.}
  \item{flag.feature}{logical value; if TRUE the process of classifier construction and validation will be repeated for each subset of features, starting with one feature and upwards.}
  \item{feat.num}{the maximal number of features to be selected.}
}
\details{
  This function's main job is to plot the barplots and boxplots to visually estimate the results of classifiers' validation.

}
\value{
  The results is visualization of the plot .
 }

\references{
   S. Dudoit, J. Fridlyand, and T. P. Speed. Comparison of discrimination methods for the
classification of tumors using gene expression data. Journal of the American Statistical Association, 97(457):77–87, 2002.
}

\seealso{
 \code{\link{select.process}}, \code{\link{classifier.loop}}
}

\examples{
# example for dataset without missing values

data(leukemia72_2)

# class label must be factor
leukemia72_2[,ncol(leukemia72_2)]<-as.factor(leukemia72_2[,ncol(leukemia72_2)])

class.method=c("svm","nn")
method="InformationGain"
disc<-"MDL"
cross.method<-"fold-crossval"

flag.feature=TRUE
thr=0.1
thr.cons=0.05
attrs.nominal=numeric()
max.f=10

out=classifier.loop(leukemia72_2,classifiers=class.method,
 feature.selection=method,disc.method=disc,
 threshold=thr, threshold.consis=thr.cons,attrs.nominal=attrs.nominal,
  no.feat=max.f,flag.feature=flag.feature,method.cross=cross.method)

plotClass.result(out$true.classified, cross.method, class.method, flag.feature, max.f)

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

class.method=c("svm","nn")
method="InformationGain"
disc<-"MDL"
cross.method<-"fold-crossval"

flag.feature=TRUE
thr=0.1
thr.cons=0.05
max.f=10

out=classifier.loop(xdata,classifiers=class.method,
 feature.selection=method,disc.method=disc,
 threshold=thr, threshold.consis=thr.cons,attrs.nominal=attrs.nominal,
 no.feat=max.f,flag.feature=flag.feature,method.cross=cross.method)

plotClass.result(out$true.classified, cross.method, class.method, flag.feature, max.f)
}
}

\keyword{feature selection}
\keyword{classification}
\keyword{missing values}
