\encoding{UTF-8}
\name{cost.curve}
\alias{cost.curve}

\title{
 Plots the RCC curve for two-class problem
}
\description{
  This function plots the Relative Cost Curves (RCC) and calculates the corresponding Area Above the RCC (AAC) value to estimate the classifier performance under unequal misclassification costs.
  It is intended for the two-class problem, but the extension to more than two classes will be produced later.
  RCC is a graphical technique for visualising the performance of binary classifiers over the full range of possible relative misclassification costs. This curve provides helpful information to choose the best set of classifiers or to estimate misclassification costs if those are not known precisely. Area Above the RCC (AAC) is a scalar measure of classifier performance under unequal misclassification costs problem.
It can be reasonably used only for two-class problem.
}
\usage{
cost.curve(data, attrs.no, pos.Class, AAC=TRUE, n=101, add=FALSE,
xlab="log2(c)",ylab="relative costs", main="RCC",lwd=2,col="black",
xlim=c(-4,4), ylim=(c(20,120)))
}
\arguments{
  \item{data}{a dataset, a matrix of feature values for several cases, the last column is for the class labels. Class labels could be numerical or character values. The function is provided for two classes.}
  \item{attrs.no}{a numerical value, containing the column number of the features to construct the RCC.}
  \item{pos.Class}{a level of the class factor to be selected as the positive class for the construction of the RCC cost curve.}
  \item{AAC}{logical value; if TRUE the AAC value will be calculated.}
  \item{n}{the  number of points for the construction of RCC curve (it corresponds to the number of cost values).}
  \item{add}{logical value; if TRUE the RCC curve can be added to the existent RCC plot.}
  \item{xlab}{name of the X axis.}
  \item{ylab}{name of the Y axis.}
  \item{main}{name of the RCC plot.}
  \item{lwd}{a positive number for line width.}
  \item{col}{the color value for the RCC plot.}
  \item{xlim}{the vector with two numeric values for the X axis limits, it defines the values for log2(cost).}
  \item{ylim}{the vector with two numeric values for the Y axis limits, it defines the values for the relative cost value.}
}
\details{
  This function's main job is to plt the RCC curve and calculates the corresponding AAC value to estimate the classifier performance under unequal misclassification costs.

  Data can be provided in matrix form, where the rows correspond to cases with feature values and class label. The columns contain the values of individual features and the last column must contain class labels with two class labels.
}
\value{
  The data can be provided with reasonable number of missing values that must be at first preprocessed with one of the imputing methods in the function  \code{\link{input_miss}}.
 }

\references{
   Olga Montvida and Frank Klawonn Relative cost curves: An alternative to AUC and an extension to 3-class problems,Kybernetika 50 no. 5, 647-660, 2014
}

\seealso{
\code{\link{compute.aucs}}, \code{\link{compute.auc.random}},
\code{\link{compute.auc.permutation}},\cr \code{\link{plotRoc.curves}}, \code{\link{input_miss}}
}

\examples{
# example for dataset without missing values
data(data_test)

# class label must be factor
data_test[,ncol(data_test)]<-as.factor(data_test[,ncol(data_test)])

xllim<--4
xulim<-4
yllim<-30
yulim<-110

attrs.no=c(1,9)
pos.Class<-levels(data_test[,ncol(data_test)])[1]
add.legend<-TRUE

aacs<-rep(0,length(attrs.no))
color<-c(1:length(attrs.no))

aacs[1] <- cost.curve(data_test, attrs.no[1], pos.Class,col=color[1],add=FALSE,
 xlim=c(xllim,xulim),ylim=c(yllim,yulim))

if(length(attrs.no)>1){
			for(i in 2:length(attrs.no)){
		      aacs[i]<- cost.curve(data_test, attrs.no[i], pos.Class,
		      col=color[i],add=TRUE,xlim=c(xllim,xulim))
		    }
	    }

if(add.legend){
	       legt <- colnames(data_test)[attrs.no]
		   for(i in 1:length(attrs.no)){
			   legt[i] <- paste(legt[i],", AAC=",round(1000*aacs[i])/1000,sep="")
}
legend("bottomright",legend=legt,col=color,lwd=2)
}
}

\keyword{feature selection}
\keyword{classification}
\keyword{AUC values}
\keyword{ROC curve}
