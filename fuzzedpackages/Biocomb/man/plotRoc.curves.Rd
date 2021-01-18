\encoding{UTF-8}
\name{plotRoc.curves}
\alias{plotRoc.curves}

\title{
 Plots the ROC curve for two-class problem
}
\description{
 This function plots the ROC curve for the two-class problem.
}
\usage{
plotRoc.curves(dattable,file.name=NULL,colours=NULL,ltys=NULL,
 add.legend=F,curve.names=NULL,include.auc=F,xaxis="",yaxis="",
 line.width=2,headline="",ispercent=F)
}
\arguments{
  \item{dattable}{a dataset, a matrix of feature values for several cases, the last column is for the class labels. Class labels could be numerical or character values.}
  \item{file.name}{the file name to save the plot.}
  \item{colours}{the color values for each plot if more than one feature or one color value in the case of one feature.}
  \item{ltys}{the  line type values for each plot if more than one feature or one  line type in the case of one feature.}
  \item{add.legend}{logical value; if TRUE the legend will be plotted at the bottom righ.}
  \item{curve.names}{a character value or vector in the case of more than one feature with curve names to be used in the legend.}
  \item{include.auc}{logical value; if TRUE the AUC value will be included in the legend .}
  \item{xaxis}{character value with the name of X axis.}
  \item{yaxis}{character value with the name of Y axis.}
  \item{line.width}{a positive number for line width.}
  \item{headline}{the character value with the name of the plot.}
  \item{ispercent}{logical value; if TRUE the true positive and false positive values of the plot are in the percents.}
}
\details{
  This function's main job is to plot the ROC curve for one or more features with the possibility to include the AUC values in the legend.

  Data can be provided in matrix form, where the rows correspond to cases with feature values and class label. The columns contain the values of individual features and the last column must contain class labels with two class labels.
}
\value{
  The data can be provided with reasonable number of missing values that must be at first preprocessed with one of the imputing methods in the function  \code{\link{input_miss}}.
 }

\references{
   David J. Hand and Robert J. Till (2001). A Simple Generalisation of the Area Under the ROC Curve for Multiple Class Classification Problems. Machine Learning 45(2), p. 171–186.
}

\seealso{
\code{\link{compute.aucs}}, \code{\link{pauclog}}, \code{\link{pauc}}, \code{\link{compute.auc.permutation}}, \code{\link{input_miss}}
}

\examples{
# example for dataset without missing values
data(leukemia72_2)

# class label must be factor
leukemia72_2[,ncol(leukemia72_2)]<-as.factor(leukemia72_2[,ncol(leukemia72_2)])

add.legend<-TRUE
include.auc<-TRUE

attrs.no=c(1,2)
xdata=leukemia72_2[,c(attrs.no,ncol(leukemia72_2))]
plotRoc.curves(dattable=xdata,add.legend=add.legend,include.auc=include.auc)
}

\keyword{feature selection}
\keyword{classification}
\keyword{AUC values}
\keyword{ROC curve}
