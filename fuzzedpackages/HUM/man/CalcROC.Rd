\encoding{UTF-8}
\name{CalcROC}
\alias{CalcROC}

\title{
 Calculate ROC points
}
\description{
  This is the auxiliary function of the HUM package. It computes a point coordinates for plotting ROC curve and returns a \dQuote{List} object,  consisting of sensitivity and specificity values for 2D-ROC curve and 3D-points for 3D-ROC curve, the optimal threshold values with the corresponding feature values and the accuracy of the classifier (feature).
}
\usage{
CalcROC(s_data, seq, thresholds)
}

\arguments{
  \item{s_data}{a list, which contains the vectors of sorted feature values for individual class labels.
  }
  \item{seq}{a numeric vector, containing the particular permutation of class labels.
  }
  \item{thresholds}{a numeric vector, containing the values of thresholds for calculating ROC curve coordinates.}
}

\details{
  This function's main job is to compute the point coordinates to plot the 2D- or 3D-ROC curve, the optimal threshold values and the accuracy of classifier. See the \dQuote{Value} section to this page for more details. The optimal threshold for two-class problem is the pair of sensitivity and specificity values for the selected feature. The optimal threshold for three-class problem is the 3D-point with the coordinates presenting the fraction of the correctly classified data objects for each class. The calculation of the optimal threshold is described in section \dQuote{Threshold}.
}

\value{
  The data must be provided without missing values in order to process. A returned list consists of the following fields: 
  \item{Sn}{a specificity values for 2D-ROC construction and the first coordinate for 3D-ROC construction}
  
  \item{Sp}{a sensitivity values for 2D-ROC construction and the second coordinate for 3D-ROC construction}
  \item{S3}{the third coordinate for 3D-ROC construction}
  \item{optSn}{the optimal specificity value for 2D-ROC construction and the first coordinate of the op-timal threshold for 3D-ROC construction}
  
  \item{optSp}{the optimal sensitivity value for 2D-ROC construction and the second coordinate of the optimal threshold for 3D-ROC construction}
  \item{optS3}{the third coordinate of the optimal threshold for 3D-ROC construction
  }
  \item{optThre}{the feature value according to the optimal threshold (optSn,optSp) for two-class problem}
  \item{optThre1}{the first feature value according to the optimal threshold (optSn,optSp,optS3) for three-class problem}
  \item{optThre2}{the second feature value according to the optimal threshold (optSn,optSp,optS3) for three-class problem}
  \item{accuracy}{the accuracy of classifier (feature) for the optimal threshold}
  }
  
\section{Threshold}{
The optimal threshold value is calculated for two-class problem as the pair \dQuote{(optSn, optSp)} corresponding to the maximal value of \dQuote{Sn+Sp}. According to \dQuote{(optSn, optSp)} the corresponding feature threshold value \dQuote{optThre} is calculated.
The optimal threshold value is calculated for three-class problem as the pair \dQuote{(optSn, optSp,optS3)} corresponding to the maximal value of \dQuote{Sn+Sp+S3}.According to \dQuote{(optSn, optSp,optS3)} the corresponding feature threshold values \dQuote{optThre1,optThre2} are calculated.
The accuracy of the classifier is the mean value of dQuote{(optSn, optSp)} for two-class problem and the mean value of \dQuote{(optSn, optSp,optS3)} for three-class problem.
}

\section{Errors}{
  If there exists NA values for features or class labels no HUM value can be calculated and an error is triggered with
  message \dQuote{Values are missing}.
}

\references{
  Li, J. and Fine, J. P. (2008): ROC Analysis with Multiple Tests and Multiple Classes: methodology and its application in microarray studies.\emph{Biostatistics}. 9 (3): 566-576.
}

\seealso{
 \code{\link{CalculateHUM_Ex}}, \code{\link{CalculateHUM_ROC}}
}

\examples{
data(sim)
indexF=names(sim[,c(3,4),drop = FALSE])
indexClass=2
label=unique(sim[,indexClass])
indexLabel=label[1:3]
out=CalculateHUM_seq(sim,indexF,indexClass,indexLabel)
HUM<-out$HUM
seq<-out$seq

indexL=NULL
for(i in 1:length(indexLabel))
{
  indexL=c(indexL,which(label==indexLabel[i]))
}
  
indexEach=NULL
indexUnion=NULL

for(i in 1:length(label))
{
  vrem=which(sim[,indexClass]==label[i])
  indexEach=c(indexEach,list(vrem))
  if(length(intersect(label[i],indexLabel))==1)
    indexUnion=union(indexUnion,vrem)
}
s_data=NULL
dataV=sim[,indexF[1]]  #single feature
prodValue=1
for (j in 1:length(indexLabel))
{
  vrem=sort(dataV[indexEach[[indexL[j]]]])

  s_data=c(s_data,list(vrem))
  prodValue = prodValue*length(vrem)
}
#calculate the threshold values for plot of 2D ROC and 3D ROC
thresholds <- sort(unique(dataV[indexUnion]))
thresholds=(c(-Inf, thresholds) + c(thresholds, +Inf))/2
  
out=CalcROC(s_data,seq[,indexF[1]], thresholds)
}

\keyword{multivariate}
\keyword{nonparametric}
\keyword{plot}
\keyword{roc}
\keyword{auc}
