\encoding{UTF-8}
\name{CalculateHUM_ROC}
\alias{CalculateHUM_ROC}

\title{
 Compute the points for ROC curve
}
\description{
  This is the  function for computing the points for ROC curve. It returns a \dQuote{List} object, consisting of sensitivity and specificity values for 2D-ROC curve and 3D-points for 3D-ROC curve. Also the optimal threshold values are returned.
  It can handle only numerical values.
}
\usage{
CalculateHUM_ROC(data,indexF,indexClass,indexLabel,seq)
}

\arguments{
  \item{data}{a dataset, a matrix of feature values for several cases, the additional column with class labels is provided. Class labels could be numerical or character values. The maximal number of classes is ten. The \code{indexClass} determines the column with class labels.
  }
  \item{indexF}{a numeric or character vector, containing the column numbers or column names of the analyzed features.
  }
  \item{indexClass}{a numeric or character value, containing the column number or column name of the class labels.
  }
  \item{indexLabel}{a character vector, containing the column names of the class labels, selected for the analysis.}
  \item{seq}{a numeric matrix, containing the permutation of the class labels for all features.}
}

\details{
  This function's main job is to compute the point coordinates to plot the 2D- or 3D-ROC curve and the optimal threshold values. See the \dQuote{Value} section to this page for more details. The function calls the \code{\link{CalcROC}} to calculate the point coordinates, optimal thresholds and accuracy of classifier (feature) in the threshold.

  Data can be provided in matrix form, where the rows correspond to cases with feature values and class label. The columns contain the values of individual features and the separate column contains class labels. The maximal number of class labels equals 10.
  In order to use all the functions of the package  it is necessary to put the class label in the last column of the dataset. The class label features must be defined as factors.
}

\value{
  The data must be provided without missing values in order to process. A returned list consists of th the following fields:
  \item{Sn}{a specificity values for 2D-ROC construction and the first coordinate for 3D-ROC construction}

  \item{Sp}{a sensitivity values for 2D-ROC construction and the second coordinate for 3D-ROC construction}
  \item{S3}{the third coordinate for 3D-ROC construction}
   \item{optSn}{the optimal specificity value for 2D-ROC construction and the first coordinate of the op-timal threshold for 3D-ROC construction}

  \item{optSp}{the optimal sensitivity value for 2D-ROC construction and the second coordinate of the optimal threshold for 3D-ROC construction}
  \item{optS3}{the third coordinate of the optimal threshold for 3D-ROC construction}
  }


\section{Errors}{
  If there exists NA values for features or class labels no HUM value can be calculated and an error is triggered with
  message \dQuote{Values are missing}.
}

\references{
  Li, J. and Fine, J. P. (2008): ROC Analysis with Multiple Tests and Multiple Classes: methodology and its application in microarray studies.\emph{Biostatistics}. 9 (3): 566-576.\cr
   Natalia Novoselova, Cristina Della Beffa, Junxi Wang, Jialiang Li, Frank Pessler, Frank Klawonn. HUM Calculator and HUM package for R: easy-to-use software tools for multicategory receiver operating characteristic analysis» / Bioinformatics. – 2014. – Vol. 30 (11): 1635-1636 doi:10.1093/ bioinformatics/btu086.
}

\seealso{
 \code{\link{CalculateHUM_Ex}}, \code{\link{CalculateHUM_seq}}
}

\examples{
data(leukemia72)
# Basic example
# class label must be factor
leukemia72[,ncol(leukemia72)]<-as.factor(leukemia72[,ncol(leukemia72)])

xdata=leukemia72
indexF=1:3
indexClass=ncol(xdata)
label=levels(xdata[,indexClass])
indexLabel=label[1:2]

out=CalculateHUM_seq(xdata,indexF,indexClass,indexLabel)
HUM<-out$HUM
seq<-out$seq
out=CalculateHUM_ROC(xdata,indexF,indexClass,indexLabel,seq)
}

\keyword{feature selection}
\keyword{classification}
\keyword{HUM values}
\keyword{AUC values}
\keyword{missing values}
