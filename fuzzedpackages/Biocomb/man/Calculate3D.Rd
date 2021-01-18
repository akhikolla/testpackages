\encoding{UTF-8}
\name{Calculate3D}
\alias{Calculate3D}

\title{
 Plot the 3D-ROC curve
}
\description{
  This function plots the 3D-ROC curve using the point coordinates, computed by the function \code{\link{CalculateHUM_ROC}}. Optionally visualizes the optimal threshold point, which gives the maximal accuracy of the classifier(feature) (see \code{\link{CalcROC}}).
}
\usage{
Calculate3D(sel,Sn,Sp,S3,optSn,optSp,optS3,thresholds,HUM,
name,print.optim=TRUE)
}

\arguments{
  \item{sel}{a character value, which is the name of the selected feature.
  }
   \item{Sn}{a numeric vector of the x-coordinates of the ROC curve..
  }
  \item{Sp}{a numeric vector of the y-coordinates of the ROC curve.
  }
  \item{S3}{a numeric vector of the z-coordinates of the ROC curve.
  }
  \item{optSn}{the first coordinate of the optimal threshold
  }
  \item{optSp}{the second coordinate of the optimal threshold
  }
  \item{optS3}{the third coordinate of the optimal threshold
  }
  \item{thresholds}{a numeric vector with threshold values to calculate point coordinates.
  }
  \item{HUM}{a numeric vector of HUM values, calculated using function.
}
  \item{name}{a character vector of class labels.
  }
  \item{print.optim}{a boolean parameter to plot the optimal threshold point on the graph. The default value is TRUE.
}
}

\details{
  This function's main job is to plot the 3D-ROC curve according to the given point coordinates.
}

\value{
  The function doesn't return any value.
  }


\section{Errors}{
   If there exists NA values for specificity or sensitivity values, or HUM values the plotting fails and an error is triggered with message \dQuote{Values are missing}
}

\references{
   Li, J. and Fine, J. P. (2008): ROC Analysis with Multiple Tests and Multiple Classes: methodology and its application in microarray studies.\emph{Biostatistics}. 9 (3): 566-576.\cr
   Natalia Novoselova, Cristina Della Beffa, Junxi Wang, Jialiang Li, Frank Pessler, Frank Klawonn. HUM Calculator and HUM package for R: easy-to-use software tools for multicategory receiver operating characteristic analysis» / Bioinformatics. – 2014. – Vol. 30 (11): 1635-1636 doi:10.1093/ bioinformatics/btu086.
}

\seealso{
 \code{\link{CalculateHUM_seq}}, \code{\link{CalculateHUM_ROC}}
}

\examples{
data(leukemia72)
# Basic example
# class label must be factor
leukemia72[,ncol(leukemia72)]<-as.factor(leukemia72[,ncol(leukemia72)])

xdata=leukemia72

indexF=names(xdata)[10]

indexClass=ncol(xdata)
label=levels(xdata[,indexClass])
indexLabel=label

out=CalculateHUM_seq(xdata,indexF,indexClass,indexLabel)
HUM<-out$HUM
seq<-out$seq
out=CalculateHUM_ROC(xdata,indexF,indexClass,indexLabel,seq)
Calculate3D(indexF,out$Sn,out$Sp,out$S3,out$optSn,out$optSp,out$optS3,
out$thresholds,HUM,indexLabel[seq])
}

\keyword{feature selection}
\keyword{classification}
\keyword{HUM values}
\keyword{AUC values}
\keyword{missing values}

