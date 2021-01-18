\encoding{UTF-8}
\name{CalculateHUM_Plot}
\alias{CalculateHUM_Plot}

\title{
 Plot 2D-ROC curve
}
\description{
  This is the main function of the HUM package. It plots the 2D-ROC curve using the point coordinates, computed by the function \code{\link{CalculateHUM_ROC}}.Optionally visualizes the optimal threshold point, which gives the maximal accuracy of the classifier(feature) (see \code{\link{CalcROC}}). 
}
\usage{
CalculateHUM_Plot(sel,Sn,Sp,optSn,optSp,HUM,print.optim=TRUE)
}

\arguments{
  \item{sel}{a character value, which is the name of the selected feature.
  }
   \item{Sn}{a numeric vector of the x-coordinates of the ROC curve, which is the specificity values of the standard ROC analysis..
  }
  \item{Sp}{a numeric vector of the y-coordinates of the ROC curve, which is the sensitivity values of the standard ROC analysis..
  }
  \item{optSn}{the optimal specificity value for 2D-ROC construction
  }
  \item{optSp}{the optimal sensitivity value for 2D-ROC construction
  }
  \item{HUM}{a numeric vector of HUM values, calculated using function \code{\link{CalculateHUM_seq}} or \code{\link{CalculateHUM_Ex}}.
  }
  \item{print.optim}{a boolean parameter to plot the optimal threshold point on the graph. The default value is TRUE. 
}
}

\details{
  This function's main job is to plot the 2D-ROC curve according to the given point coordinates.
}

\value{
  The function doesn't return any value.
  }


\section{Errors}{
  If there exists NA values for specificity or sensitivity values, or HUM values the plotting fails and an error is triggered with message \dQuote{Values are missing}.
}

\references{
  Li, J. and Fine, J. P. (2008): ROC Analysis with Multiple Tests and Multiple Classes: methodology and its application in microarray studies.\emph{Biostatistics}. 9 (3): 566-576.
}

\seealso{
 \code{\link{CalculateHUM_seq}}, \code{\link{CalculateHUM_ROC}}
}

\examples{
data(sim)
# Basic example
indexF=names(sim[,c(3),drop = FALSE])
indexClass=2
label=unique(sim[,indexClass])
indexLabel=label[1:2]
out=CalculateHUM_seq(sim,indexF,indexClass,indexLabel)
HUM<-out$HUM
seq<-out$seq
out=CalculateHUM_ROC(sim,indexF,indexClass,indexLabel,seq)
CalculateHUM_Plot(indexF,out$Sn,out$Sp,out$optSn,out$optSp,HUM)
}

\keyword{multivariate}
\keyword{nonparametric}
\keyword{plot}
\keyword{roc}
\keyword{auc}
