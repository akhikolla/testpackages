\name{HUM-package}
\alias{HUM-package}
\alias{HUM}
\docType{package}
\title{
HUM calculation
}
\description{
  Functions to calculate AUC (area under curve) value  for two classes and HUM (hypervolume under manifold) for more class labels in order to estimate the informativity of features to outcome. Tools for visualizing ROC curve in 2D- and 3D-space.
}
\details{
\tabular{ll}{
Package: \tab HUM\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-10-25\cr
License: \tab GPL (>= 3)\cr
}
The basic unit of the HUM package is the \code{\link{CalculateHUM_seq}} function. It will calculate the AUC in case of two class labels and HUM for more than two class labels for the input features. Function \code{\link{CalculateHUM_Ex}} is the extension of main function and provides the possibility to calculate all the combinations of \code{amountL} from all the class labels. Function  \code{\link{CalculateHUM_ROC}} calculates the point coordinates in order to plot the 2D- and 3D-ROC curve, accuracy and the optimal threshold for the classifier (feature). The Functions \code{\link{CalcGene}} and \code{\link{CalcROC}} are the auxiliar function to perform the calculation. Function \code{\link{CalcROC}} calculates the point coordinates of a single feature for two-class or three-class problem, the optimal threshold for the 2-D and 3-D ROC curve and the corresponding feature values, the accuracy of the classifier (feature) for the optimal threshold.
}
\section{Functions}{
  \tabular{ll}{
    \code{\link{CalculateHUM_seq}} \tab Calculate a maximal HUM value amd the corresponding permutation of class labels\cr
    \code{\link{CalculateHUM_Ex}} \tab Calculate the HUM values with exaustive serach for specified number of class labels \cr
    \code{\link{CalculateHUM_ROC}} \tab Function to construct and plot the 2D- or 3d-ROC curve \cr
    \code{\link{CalcGene}} \tab Compute the HUM value for one feature \cr
    \code{\link{CalcROC}} \tab Compute the point coordinates to plot the 2D- or 3D-ROC curve \cr
    \code{\link{CalculateHUM_Plot}} \tab Plot the 2D-ROC curve \cr
    \code{\link{Calculate3D}} \tab Plot the 3D-ROC curve \cr
  }
}

\section{Dataset}{
  This package comes with one simulated dataset and a real dataset of 92 patients with 11 features with disease.
}

\section{Installing and using}{
  To install this package, make sure you are connected to the internet and issue the following command in the R prompt:
  
  \preformatted{
    install.packages("HUM")
  }
  
  To load the package in R:
 
  \preformatted{
    library(HUM)
  }
}

\author{
Natalia Novoselova, Frank Pessler

Maintainer: Natalia Novoselova <novos65@mail.ru>
}
\references{
  Li, J. and Fine, J. P. (2008): ROC Analysis with Multiple Tests and Multiple Classes: methodology and its application in microarray studies.\emph{Biostatistics}. 9 (3): 566-576.
}

\keyword{package}
\keyword{multivariate}
\keyword{nonparametric}
\keyword{plot}
\keyword{roc}
\keyword{auc}

\seealso{
CRAN packages \pkg{pROC}, or Bioconductor's \pkg{roc}
for ROC curves.

CRAN packages \pkg{Rcpp}, \pkg{gtools}, \pkg{rgl} employed in this package.
}
\examples{
data(sim)

# Compute the HUM value with all possible class label permutation
indexF=c(3,4);
indexClass=2;
label=unique(sim[,indexClass])
indexLabel=label[1:3]
out=CalculateHUM_seq(sim,indexF,indexClass,indexLabel)
# Compute the HUM value with exaustive search of all class label combinations
\dontrun{data(sim)
indexF=c(3,4);
indexClass=2;
labels=unique(sim[,indexClass])
amountL=4;
out=CalculateHUM_Ex(sim,indexF,indexClass,labels,amountL)
}
# Calculate the coordinates for 2D- or 3D- ROC curve and the optimal threshold point
\dontrun{data(sim)
indexF=names(sim[,c(3),drop = FALSE])
indexClass=2
label=unique(sim[,indexClass])
indexLabel=label[1:3]
out=CalculateHUM_seq(sim,indexF,indexClass,indexLabel)
HUM<-out$HUM
seq<-out$seq
out=CalculateHUM_ROC(sim,indexF,indexClass,indexLabel,seq)
}
}
