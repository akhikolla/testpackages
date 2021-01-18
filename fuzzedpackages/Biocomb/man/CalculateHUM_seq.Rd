\encoding{UTF-8}
\name{CalculateHUM_seq}
\alias{CalculateHUM_seq}

\title{
 Calculate HUM value
}
\description{
  This function calculates the features weights using the HUM (hypervolume under manifold) values criterion measure and is used for ranking the features (in decreasing order of HUM values). HUM values are the extension of the AUC values for more than two classes.
  It can handle only numerical values.
  It computes a HUM value and returns a \dQuote{List} object, consisting of HUM value and the best permutation of class labels in \dQuote{seq} vector. This \dQuote{seq} vector can be passed to the function \code{\link{CalculateHUM_ROC}} for the calculating the coordinates of the 2D or 3D ROC.
  This function is used internally to perform the classification with feature selection using the function \dQuote{classifier.loop} with argument \dQuote{HUM} for feature selection.
}
\usage{
CalculateHUM_seq(data,indexF,indexClass,indexLabel)
}
\arguments{
  \item{data}{a dataset, a matrix of feature values for several cases, the additional column with class labels is provided. Class labels could be numerical or character values. The maximal number of classes is ten. The \code{indexClass} determines the column with class labels.}
  \item{indexF}{a numeric or character vector, containing the column numbers or column names of the analyzed features.}
  \item{indexClass}{a numeric or character value, containing the column number or column name of the class labels.}
  \item{indexLabel}{a character vector, containing the column names of the class labels, selected for the analysis.}
}
\details{
  This function's main job is to compute the maximal HUM value between the all possible permutations of class labels, selected for analysis. See the
  \dQuote{Value} section to this page for more details. Before
  returning, it will call the \code{\link{CalcGene}} function to calculate the HUM value for each feature (object).

  Data can be provided in matrix form, where the rows correspond to cases with feature values and class label. The columns contain the values of individual features and the separate column contains class labels. The maximal number of class labels equals 10. The computational efficiency of the function descrease in the case of more than 1000 cases with more than 6 class labels.
  In order to use all the functions of the package  it is necessary to put the class label in the last column of the dataset.The class label features must be defined as factors.
}
\value{
  The data can be provided with reasonable number of missing values that must be at first preprocessed with one of the imputing methods in the function  \code{\link{input_miss}}.
  A returned list consists of th the following fields:
  \item{HUM}{a list of HUM values for the specified number of analyzed features}
  \item{seq}{a list of vectors, each containing the sequence of class labels}
  }

\references{
   Li, J. and Fine, J. P. (2008): ROC Analysis with Multiple Tests and Multiple Classes: methodology and its application in microarray studies.\emph{Biostatistics}. 9 (3): 566-576.\cr
   Natalia Novoselova, Cristina Della Beffa, Junxi Wang, Jialiang Li, Frank Pessler, Frank Klawonn. HUM Calculator and HUM package for R: easy-to-use software tools for multicategory receiver operating characteristic analysis» / Bioinformatics. – 2014. – Vol. 30 (11): 1635-1636 doi:10.1093/ bioinformatics/btu086.
}

\seealso{
 \code{\link{CalculateHUM_Ex}}, \code{\link{CalculateHUM_ROC}}
}

\examples{
data(leukemia72)
# Basic example
# class label must be factor
leukemia72[,ncol(leukemia72)]<-as.factor(leukemia72[,ncol(leukemia72)])

xdata=leukemia72
indexF=1:2
indexClass=ncol(xdata)
label=levels(xdata[,indexClass])
indexLabel=label[1:2]

out=CalculateHUM_seq(xdata,indexF,indexClass,indexLabel)

}

\keyword{feature selection}
\keyword{classification}
\keyword{HUM values}
\keyword{AUC values}
\keyword{missing values}
