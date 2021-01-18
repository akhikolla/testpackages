\encoding{UTF-8}
\name{CalculateHUM_Ex}
\alias{CalculateHUM_Ex}

\title{
 Calculate HUM value
}
\description{
  This is the main function of the HUM package. It computes a HUM value and returns a \dQuote{List} object, consisting of HUM value and the best permutation of class labels in \dQuote{seq} vector. This \dQuote{seq} vector can be passed to the function \code{\link{CalculateHUM_ROC}}. 
}
\usage{
CalculateHUM_Ex(data,indexF,indexClass,allLabel,amountL)
}

\arguments{
  \item{data}{a dataset, a matrix of feature values for several cases, the additional column with class labels is provided. Class labels could be numerical or character values. The maximal number of classes is ten. The \code{indexClass} determines the column with class labels.
  }
  \item{indexF}{a numeric or character vector, containing the column numbers or column names of the analyzed features.
  }
  \item{indexClass}{a numeric or character value, containing the column number or column name of the class labels.
  }
  \item{allLabel}{a character vector, containing the column names of the class labels, selected for the analysis.}
  \item{amountL}{a character vector, containing the column names of the class labels, selected for the analysis.}
}

\details{
  This function's main job is to compute the maximal HUM value between the all possible permutations of class labels, selected for analysis. See the
  \dQuote{Value} section to this page for more details. Before
  returning, it will call the \code{\link{CalcGene}} function to calculate the HUM value for each case (object).

  Data can be provided in matrix form, where the rows correspond to cases with feature values and class label. The columns contain the values of individual features and the separate column contains class labels. The maximal number of class labels equals 10. The computational efficiency of the function descrease in the case of more than 1000 cases with more than 6 class labels..
}

\value{
  The data must be provided without missing values in order to process. A returned list consists of th the following fields: 
  \item{HUM}{a list of HUM values for the specified number of analyzed features}
  
  \item{seq}{a list of vectors, each containing the sequence of class labels}
  }


\section{Errors}{
  If there exists NA values for features or class labels no HUM value can be calculated and an error is triggered with
  message \dQuote{Values are missing}.
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
indexF=c(3,4);
indexClass=2;
allLabel=c("Normal","OrthArthr","OA","Early")
amountL=2
out=CalculateHUM_Ex(sim,indexF,indexClass,allLabel,amountL)
}

\keyword{multivariate}
\keyword{nonparametric}
\keyword{plot}
\keyword{roc}
\keyword{auc}
