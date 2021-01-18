\encoding{UTF-8}
\name{generate.data.miss}
\alias{generate.data.miss}

\title{
 Generate the dataset with missing values
}
\description{
  The function for the generation the dataset with missing values from the input dataset with all the values. It is mainly intended for the testing purposes.
  The results is in the form of \dQuote{data.frame} which corresponds to the input data.frame or matrix, where missing values are inserted. The percent of missing values is supplied as the input parameters. The processed dataset can be used in the algorithms for missing value imputation \dQuote{input_miss} or for any other purposes.
}
\usage{
generate.data.miss(data,percent=5,filename=NULL)
}
\arguments{
  \item{data}{a dataset, a matrix of feature values for several cases, the last column is for the class labels. Class labels could be numerical or character values. This data set has not missing values}
  \item{percent}{a numerical value for the percent of the missing values to be inserted into the dataset.}
  \item{filename}{a character name of the output file to save the dataset with missing values.}
}
\details{
  This function's main job is to generate the dataset with missing values from the input dataset with all the values. See the \dQuote{Value} section to this page for more details.

  Data can be provided in matrix form, where the rows correspond to cases with feature values and class label. The columns contain the values of individual features and the last column must contain class labels. The maximal number of class labels equals 10.
  The class label features and all the nominal features must be defined as factors.
}
\value{
  A returned data.frame corresponds to the input dataset with inserted missing values.
  }

\references{
 McShane LM, Radmacher MD, Freidlin B, Yu R, Li MC, Simon R. Methods for assessing reproducibility of clustering patterns observed in analyses of microarray data. Bioinformatics. 2002 Nov;18(11):1462-9.
}

\seealso{
 \code{\link{input_miss}}, \code{\link{select.process}},
 \code{\link{classifier.loop}}
}

\examples{
# example

data(leukemia72_2)

percent =5
f.name=NULL #file name to include
out=generate.data.miss(data=leukemia72_2,percent=percent,filename=f.name)
}

\keyword{feature selection}
\keyword{classification}
\keyword{missing values}
