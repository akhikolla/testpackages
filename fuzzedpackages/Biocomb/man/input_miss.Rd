\encoding{UTF-8}
\name{input_miss}
\alias{input_miss}

\title{
 Process the dataset with missing values
}
\description{
  The main function for handling with missing values.
  It performs the missing values imputation using two different approachs: imputation with mean values and using the nearest neighbour algorithm. It can handle both numerical and nominal values. The function also delete the features with the number of missing values more then specified threshold.
  The results is in the form of \dQuote{list} with the processed dataset and the logical value, which indicates the success or failure of processing. The processed dataset can be used in the algorithms for feature selection \dQuote{select.process} and classification \dQuote{classifier.loop}.
}
\usage{
input_miss(matrix,method.subst="near.value",
attrs.nominal=numeric(),delThre=0.2)
}
\arguments{
  \item{matrix}{a dataset, a matrix of feature values for several cases, the last column is for the class labels. Class labels could be numerical or character values. The maximal number of classes is ten.}
  \item{method.subst}{a method of missing value processing.
  There are two realized methods: substitution with mean value
  ('mean.value') and nearest neighbour algorithm\cr  ('near.value').}
  \item{attrs.nominal}{a numerical vector, containing the column numbers of the nominal features, selected for the analysis.}
  \item{delThre}{the minimal threshold for the deletion of features with missing values. It is in the interval [0,1], where for delThre=0 all features having at least one missing value will be deleted.}
}
\details{
  This function's main job is to handle the missing values in the dataset. See the
  \dQuote{Value} section to this page for more details.

  Data can be provided in matrix form, where the rows correspond to cases with feature values and class label. The columns contain the values of individual features and the last column must contain class labels. The maximal number of class labels equals 10.
  The class label features and all the nominal features must be defined as factors.
}
\value{
  The data are provided with reasonable number of missing values that is preprocessed with one of the imputing methods.

  A returned list consists of the the following fields:
  \item{data}{a processed dataset}
  \item{flag.miss}{logical value; if TRUE the processing is successful, if FALSE the input dataset is returned without processing.}
  }

\references{
 McShane LM, Radmacher MD, Freidlin B, Yu R, Li MC, Simon R. Methods for assessing reproducibility of clustering patterns observed in analyses of microarray data. Bioinformatics. 2002 Nov;18(11):1462-9.
}

\seealso{
 \code{\link{select.process}}, \code{\link{classifier.loop}}
}

\examples{
# example for dataset with missing values
data(leukemia_miss)
xdata=leukemia_miss

# class label must be factor
xdata[,ncol(xdata)]<-as.factor(xdata[,ncol(xdata)])

# nominal features must be factors
attrs.nominal=101
xdata[,attrs.nominal]<-as.factor(xdata[,attrs.nominal])

delThre=0.2
out=input_miss(xdata,"mean.value",attrs.nominal,delThre)
if(out$flag.miss)
{
 xdata=out$data
}
}

\keyword{feature selection}
\keyword{classification}
\keyword{missing values}
