\encoding{UTF-8}
\name{CalcGene}
\alias{CalcGene}

\title{
 Calculate HUM value
}
\description{
  This is the auxiliary function of the Biocomb package. It computes a HUM value for individual feature and returns a \dQuote{List} object, consisting of HUM value and the best permutation of class labels in \dQuote{seq} vector. This \dQuote{seq} vector can be passed to the function \code{\link{CalculateHUM_ROC}}.
}
\usage{
CalcGene(s_data, seqAll, prodValue,thresholds)
}

\arguments{
  \item{s_data}{a list, which contains the vectors of sorted feature values for individual class labels.
  }
  \item{seqAll}{a numeric matrix of all the permutations of the class labels, where each row corresponds to individual permutation vector.
  }
  \item{prodValue}{a numeric value, which is the product of the sizes of feature vectors, corresponding to analized class labels.
  }
  \item{thresholds}{a numeric vector, containing the values of thresholds for calculating ROC curve coordinates.}
 }

\details{
  This function's main job is to compute the maximal HUM value between the all possible permutations of class labels for individual feature, selected for analysis. See the \dQuote{Value} section to this page for more details.

}

\value{
  The data must be provided without missing values in order to process. A returned list consists of the following fields:
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
 \code{\link{CalculateHUM_Ex}}, \code{\link{CalculateHUM_ROC}}
}

\examples{

data(leukemia72)
# Basic example
# class label must be factor
leukemia72[,ncol(leukemia72)]<-as.factor(leukemia72[,ncol(leukemia72)])

xdata=leukemia72
indexF=3
indexClass=ncol(xdata)
label=levels(xdata[,indexClass])

indexLabel=label[1:2]

indexL=NULL
for(i in 1:length(indexLabel))
{
  indexL=c(indexL,which(label==indexLabel[i]))
}

indexEach=NULL
indexUnion=NULL
for(i in 1:length(label))
{
  vrem=which(xdata[,indexClass]==label[i])
  indexEach=c(indexEach,list(vrem))
  if(length(intersect(label[i],indexLabel))==1)
  indexUnion=union(indexUnion,vrem)
}

s_data=NULL
dataV=xdata[,indexF]
prodValue=1

for (j in 1:length(indexLabel))
{
  vrem=sort(dataV[indexEach[[indexL[j]]]])
  s_data=c(s_data,list(vrem))
  prodValue = prodValue*length(vrem)
}

len=length(indexLabel)
seq=permutations(len,len,1:len)

#claculate the threshold values
thresholds <- sort(unique(dataV[indexUnion]))
thresholds=(c(-Inf, thresholds) + c(thresholds, +Inf))/2

out=CalcGene(s_data,seq,prodValue,thresholds)
}

\keyword{feature selection}
\keyword{classification}
\keyword{HUM values}
\keyword{AUC values}
\keyword{missing values}
