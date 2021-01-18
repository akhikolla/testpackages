\encoding{UTF-8}
\name{CalcGene}
\alias{CalcGene}

\title{
 Calculate HUM value
}
\description{
  This is the auxiliary function of the HUM package. It computes a HUM value for individual feature and returns a \dQuote{List} object, consisting of HUM value and the best permutation of class labels in \dQuote{seq} vector. This \dQuote{seq} vector can be passed to the function \code{\link{CalculateHUM_ROC}}. 
}
\usage{
CalcGene(s_data, seqAll, prodValue)
}

\arguments{
  \item{s_data}{a list, which contains the vectors of sorted feature values for individual class labels.
  }
  \item{seqAll}{a numeric matrix of all the permutations of the class labels, where each row corresponds to individual permutation vector.
  }
  \item{prodValue}{a numeric value, which is the product of the sizes of feature vectors, corresponding to analized class labels.
  }
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
data(sim)
# Basic example
indexF=3;
indexClass=2;
indexLabel=c("Normal","OrthArthr")
s_data=NULL;
prodValue=1;
for(i in 1:length(indexLabel))
{
  index=which(sim[,indexClass]==indexLabel[i])
  vrem=sort(sim[index,indexF])
  s_data=c(s_data,list(vrem))
  prodValue=prodValue*length(index)
}
len=length(indexLabel)
seqAll=permutations(len,len,1:len)
out=CalcGene(s_data, seqAll, prodValue)
}

\keyword{multivariate}
\keyword{nonparametric}
\keyword{plot}
\keyword{roc}
\keyword{auc}
