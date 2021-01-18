\name{ExampleBSseq}
\docType{data}
\alias{ExampleBSseq}
\title{BSseq example dataset}
\description{
  A simulated example dataset of BSseq for PQLseq.
}
\usage{data(ExampleBSseq)}
\format{
  Contains the following objects:
  \describe{
  	\item{mcount}{a data frame containing the read counts for 5 sites.}
    \item{predictor}{a vector of 100 observations on a continuous variable.}
    \item{relatednessmatrix}{a genetic relationship matrix for 100 individuals.}
    \item{totalcount}{a data frame containing the total read counts for 5 sites.}
  }
}
\examples{
data(ExampleBSseq)
attach(ExampleBSseq)
model_DNA=pqlseq(RawCountDataSet=mcount, Phenotypes=predictor, 
  RelatednessMatrix=relatednessmatrix, LibSize=totalcount,
  fit.model="BMM",numCore=1)
head(model_DNA)
detach(ExampleBSseq)
}
\keyword{datasets}

