\name{sample_data}
\docType{data}
\alias{sample_data}
\title{Preprocessed 293T–Jurkat scRNA-seq dataset}
\description{
  This data set has been generated from 293T and Jurkat cell data containing a total of ~3200 cells,
  with an almost equal number of representative transcriptomes of each type.
  The cells were mixed in vitro at equal proportions.
  Authors of the study resolved the cell types bioinformatically exploiting their SNV profiles.
  With this data, we mimicked the rare cell phenomenon by bioinformatically diluting
  Jurkat cell proportion to 2.5\% in the data.
}
\usage{data(sample_data)}
\format{A matrix containing 1580 observations (in rows), with 1000 features per observation (in columns).}
\source{10xgenomics, https://support.10xgenomics.com/single-cell-gene-expression/datasets}
\references{
  Zheng, G. X. et al.
  Massively parallel digital transcriptional profiling of single cells.
  \emph{Nat. Commun.} 8, 14049 (2017).
}
\keyword{datasets}
