\name{Kryukov}
\docType{data}
\alias{Kryukov}

\title{Kryukov data set}

\description{
  The data from \emph{Kryukov et al, 2009}, contain simulated site frequency spectrum data using European demographic models with purifying selection.
}

\usage{data(Kryukov) }

\format{
  \code{Kryukov} is a dataframe with four columns:
  \describe{
    \item{gene}{The unit of each variant}
    \item{maf}{The maf of each variant in the European population}
    \item{selection.coefficient}{The selction coefficient of each variant in the European population}
    \item{position}{The position of each variant}
  }
}

\details{
  200 units are present corresponding to 200 genes. For each unit, the data set contains the maf in the European population, the selection coefficient and the position of each variant.
}

\source{ The data were obtained from the SeqPower software (see also \url{http://www.bioinformatics.org/spower/input#data_download}). }

\references{Kryukov et al, 2009, \emph{Power of deep, all-exon resequencing for discovery of human trait genes}, Proceedings of the National Academy of Sciences, DOI:10.1073/pnas.0812824106}

