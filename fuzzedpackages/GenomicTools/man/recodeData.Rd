\name{recodeData}
\alias{recodeData}

\title{
  Recode a Genotype Data Matrix.
}

\description{
  This function recodes a genotype data matrix in a format as expected from \code{eqtl} or \code{mdr}.
}
\usage{
  recodeData(X)
}

\arguments{
  \item{X}{Matrix with genotype information.}
}

\details{
  This function recodes the values given in the data matrix \code{X} (typically \code{AA}, \code{AB} and \code{BB}) and
  substitutes it with 0,1 and 2. Missing values are encoded as 3.
}
\value{
 A matrix with same dimension as the input, but in 0,1,2 encoding.
}

\seealso{
  \code{\link{eQTL}} , \code{\link{mdr}}
}

\author{
  Daniel Fischer
}

\examples{
# genotData <- read.table("MDR_format_ready_BCR.txt",header=T)

# temp <- recodeData(genotData)

}

\keyword{
  multivariate
}
 
