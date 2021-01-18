\name{geneEXP}
\alias{geneEXP}
\docType{data}
\title{Simulated Expression Data}
\description{
  Simulated expression data.
}
\usage{data(geneEXP)}

\format{
  A matrix with 50 rows representing individuals and 1000 genes arranged in the columns.
}

\details{
  Just simulated, numerical expression data for testing purposes. The genes are the same names as they can be found in the \code{annotTrack} object,
  the data consists out of 50 samples. The first 30 rows represent cases, the last 20 the controls. The expression data follows simply a N(0,1)
  distribution, however, the first 100 genes were shifted by +1. 
}

\source{
  Own simulation, code can be given upon request.
}

\examples{
data(Xgene)
}
\keyword{datasets}
 
