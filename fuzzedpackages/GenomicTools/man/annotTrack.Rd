\name{annotTrack}
\alias{annotTrack}
\docType{data}
\title{Example Annotation Track}
\description{
  An example for a typical annotation track.
}
\usage{data(annotTrack)}

\format{
  A data table with 1000 rows, each representing one annotation with 11 columns as provided from standard gtf format from Ensembl. 
}

\details{
  This is an example homo sapiens annotation track as it was downloaded from the Ensembl ftp download page. In total there are 1000
  annotations from the human genome in release 85.
}

\source{
  http://www.ensembl.org/info/data/ftp/index.html
}

\examples{
# The object was created from the downloaded Ensembl file as follows
\dontrun{
  ensGTF <- importGTF(file="Homo_sapiens.GRCh38.85.gtf.gz")
  annotTrack <- ensGTF[1:1000,]
  save(annotTrack,file="annotTrack.rda")
}

data(annotTrack)
annotTrack
}
\keyword{datasets}
 
