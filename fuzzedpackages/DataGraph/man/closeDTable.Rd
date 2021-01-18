\name{closeDTable}
\alias{closeDTable}
\title{Save a data frame to a .dtable file}
\usage{
closeDTable(path)
}
\arguments{
\item{path}{A path name}
}
\value{
Nothing returned
}
\description{
Closes a .dtable file
}
\examples{
\dontrun{
openDTBin("/tmp/test")
addDTBin("/tmp/test","A",Orange)
closeDTBin("/tmp/test")
}
}