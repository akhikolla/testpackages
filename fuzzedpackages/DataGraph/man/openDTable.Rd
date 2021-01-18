\name{openDTable}
\alias{openDTable}
\title{Open a .dtable file so you can write into it using addDTable}
\usage{
openDTable(path)
}
\arguments{
\item{path}{A path name}
}
\value{
Nothing returned
}
\description{
Saves the data frame to a file
}
\examples{
\dontrun{
openDTBin("/tmp/test")
addDTBin("/tmp/test","A",Orange)
closeDTBin("/tmp/test")
}
}