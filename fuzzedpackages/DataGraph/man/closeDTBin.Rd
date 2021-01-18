\name{closeDTBin}
\alias{closeDTBin}
\title{Close a .dtbin file}
\usage{
closeDTBin(path)
}
\arguments{
\item{path}{A path name}
}
\value{
Nothing returned
}
\description{
Closes a .dtbin file
}
\examples{
\dontrun{
openDTBin("/tmp/test")
addDTBin("/tmp/test","A",Orange)
closeDTBin("/tmp/test")
}
}