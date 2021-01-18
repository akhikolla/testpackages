\name{infoDTBin}
\alias{infoDTBin}
\title{Save a data frame to a .dtable file}
\usage{
infoDTBin(path)
}
\arguments{
\item{path}{A path name}
}
\value{
Nothing returned, but data is printed to the console.  Intended as a debugging aid.
}
\description{
Displays the content of a .dtbin file
}
\examples{
\dontrun{
openDTBin("/tmp/test")
addDTBin("/tmp/test","A",Orange)
infoDTBin("/tmp/test")
closeDTBin("/tmp/test")
}
}