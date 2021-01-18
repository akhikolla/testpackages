\name{openDTBin}
\alias{openDTBin}
\title{Open a new .dtbin file}
\usage{
openDTBin(path)
}
\arguments{
\item{path}{A path name}
}
\value{
Nothing returned
}
\description{
Opens a .dtbin file.  Prints out an error message if the file is already open.
}
\examples{
\dontrun{
openDTBin("/tmp/test")
addDTBin("/tmp/test","A",Orange)
infoDTBin("/tmp/test")
closeDTBin("/tmp/test")
}
}