\name{syncDTBin}
\alias{syncDTBin}
\title{Adds a .sync file with the same name as the .dtbin file.}
\usage{
syncDTBin(path)
}
\arguments{
\item{path}{A path name}
}
\value{
Nothing returned
}
\description{
Used so that DataGraph can read up to an internal boundary inside the data file and avoids reading a partial table.  DataGraph monitors this file and reloads the file when it sees a new sync file.
}
\examples{
\dontrun{
openDTBin("/tmp/test")
addDTBin("/tmp/test","Orange",Orange)
syncDTBin("/tmp/test")
addDTBin("/tmp/test","Orange",Orange)
closeDTBin("/tmp/test")
}
}