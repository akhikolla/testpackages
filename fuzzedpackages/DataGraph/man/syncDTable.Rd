\name{syncDTable}
\alias{syncDTable}
\title{Adds a .sync file with the same name as the .dtable file.}
\usage{
syncDTable(path)
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
openDTable("/tmp/test")
addDTable("/tmp/test",Orange)
syncDTable("/tmp/test")
addDTable("/tmp/test",Orange)
closeDTable("/tmp/test")
}
}