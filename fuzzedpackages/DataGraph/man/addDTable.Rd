\name{addDTable}
\alias{addDTable}
\title{Save a data frame to a .dtable file}
\usage{
addDTable(path,data)
}
\arguments{
\item{path}{A path name}
\item{data}{A data frame}
}
\value{
Nothing returned
}
\description{
Adds a table to an already open table file.  Need to open it first by using the openDTable function.
}
\examples{
\dontrun{
openDTable("/tmp/test")
addDTable("/tmp/test",Orange)
closeDTable("/tmp/test")
}
}