\name{writeDTable}
\alias{writeDTable}
\title{Save a data frame to a .dtable file}
\usage{
writeDTable(path,data)
}
\arguments{
\item{path}{A path name}
\item{data}{A data frame}
}
\value{
Nothing returned
}
\description{
Saves the data frame to a file
}
\examples{
\dontrun{
writeDTable("/tmp/test",Orange)
}
}