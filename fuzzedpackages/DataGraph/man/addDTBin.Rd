\name{addDTBin}
\alias{addDTBin}
\title{Save a data frame to a .dtable file}
\usage{
addDTBin(path,name,data,time)
}
\arguments{
\item{path}{A path name}
\item{name}{The name of the variable}
\item{data}{Either a data frame or a vector}
\item{time}{optional: What time value this variable is at}
}
\value{
Nothing returned
}
\description{
Adds the data into a .dtbin file.
}
\examples{
\dontrun{
openDTBin("/tmp/test")
addDTBin("/tmp/test","A",Orange)
closeDTBin("/tmp/test")
}
}