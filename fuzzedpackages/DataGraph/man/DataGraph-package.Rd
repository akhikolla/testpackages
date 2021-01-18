\name{DataGraph-package}
\alias{DataGraph-package}
\alias{DataGraph}
\docType{package}
\title{
  Export Data from R so DataGraph can Read it
}
\description{
Functions to save either .dtable or .dtbin files that can be read by DataGraph,
a graphing and analysis application for MacOS.  Can save a data frame,
collection of data frames and sequences of data frames and individual vectors.
}
\details{
There is a total of eight functions in this package
One step to save a data frame into a .dtable file
writeDTable(path name,data frame)

To save multiple data frames into a single .dtable file
openDTable(path name)
addDTable(path name,data frame)
closeDTable(path name)

To save a .dtbin file, which can contain multiple lists and data frames by name and
each one can be a sequence saved by "time"
openDTBin(path name)
addDTBin(path name,entry name,data (frame or column),time (optional))
infoDTBin(path name)  : Prints out the current content of the file
closeDTBin(path name)
}
\author{
David Adalsteinsson

Maintainer: David Adalsteinsson <david@visualdatatools.com>
}
\keyword{ package }
\examples{
  \dontrun{
     writeDTable("/tmp/test",mtcars)

     openDTable("/tmp/test2" name)
     addDTable("/tmp/test2",data frame)
     ...  modify the table
     addDTable("/tmp/test2",data frame)
     closeDTable("/tmp/test2")

     openDTBin("/tmp/test3")
     addDTBin("/tmp/test3","Cars",mtcars)
     addDTBin("/tmp/test3","Orange",Orange)
     infoDTBin("/tmp/test3")
     closeDTBin("/tmp/test3")
  }
}
