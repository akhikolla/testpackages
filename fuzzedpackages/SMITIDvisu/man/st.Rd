\name{st}
\alias{st}
\docType{data}
\title{
  A SMITIDstruct package variable.
}
\description{
  A SMITIDstruct package variable from simul-chain as a list.
  The list is a set of HostSet, ViralPopset and an index
}
\usage{data("st")}
\format{
  The format is:
List of 3
 $ :List of 21
  ..
  ..- attr(*, "class")= chr "HostSet"
 $ :List of 20
  ..
  ..- attr(*, "class")= chr "ViralPopSet"
 $ :'data.frame':	79 obs. of  3 variables:
  ..$ TIME     : chr [1:79] "0" "0" "1.26" "1.35" ...
  ..$ ID_HOST  : chr [1:79] "1" "2" "2" "2" ...
  ..$ EVENTCODE: chr [1:79] "000011" "000110" "001000" "001000" ...
}
\examples{
data(st)
## maybe str(st) ; plot(st) ...
}
\keyword{datasets}
