\name{kvh-package}
\alias{kvh-package}
\alias{kvh}
\docType{package}
\title{
  Write/read KVH (key-value hierarchy) file
}
\description{
   The format KVH is a lightweight format that can be read/written both by humans and machines.
   It can be useful in situations where XML or alike formats seem to be an overkill.
   We provide an ability to manipulate kvh files in R with a good efficiency due to Rcpp use.
   The content read in kvh file is hierarchically organized in nested lists. The key-values are
   always returned as character strings. It's up to user to convert them further in useful types
   (numeric vectors, matrices and so on).
}
\author{
  Serguei Sokol.

  Maintainer: Serguei Sokol <sokol@insa-toulouse.fr>
}
\references{
  http://serguei.sokol.free.fr/kvh-format/
}
\keyword{ kvh-format }
\examples{
  \dontrun{
     # prepare object to write to kvh file
     obj=list(x=structure(1:3, names=letters[1:3]), R=R.version)
     # write it
     obj2kvh(obj, "test", "test.kvh") # will create test.kvh file
     # read it back
     l=kvh_read("test.kvh")
     # check a field
     l$test$x # NB. it has a character values put in a list not a numeric vector as it was in obj.
  }
}
