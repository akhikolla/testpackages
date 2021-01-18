\name{readkvs}
\alias{readkvs}
\title{Read key/value pair observations.}
\arguments{
  \item{path}{The path to the input file.}

  \item{delim_tuple}{The tuple deliminator.}

  \item{delim_field}{The field deliminator.}
}
\value{
A list of observations. Each observation is a named list of
key/value pairs.
}
\description{
Provides the workhorse function to read key/value pair
observations from a flat file. Note that, same keys and
corresponding values in the same observation are returned
as vectors.
}
\examples{
\dontrun{
## Assume that we have a file at path "file.dat" with the
## following contents:
##
## key1=value11;key2=value12
## key1=value21;key2=value22
## key1=value31;key2=value320;key2=value321
##
## Let's read this in and display as a matrix:
observations <- readkvs("file.dat", ";", "=")
do.call(rbind, observations)
}
}

