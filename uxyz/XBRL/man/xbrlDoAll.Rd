\name{xbrlDoAll}
\alias{xbrlDoAll}
\title{Function to do all required work on an XBRL instance and its
  associated DTS.}
\description{
  \code{xbrlDoAll} performs all the steps necessary to extract a list
  of dataframes, starting from an XBRL instance file. \code{xbrlDoAll}
  offers the highest level of abstraction from the required steps. A
  medium-level access can be achieved by using \code{\link{XBRL}}. A
  low-level access can be achieved by using the functions related to
  \code{\link{xbrlParse}}.
}
\usage{
xbrlDoAll(file.inst, cache.dir = "xbrl.Cache", prefix.out = NULL,
          verbose = FALSE, delete.cached.inst = TRUE)
}
\arguments{
  \item{file.inst}{the name of the instance file}
  \item{cache.dir}{the name of a local directory to be used as a cache
    of discovered taxonomies. if NULL, no cache is performed}
  \item{prefix.out}{if given, prefix used to create csv files with the
    content of the data frames produced}
  \item{verbose}{Set to TRUE to see output of the process of discovery
  and analysis if documents}
  \item{delete.cached.inst}{Set to TRUE to delete downloaded instance,
  schema and linkbase files}
}
\value{a list of data frames with extracted information}
\author{
Roberto Bertolusso and Marek Kimmel
}
\seealso{
  \code{\link{XBRL}}, \code{\link{xbrlParse}}
}
\examples{
\dontrun{
## Setting stringsAsFactors = FALSE is highly recommended
## to avoid data frames to create factors from character vectors.
options(stringsAsFactors = FALSE)

## XBRL instance file to be analyzed, accessed
## directly from SEC website:
inst <- "https://www.sec.gov/Archives/edgar/data/21344/000002134413000050/ko-20130927.xml"

xbrl.vars <- xbrlDoAll(inst, cache.dir="XBRLcache", prefix.out="out", verbose=TRUE)
}
}
\keyword{XBRL}
\keyword{high-level}

