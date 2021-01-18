\name{XBRL}
\alias{XBRL}
\title{"Mutable state" function that exposes "methods" that analyze an
  XBRL instance an its associated DTS.}
\description{
  \code{XBRL} is a "mutable state" function that offers a object-oriented
  programming-like interface, by exposing "methods" that perform accions
  on XBRL files. It keeps a list of data frames that are populated by
  the "methods" by calls to the low-level functions (see
  \code{\link{xbrlParse}}).
  
  The approach used is derived from an example included in *10.7 Scope*
  of *An Introduction to R*. The employed terminology may be inaccurate,
  hence the quotations. \code{XBRL} "methods" are called by
  \code{\link{xbrlDoAll}}. See examples of use below. 
}
\usage{
XBRL()
}
\arguments{no arguments}
\value{a list of "methods"}
\author{
Roberto Bertolusso and Marek Kimmel
}
\seealso{
  \code{\link{xbrlDoAll}}, \code{\link{xbrlParse}}
}
\examples{
\dontrun{
## Setting stringsAsFactors = FALSE is highly recommended
## to avoid data frames to create factors from character vectors.
options(stringsAsFactors = FALSE)

## XBRL instance file to be analyzed, accessed
## directly from SEC website:
inst <- "https://www.sec.gov/Archives/edgar/data/21344/000002134413000050/ko-20130927.xml"

xbrl <- XBRL()
xbrl$setCacheDir("XBRLcache")
xbrl$openInstance(inst)
## Perform a discovery of the taxonomy:
xbrl$processSchema(xbrl$getSchemaName())
## Process instance file:
xbrl$processContexts()
xbrl$processUnits()
xbrl$processFacts()
xbrl$processFootnotes()
xbrl$closeInstance()
xbrl.vars <- xbrl$getResults()
}
}
\keyword{XBRL}
\keyword{mid-level}

