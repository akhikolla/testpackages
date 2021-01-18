\name{XBRL-package}
\alias{XBRL-package}
\docType{package}
\title{Extraction of business financial information from XBRL documents.}
\description{
  Functions to extract business financial information from
  an Extensible Business Reporting Language (XBRL) instance file and the
  associated collection of files that defines its Discoverable Taxonomy
  Set (DTS), usually disseminated across different remote locations.

  XBRL documents employ many technologies defined by XML standards, such
  as XML Schema, Namespaces, XLink, and XPath, which make the extraction
  of information involved. The functions provided by this package
  address this complexity, returning data frames containing the enclosed
  information readily available to perform analyses.

  XBRL has been successfully tested in analyzing 10-K and 10-Q
  submissions from USA filing companies that are publicly available at
  the EDGAR system of the Securities and Exchange Commission (SEC).
  It has been reported to work with inline XBRL documents.

  As published taxonomy files (such as
  "http://xbrl.fasb.org/us-gaap/2013/elts/us-gaap-2013-01-31.xsd") 
  are immutable and are used by most filers, XBRL offers the option of
  downloading them only the first time they are referred, keeping a
  local cache copy that can be used from then on. This speeds-up the
  process, specially on slow connections.

  XBRL is not and does not aim to be (at least yet) a validating XBRL
  processor.

  XBRL is still a work in progress. As such, functions and produced data
  frames may change structure in future versions until stabilization.
}
\details{
\tabular{ll}{
Package: \tab XBRL\cr
Type: \tab Package\cr
Version: \tab 0.99.11\cr
Date: \tab 2014-01-01\cr
License: \tab GPL (>=2)\cr
}
The package offers 3 levels of access:

1) A function that "does it all" (see \code{\link{xbrlDoAll}}).

2) A "mutable state" function that exposes "methods" to be called
individually (see \code{\link{XBRL}}).

3) Individual specialized functions in C++, glued to a corresponding R
calling function (see \code{\link{xbrlParse}}).

See examples of use below.
}
\author{
Roberto Bertolusso and Marek Kimmel

Maintainer: Roberto Bertolusso <rbertolusso@rice.edu>
}
\references{
  http://www.xbrl.org
}
\seealso{
  \code{\link{xbrlDoAll}}, \code{\link{XBRL}}, \code{\link{xbrlParse}}
}
\keyword{ package XBRL XML }
\examples{
\dontrun{
## Setting stringsAsFactors = FALSE is highly recommended
## to avoid data frames to create factors from character vectors.
options(stringsAsFactors = FALSE)

## Load the library
library(XBRL)

## XBRL instance file to be analyzed, accessed
## directly from SEC website:
inst <- "https://www.sec.gov/Archives/edgar/data/21344/000002134413000050/ko-20130927.xml"

## Level 1: Function that does all work and returns
## a list of data frames with extracted information:
xbrl.vars <- xbrlDoAll(inst, verbose=TRUE)

## Level 2: Using the XBRL() "mutable state" function:
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

## Level 3: Using specialized functions that call C++ code directly:
## Parse the instance (doc is an pointer to external memory that needs to be freed):
## NOTE: in this case, inst needs to be a local file, or accessible
##       as http (not https).
inst <- "ko-20130927.xml"
doc <- xbrlParse(inst)
## Get a data frame with facts:
fct <- xbrlProcessFacts(doc)
## Get a data frame with contexts:
cts <- xbrlProcessContexts(doc)
## Get a data frame with units:
unt <- xbrlProcessUnits(doc)
## Free the external memory used:
xbrlFree(doc)
}
}

