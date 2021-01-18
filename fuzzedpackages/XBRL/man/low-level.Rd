\name{Low-level functions}
\alias{xbrlParse}
\alias{xbrlFree}
\alias{xbrlGetSchemaName}
\alias{xbrlGetLinkbaseNames}
\alias{xbrlGetImportNames}
\alias{xbrlProcessElements}
\alias{xbrlProcessLabels}
\alias{xbrlProcessArcs}
\alias{xbrlProcessRoles}
\alias{xbrlProcessContexts}
\alias{xbrlProcessFacts}
\alias{xbrlProcessUnits}
\alias{xbrlProcessFootnotes}
\title{Set of low-level functions to parse an XBRL file and extract
  information from it.}
\description{
  \code{xbrlParse} performs the parsing of the document. It returns a
  pointer to external memory used by the low-level
  functions. \code{xbrlFree} releases the external memory used.
  
  \code{xbrlGetSchemaName} is used on an instance document to extract the
  name of the associated schema. \code{xbrlGetLinkbaseNames} and
  \code{xbrlGetImportNames} are used on a taxonomy document to extract the
  names of linked documents during the process of discovery.
  
  \code{xbrlProcessXXX} functions are used to extract information from
  taxonomy and instance files in the form of data frames.
}
\usage{
## Use on any XBRL file
xbrlParse(file)
xbrlFree(doc)

## Use on an instance file
xbrlGetSchemaName(doc)

## Use on a taxonomy file
xbrlGetLinkbaseNames(doc)
xbrlGetImportNames(doc)

xbrlProcessElements(doc)
xbrlProcessLabels(doc)
xbrlProcessArcs(doc, arcType)

## Use on a schema file
xbrlProcessRoles(doc)

## Use on an instance file
xbrlProcessContexts(doc)
xbrlProcessFacts(doc)
xbrlProcessUnits(doc)
xbrlProcessFootnotes(doc)
}
\arguments{
  \item{file}{the name of the XBRL document. It can be a URL or local file name}
  \item{doc}{pointer to external memory, returned by \code{xbrlParse}}
  \item{arcType}{Either "presentation", "calculation", or "definition"}
}
\value{
  \code{xbrlParse} returns an external pointer.

  \code{xbrlFree} returns NULL.

  \code{xbrlGetSchemaName}, \code{xbrlGetLinkbaseNames} and
  \code{xbrlGetImportNames} return a character vector

  \code{xbrlProcessXXX} returns a data.frame
}
\author{
Roberto Bertolusso and Marek Kimmel
}
\seealso{
  \code{\link{xbrlDoAll}}
}
\examples{
\dontrun{
## Setting stringsAsFactors = FALSE is highly recommended
## to avoid data frames to create factors from character vectors.
options(stringsAsFactors = FALSE)

## XBRL instance file to be analyzed.
## NOTE: in this case, inst needs to be a local file, or accessible
##       as http (not https).
inst <- "ko-20130927.xml"
## Parse the instance (doc is an pointer to
## external memory that needs to be freed after use):
doc <- xbrlParse(inst)
## Get a data frame with facts:
fct <- xbrlProcessFacts(doc)
## Get a data frame with contexts:
cts <- xbrlProcessContexts(doc)
## Get a data frame with units:
unt <- xbrlProcessUnits(doc)
sch <- xbrlGetSchemaName(doc)
## Free the external memory used:
xbrlFree(doc)

dname <- dirname(inst)
## Parse the schema file:
docS <- xbrlParse(paste0(dname,"/",sch))
## Get roles:
rls <- xbrlProcessRoles(docS)
head(rls)
## Free the external memory used:
xbrlFree(docS)
}
}
\keyword{XBRL}
\keyword{low-level}

