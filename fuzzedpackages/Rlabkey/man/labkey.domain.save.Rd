\name{labkey.domain.save}
\alias{labkey.domain.save}
\title{Updates an existing LabKey domain}
\description{
Modify an existing domain with the specified domain design.
}
\usage{
labkey.domain.save(baseUrl=NULL, folderPath, schemaName, queryName, domainDesign)
}
\arguments{
  \item{baseUrl}{a string specifying the \code{baseUrl}for the labkey server}
  \item{folderPath}{a string specifying the \code{folderPath} }
  \item{schemaName}{a string specifying the name of the schema of the domain }
  \item{queryName}{a string specifying the query name }
  \item{domainDesign}{a list data structure with the domain design to update to }
}
\value{
A list containing elements describing the domain after the update. The structure is the same as a domain design created by \code{labkey.createDomainDesign}
}

\author{Karl Lum}
\seealso{
\code{\link{labkey.domain.get}},
\code{\link{labkey.domain.inferFields}},
\code{\link{labkey.domain.createDesign}},
\code{\link{labkey.domain.createIndices}},
\code{\link{labkey.domain.create}},
\code{\link{labkey.domain.drop}},
\code{\link{labkey.domain.createConditionalFormat}},
\code{\link{labkey.domain.createConditionalFormatQueryFilter}},
\code{\link{labkey.domain.FILTER_TYPES}}
}
\examples{
\donttest{

library(Rlabkey)
## change the type of one of the columns
domain <- labkey.domain.get(baseUrl="http://labkey/", folderPath="home",
    schemaName="lists", queryName="test list")

domain$fields[3,]$rangeURI = "xsd:string"
domain$fields[3,]$name = "changed to string"

labkey.domain.save(baseUrl="http://labkey/", folderPath="home",
    schemaName="lists", queryName="test list", domainDesign=domain)
}
}
\keyword{IO}
