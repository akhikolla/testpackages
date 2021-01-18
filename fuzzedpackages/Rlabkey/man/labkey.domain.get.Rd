\name{labkey.domain.get}
\alias{labkey.domain.get}
\title{Returns the metadata for an existing LabKey domain}
\description{
Get the data structure for a domain.
}
\usage{
labkey.domain.get(baseUrl=NULL, folderPath, schemaName, queryName)
}
\arguments{
  \item{baseUrl}{a string specifying the \code{baseUrl}for the labkey server}
  \item{folderPath}{a string specifying the \code{folderPath} }
  \item{schemaName}{a string specifying the name of the schema of the domain }
  \item{queryName}{a string specifying the query name }
}
\details{
Returns the domain design of an existing domain. The returned domain design can be used for reporting purposes or it can
be modified and used to create a new domain or update the domain source.
}
\value{
A list containing elements describing the domain. The structure is the same as a domain design created by \code{labkey.createDomainDesign}
}

\author{Karl Lum}
\seealso{
\code{\link{labkey.domain.create}},
\code{\link{labkey.domain.inferFields}},
\code{\link{labkey.domain.createDesign}},
\code{\link{labkey.domain.createIndices}},
\code{\link{labkey.domain.save}},
\code{\link{labkey.domain.drop}}
}
\examples{
\donttest{

## retrieve an existing domain
library(Rlabkey)

labkey.domain.get(baseUrl="http://labkey/", folderPath="home",
    schemaName="lists", queryName="test list")
}
}
\keyword{IO}
