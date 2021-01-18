\name{getLookups}
\alias{getLookups}
\title{ Get related data fields that are available to include in a query on a given query object}
\description{ Retrieve a related query object referenced by a lookup column in the current query }

\usage{
getLookups(session, lookupField) 
}
\arguments{
  \item{session}{the session key returned from getSession}
  \item{lookupField}{an object representing a lookup field on LabKey Server, a named member of a query object. }
}

\details{  Lookup fields in LabKey Server are the equivalent of declared foreign keys
}

\value{ A query object representing the related data set.  The fields of a lookup query object 
are usually added to the colSelect parameter in getRows, If a lookup query object is used as the query parameter in
getRows, the call will return all of the base query columns and all of the lookup query columns.
A lookup query object is very similar to base
query objects that are named elemenets of a schema object, A lookup query
object, however, does not have a parent schema object, it is only returned by getLookups.  Also, the field names in
a lookup query object are compound names relative to the base query object used in getLookups.
}
\references{https://www.labkey.org/Documentation/wiki-page.view?name=propertyFields}
\author{Peter Hussey}
\seealso{\code{\link{getSession}}, \code{\link{getRows}} \code{\link{getSchema}}}
\examples{
\donttest{

##  get fields from lookup tables and add to query
# library(Rlabkey)

s<- getSession(baseUrl="http://localhost:8080/labkey",
    folderPath="/apisamples")
 
scobj <- getSchema(s, "lists")

# can add fields from related queries
lucols <- getLookups(s, scobj$AllTypes$Category)

# keep going to other tables
lucols2 <- getLookups(s, lucols[["Category/Group"]])

cols <- c(names(scobj$AllTypes)[2:6], names(lucols)[2:4])

getRows(s, scobj$AllTypes, colSelect=paste(cols, sep=","))

}
}
\keyword{file}
