\name{getSchema}
\alias{getSchema}
\title{ Returns an object representing a LabKey schema }
\description{ Creates and returns an object representing a LabKey schema, containing child objects representing LabKey queries}

\usage{
getSchema(session, schemaIndex) 
}
\arguments{
  \item{session}{the session key returned from getSession}
  \item{schemaIndex}{the name of the schema that contains the table on which you want to base a query, or the number of that schema as displayed by print(session)}
}
\details{Creates and returns an object representing a LabKey schema, containing child objects representing LabKey queries.  This compound object is created
by calling \code{labkey.getQueries} on the requested schema and \code{labkey.getQueryDetails} on each returned query.  The information returned in the
schema objects is essentially the same as the schema and query objects shown in the Schema Browser on LabKey Server.
}

\value{an object representing the schema.  The named elements of a schema are the queries within that schema.  
}
\references{https://www.labkey.org/Documentation/wiki-page.view?name=querySchemaBrowser}
\author{Peter Hussey}
\seealso{\code{\link{getSession}}}
\examples{
\donttest{

## the basics of using session, schema, and query objects
# library(Rlabkey)

s<- getSession(baseUrl="http://localhost:8080/labkey",
    folderPath="/apisamples")

sch<- getSchema(s, "lists")

# can walk down the populataed schema tree from schema node or query node
sch$AllTypes$Category
sch$AllTypes$Category$caption
sch$AllTypes$Category$type

# can add fields from related queries
lucols <- getLookups(s, sch$AllTypes$Category)

cols <- c(names(sch$AllTypes[2:6]), names(lucols)[2:4])

getRows(s, sch$AllTypes, colSelect=cols)

}
}
\keyword{file}
