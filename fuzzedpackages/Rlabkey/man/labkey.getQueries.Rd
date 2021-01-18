\name{labkey.getQueries}
\alias{labkey.getQueries}
\title{Retrieve a list of available queries for a specified LabKey schema}
\description{
Fetch a list of queries available to the current user within in a specified folder context and specified schema
}
\usage{
labkey.getQueries(baseUrl, folderPath, schemaName)
}
\arguments{
  \item{baseUrl}{a string specifying the address of the LabKey Server, including the context root}
  \item{folderPath}{a string specifying the hierarchy of folders to the current folder (container) for the operation, starting with the project folder }
  \item{schemaName}{a string specifying the schema name in which the query object is defined}
}  
\details{
\dQuote{Query} is the LabKey term for a data container that acts like a relational table within LabKey Server. Queries 
include lists, assay data results, user-defined queries, built-in SQL tables in individual modules, and 
tables or table-like objects in external schemas,  For a specific queriable object,
the data that is visible depends on the current user's permissions in a given folder.
Function arguments identify the location of the server and the folder path.
}
\value{
The available queries are returned as a three-column data frame containing {one row for each field} for each query 
in the specified schema.  The three columns are

\code{queryName}{the name of the query object, repeated once for every field defined as output of the query.}  \cr
\code{fieldName} {the name of a query output field} \cr 
\code{caption}{the caption of the named field as shown in the column header of a data grid, also known as a label}\cr 

}
\references{http://www.omegahat.net/RCurl/,\cr
https://www.labkey.org/project/home/begin.view}
\author{Peter Hussey, peter@labkey.com}
\seealso{
{Retrieve data:  } \code{\link{labkey.selectRows}}, \code{\link{makeFilter}}, \code{\link{labkey.executeSql}} \cr
{Modify data:  } \code{\link{labkey.updateRows}}, \code{\link{labkey.insertRows}}, \code{\link{labkey.importRows}}, \code{\link{labkey.deleteRows}}\cr
{List available data: } \code{\link{labkey.getSchemas}}, \code{\link{labkey.getQueryViews}}, 
\code{\link{labkey.getQueryDetails}}, \code{\link{labkey.getDefaultViewDetails}}, \code{\link{labkey.getLookupDetails}}

}
\examples{
\donttest{

## List of queries in a schema
# library(Rlabkey)

queriesDF <- labkey.getQueries(
	baseUrl="https://www.labkey.org",
	folderPath="/home",
	schemaName="lists"
)
queriesDF

}
}
\keyword{IO}
