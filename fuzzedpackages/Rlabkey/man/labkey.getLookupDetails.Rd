\name{labkey.getLookupDetails}
\alias{labkey.getLookupDetails}
\title{Retrieve detailed information on a LabKey query}
\description{
Fetch a list of output columns and their attributes from the query referenced by a lookup field}
\usage{
labkey.getLookupDetails(baseUrl, folderPath,
    schemaName, queryName, lookupKey)
}
\arguments{
  \item{baseUrl}{a string specifying the address of the LabKey Server, including the context root}
  \item{folderPath}{a string specifying the hierarchy of folders to the current folder (container) for the operation, starting with the project folder }
  \item{schemaName}{a string specifying the schema name in which the query object is defined}
  \item{queryName}{a string specifying the name the query}
  \item{lookupKey}{a string specifying the qualified name of a lookup field (foreign key) relative to the query specified by queryName }
}
\details{
When \code{getQueryDetails} returns non-NA values for the lookupQueryName, the \code{getLookupDetails} function can be called
to enumerate the fields from the query referenced by the lookup.  These lookup fields can be added to the \code{colSelect} list of \code{selectRows}.
}
\value{
The available schemas are returned as a data frame, with the same columns
as detailed in \code{\link{labkey.getQueryDetails}}
}
\author{Peter Hussey, peter@labkey.com}
\seealso{
{Retrieve data:  } \code{\link{labkey.selectRows}},\code{\link{makeFilter}}, \code{\link{labkey.executeSql}} \cr
{Modify data:  } \code{\link{labkey.updateRows}}, \code{\link{labkey.insertRows}}, \code{\link{labkey.importRows}}, \code{\link{labkey.deleteRows}}\cr
{List available data: } \code{\link{labkey.getSchemas}}, \code{\link{labkey.getQueries}}, \code{\link{labkey.getQueryViews}}, 
\code{\link{labkey.getQueryDetails}}, \code{\link{labkey.getDefaultViewDetails}}
}
\examples{
\donttest{

## Details  of fields of a query referenced by a lookup field
# library(Rlabkey)

lu1  <- labkey.getLookupDetails(
	baseUrl="http://localhost:8080/labkey",
	folderPath="/apisamples",
	schemaName="lists",
	queryName="AllTypes",
	lookupKey="Category"
)
lu1

## When a lookup field points to a query object that itself has a lookup
## field, use a compound fieldkey consisting of the lookup fields from
## the base query object to the target lookupDetails, separated by
## forward slashes
lu2<- labkey.getLookupDetails(
	baseUrl="http://localhost:8080/labkey",
	folderPath="/apisamples",
	schemaName="lists",
	queryName="AllTypes",
	lookupKey="Category/Group"
)
lu2

## Now select a result set containing a field from the base query, a
## field from the 1st level of lookup, and one from the 2nd
rows<- labkey.selectRows(
	baseUrl="http://localhost:8080/labkey",
	folderPath="/apisamples",
	schemaName="lists",
	queryName="AllTypes",
	colSelect=c("DisplayFld","Category/Category","Category/Group/GroupName"), 
	colFilter = makeFilter(c("Category/Group/GroupName",
	   "NOT_EQUALS","TypeRange")), maxRows=20
 )
 rows

}
}
\keyword{IO}
