\name{labkey.getQueryViews}
\alias{labkey.getQueryViews}
\title{Retrieve a list of available named views defined on a query in a schema}
\description{
Fetch a list of named query views available to the current user in a specified folder context, schema and query
}
\usage{
labkey.getQueryViews(baseUrl, folderPath, schemaName, queryName)
}
\arguments{
  \item{baseUrl}{a string specifying the address of the LabKey Server, including the context root}
  \item{folderPath}{a string specifying the hierarchy of folders to the current folder (container) for the operation, starting with the project folder }
  \item{schemaName}{a string specifying the schema name in which the query object is defined}
  \item{queryName}{a string specifying the name the query}
}
\details{
Queries have a default \dQuote{view} associeated with them, and can also have any number of named views.
A named query view is created by using the \dQuote{Customize View} button option on a LabKey data grid page.  
Use \code{getDefaultViewDetails} to get inforation about the default (unnamed) view.
}
\value{
The available views for a query are returned as a three-column data frame, with one row per view output field.

\code{viewName} {The name of the view, or NA for the default view.}\cr
\code{fieldName} {The name of a field within the view, as defined in the query object to which the field belongs}\cr 
\code{key} {The name of the field relative to the base query,  Use this value in the colSelect parameter of \code{labkey.selectRows()} .}\cr 

}
\references{https://www.labkey.org/Documentation/wiki-page.view?name=savingViews}
\author{Peter Hussey, peter@labkey.com}
\seealso{
{Retrieve data:  } \code{\link{labkey.selectRows}},\code{\link{makeFilter}}, \code{\link{labkey.executeSql}} \cr
{Modify data:  } \code{\link{labkey.updateRows}}, \code{\link{labkey.insertRows}}, \code{\link{labkey.importRows}}, \code{\link{labkey.deleteRows}} \cr
{List available data: } \code{\link{labkey.getSchemas}}, \code{\link{labkey.getQueries}},
\code{\link{labkey.getQueryDetails}},\code{\link{labkey.getDefaultViewDetails}},\code{\link{labkey.getLookupDetails}}


}
\examples{
\donttest{

## List of views defined for a query in a schema
# library(Rlabkey)

viewsDF <- labkey.getQueryViews(
	baseUrl="http://localhost:8080/labkey",
	folderPath="/apisamples",
	schemaName="lists",
	queryName="AllTypes"
)

}
}
\keyword{IO}
