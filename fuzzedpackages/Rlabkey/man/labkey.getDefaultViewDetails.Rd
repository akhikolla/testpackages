\name{labkey.getDefaultViewDetails}
\alias{labkey.getDefaultViewDetails}
\title{Retrieve the fields of a LabKey query view}
\description{
Fetch a list of output fields and their attributes that are avaialble from the default view of a given query}
\usage{
labkey.getDefaultViewDetails(baseUrl, folderPath,
    schemaName, queryName)
}
\arguments{
  \item{baseUrl}{a string specifying the \code{baseUrl}for the labkey server}
  \item{folderPath}{a string specifying the \code{folderPath} }
  \item{schemaName}{a string specifying the  \code{schemaName} for the query}
  \item{queryName}{a string specifying the \code{queryName}}
}
\details{
Queries have a default \dQuote{views} associeated with them.  A query view can describe a subset or superset of the fields defined by the query.  
A query view is defined by using the \dQuote{Customize View} button option on a LabKey data grid page.  \code{getDefaultViewDetails} has the same arguments 
and returns the same shape of result data frame as \code{getQueryDetails}.The default view is the what you will get back on calling 
\code{labkey.selectRows} or \code{getRows}.
}
\value{
The output field attributes of the default view are returned as a data frame.  See \code{\link{labkey.getQueryDetails}} for a description.
}
\author{Peter Hussey, peter@labkey.com}
\seealso{
{Retrieve data:} \code{\link{labkey.selectRows}}, \code{\link{makeFilter}}, \code{\link{labkey.executeSql}} \cr
{Modify data:  } \code{\link{labkey.updateRows}}, \code{\link{labkey.insertRows}}, \code{\link{labkey.importRows}}, \code{\link{labkey.deleteRows}}\cr
{List available data: } \code{\link{labkey.getSchemas}}, \code{\link{labkey.getQueries}}, \code{\link{labkey.getQueryViews}}, 
\code{\link{labkey.getQueryDetails}}, \code{\link{labkey.getLookupDetails}}
}
\examples{
\donttest{

## Details  of fields of a default query view
# library(Rlabkey)

queryDF <- labkey.getDefaultViewDetails(
	baseUrl="http://localhost:8080/labkey",
	folderPath="/apisamples",
	schemaName="lists",
	queryName="AllTypes")
	
queryDF

}
}
\keyword{IO}
