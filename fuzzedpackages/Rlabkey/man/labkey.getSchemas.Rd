\name{labkey.getSchemas}
\alias{labkey.getSchemas}
\title{Retrieve a list of available schemas from a labkey database}
\description{
Fetch a list of schemas available to the current user in a specified folder context
}
\usage{
labkey.getSchemas(baseUrl, folderPath)
}
\arguments{
  \item{baseUrl}{a string specifying the address of the LabKey Server, including the context root}
  \item{folderPath}{a string specifying the hierarchy of folders to the current folder (container) for the operation, starting with the project folder }
}
\details{
Schemas act as the name space for query objects in LabKey Server. Schemas are generatlly
associated with a LabKey Server "module" that provides some specific functionality. Within a queriable object,
the specific data that is visible depends on the current user's permissions in a given folder.
Function arguments are the components of the url that identify the location of the server and the folder path.

}
\value{
The available schemas are returned as a single-column data frame.
}
\references{http://www.omegahat.net/RCurl/,\cr
https://www.labkey.org/project/home/begin.view}
\author{Peter Hussey, peter@labkey.com}
\seealso{
{Retrieve data:  } \code{\link{labkey.selectRows}}, \code{\link{makeFilter}}, \code{\link{labkey.executeSql}} \cr
{Modify data:  } \code{\link{labkey.updateRows}}, \code{\link{labkey.insertRows}}, \code{\link{labkey.importRows}}, \code{\link{labkey.deleteRows}}\cr
{List available data: } \code{\link{labkey.getQueries}}, \code{\link{labkey.getQueryViews}}, 
\code{\link{labkey.getQueryDetails}}, \code{\link{labkey.getDefaultViewDetails}},
\code{\link{labkey.getLookupDetails}},

}
\examples{
\donttest{

## List of schemas
# library(Rlabkey)

schemasDF <- labkey.getSchemas(
	baseUrl="http://localhost:8080/labkey",
	folderPath="/apisamples"
)

}
}
\keyword{IO}
