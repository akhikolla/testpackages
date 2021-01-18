\name{labkey.importRows}
\alias{labkey.importRows}
\title{Import rows of data into a LabKey Server}
\description{
Bulk import rows of data into the database.
}
\usage{
labkey.importRows(baseUrl, folderPath,
    schemaName, queryName, toImport, na)
}
\arguments{
  \item{baseUrl}{a string specifying the \code{baseUrl}for the labkey server}
  \item{folderPath}{a string specifying the \code{folderPath} }
  \item{schemaName}{a string specifying the  \code{schemaName} for the query}
  \item{queryName}{a string specifying the  \code{queryName} }
  \item{toImport}{a data frame containing rows of data to be imported}
  \item{na}{(optional) the value to convert NA's to, defaults to NULL}
}
\details{
Multiple rows of data can be imported in bulk.  The \code{toImport} data frame must contain
values for each column in the dataset and must be created with the \code{stringsAsFactors} option
set to FALSE. The names of the data in the data frame must be the column names from the
LabKey Server.To import a value of NULL, use an empty string ("") in the data frame (regardless of the database column type).
Also, when importing data into a study dataset, the sequence number must be specified. \cr \cr
Note: requires LabKey server version 13.3 or greater.
}

\value{
A list is returned with named categories of \bold{command}, \bold{rowsAffected}, \bold{queryName}, \bold{containerPath} and \bold{schemaName}.
The \bold{schemaName}, \bold{queryName} and \bold{containerPath} properties contain the same schema, query
and folder path used in the request.  The
\bold{rowsAffected} property indicates the number of rows affected by the API action. This will typically be the same
number as passed in the request.
}
\author{Cory Nathe}
\seealso{
\code{\link{labkey.selectRows}}, \code{\link{labkey.executeSql}}, \code{\link{makeFilter}},
\code{\link{labkey.insertRows}}, \code{\link{labkey.updateRows}}, \cr
\code{\link{labkey.deleteRows}}
}
\examples{
\donttest{

## Note that users must have the necessary permissions in the database
## to be able to modify data through the use of these functions
# library(Rlabkey)

newrows <- data.frame(
	DisplayFld="Imported from R"
	, RequiredText="abc"
	, RequiredInt=1
	, stringsAsFactors=FALSE)
newrows = newrows[rep(1:nrow(newrows),each=5),]

importedInfo <- labkey.importRows("http://localhost:8080/labkey",
    folderPath="/apisamples", schemaName="lists", queryName="AllTypes",
    toImport=newrows)

importedInfo$rowsAffected

}
}
\keyword{IO}