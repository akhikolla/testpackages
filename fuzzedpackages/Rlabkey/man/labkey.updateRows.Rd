\name{labkey.updateRows}
\alias{labkey.updateRows}
\title{Update existing rows of data in a labkey database}
\description{
Send data from an R session to update existing rows of data in the database.
}
\usage{
labkey.updateRows(baseUrl, folderPath,
    schemaName, queryName, toUpdate, provenanceParams=NULL)
}
\arguments{
  \item{baseUrl}{a string specifying the \code{baseUrl}for the labkey server}
  \item{folderPath}{a string specifying the \code{folderPath} }
  \item{schemaName}{a string specifying the  \code{schemaName}for the query}
  \item{queryName}{a string specifying the  \code{queryName}} 
  \item{toUpdate}{a data frame containing the row(s) of data to be updated}
  \item{provenanceParams}{the provenance parameter object which contains the options to include as part of a provenance recording.
    This is a premium feature and requires the Provenance LabKey module to function correctly, if it is not present this parameter will be ignored.
  }
}
\details{
A single row or multiple rows of data can be updated.  The \code{toUpdate} data frame should contain 
the rows of data to be updated and must be created with the \code{stringsAsFactors} option
set to FALSE. The names of the data in the data frame must be the column names from the
labkey database.   To update a row/column to a value of NULL, use an empty string ("") in the data frame (regardless of the database column type).

}
\value{
A list is returned with named categories of \bold{command}, \bold{rowsAffected}, \bold{rows}, \bold{queryName}, \bold{containerPath} and \bold{schemaName}.
The \bold{schemaName}, \bold{queryName} and \bold{containerPath} properties contain the same schema, query 
and folder path used in the request.  The
\bold{rowsAffected} property indicates the number of rows affected by the API action. This will typically be the same
number as passed in the request. The \bold{rows} property contains a list of row objects corresponding to the rows
updated.
}
\author{Valerie Obenchain}
\seealso{
\code{\link{labkey.selectRows}}, \code{\link{labkey.executeSql}}, \code{\link{makeFilter}}, 
\code{\link{labkey.insertRows}}, \code{\link{labkey.importRows}}, \cr
\code{\link{labkey.deleteRows}},
\code{\link{labkey.provenance.createProvenanceParams}},
\code{\link{labkey.provenance.startRecording}},
\code{\link{labkey.provenance.addRecordingStep}},
\code{\link{labkey.provenance.stopRecording}}
}
\examples{
\donttest{

## Insert, update and delete
## Note that users must have the necessary permissions in the database
## to be able to modify data through the use of these functions
# library(Rlabkey)

newrow <- data.frame(
	DisplayFld="Inserted from R"
	, TextFld="how its done"
	, IntFld= 98 
	, DoubleFld = 12.345
	, DateTimeFld = "03/01/2010"
	, BooleanFld= FALSE
	, LongTextFld = "Four score and seven years ago"
#	, AttachmentFld = NA    		#attachment fields not supported 
	, RequiredText = "Veni, vidi, vici"
	, RequiredInt = 0
	, Category = "LOOKUP2"
	, stringsAsFactors=FALSE)

insertedRow <- labkey.insertRows("http://localhost:8080/labkey",
    folderPath="/apisamples", schemaName="lists", queryName="AllTypes",
    toInsert=newrow)
newRowId <- insertedRow$rows[[1]]$RowId

selectedRow<-labkey.selectRows("http://localhost:8080/labkey",
    folderPath="/apisamples", schemaName="lists", queryName="AllTypes",
    colFilter=makeFilter(c("RowId", "EQUALS", newRowId)))
selectedRow

updaterow=data.frame(
	RowId=newRowId
	, DisplayFld="Updated from R"
	, TextFld="how to update"
	, IntFld= 777 
	, stringsAsFactors=FALSE)

updatedRow <- labkey.updateRows("http://localhost:8080/labkey",
    folderPath="/apisamples", schemaName="lists", queryName="AllTypes",
    toUpdate=updaterow)
selectedRow<-labkey.selectRows("http://localhost:8080/labkey",
    folderPath="/apisamples", schemaName="lists", queryName="AllTypes",
    colFilter=makeFilter(c("RowId", "EQUALS", newRowId)))
selectedRow

deleterow <- data.frame(RowId=newRowId, stringsAsFactors=FALSE)
result <- labkey.deleteRows(baseUrl="http://localhost:8080/labkey",
    folderPath="/apisamples", schemaName="lists", queryName="AllTypes",
    toDelete=deleterow)
str(result)

}
}
\keyword{IO}
