\name{labkey.deleteRows}
\alias{labkey.deleteRows}
\title{Delete rows of data from a LabKey database}
\description{
Specify rows of data to be deleted from the LabKey Server
}
\usage{
labkey.deleteRows(baseUrl, folderPath,
    schemaName, queryName, toDelete, provenanceParams=NULL)
}
\arguments{
  \item{baseUrl}{a string specifying the \code{baseUrl}for LabKey server}
  \item{folderPath}{a string specifying the \code{folderPath}}
  \item{schemaName}{a string specifying the  \code{schemaName} for the query}
  \item{queryName}{a string specifying the  \code{queryName}}
  \item{toDelete}{a data frame containing a single column of data containing the data identifiers of the rows to be deleted} 
  \item{provenanceParams}{the provenance parameter object which contains the options to include as part of a provenance recording.
    This is a premium feature and requires the Provenance LabKey module to function correctly, if it is not present this parameter will be ignored.
  }
}
\details{
A single row or multiple rows of data can be deleted.  For the \code{toDelete} data frame, version 0.0.5 or later accepts either a single column 
of data containing the data identifiers  (e.g., key or lsid) or the entire row of data to be deleted.
The names of the data in the data frame must be the column names from the LabKey Server.
The data frame must be created with the \code{stringsAsFactors} set to FALSE.

NOTE: Each variable in a dataset has both a column label and a column name. The column label is visible at the top
of each column on the web page and is longer and more descriptive. The column name is shorter and is 
used \dQuote{behind the scenes} for database manipulation. It is the column name that must be used in
the Rlabkey functions when a column name is expected. To identify a particular column name in a dataset on
a web site, use the \dQuote{export to R script} option available as a drop down option under the \dQuote{views} 
tab for each dataset.
}
\value{
A list is returned with named categories of \bold{command}, \bold{rowsAffected}, \bold{rows}, \bold{queryName}, \bold{containerPath} and \bold{schemaName}.
The \bold{schemaName}, \bold{queryName} and \bold{containerPath} properties contain the same schema, query 
and folder path used in the request.  The
\bold{rowsAffected} property indicates the number of rows affected by the API action. This will typically be the same
number as passed in the request. The \bold{rows} property contains a list of rows corresponding to the rows
deleted.
}
\author{Valerie Obenchain}
\seealso{
\code{\link{labkey.selectRows}}, \code{\link{labkey.executeSql}}, \code{\link{makeFilter}}, 
\code{\link{labkey.insertRows}}, \code{\link{labkey.importRows}}, \cr
\code{\link{labkey.updateRows}},
\code{\link{labkey.provenance.createProvenanceParams}},
\code{\link{labkey.provenance.startRecording}},
\code{\link{labkey.provenance.addRecordingStep}},
\code{\link{labkey.provenance.stopRecording}}
}

\examples{
\donttest{

## Insert, update and delete
## Note that users must have the necessary permissions in the LabKey Server
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
#	, AttachmentFld = NA    #attachment fields not supported 
	, RequiredText = "Veni, vidi, vici"
	, RequiredInt = 0
	, Category = "LOOKUP2"
	, stringsAsFactors=FALSE)

insertedRow <- labkey.insertRows("http://localhost:8080/labkey",
    folderPath="/apisamples", schemaName="lists",
    queryName="AllTypes", toInsert=newrow)
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
    folderPath="/apisamples", schemaName="lists",
    queryName="AllTypes", toUpdate=updaterow)
selectedRow<-labkey.selectRows("http://localhost:8080/labkey",
    folderPath="/apisamples", schemaName="lists", queryName="AllTypes",
    colFilter=makeFilter(c("RowId", "EQUALS", newRowId)))
selectedRow

deleterow <- data.frame(RowId=newRowId, stringsAsFactors=FALSE)
result <- labkey.deleteRows(baseUrl="http://localhost:8080/labkey",
    folderPath="/apisamples", schemaName="lists",
    queryName="AllTypes",  toDelete=deleterow)
result

}
}
\keyword{IO}
