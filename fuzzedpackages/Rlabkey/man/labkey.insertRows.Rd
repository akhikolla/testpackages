\name{labkey.insertRows}
\alias{labkey.insertRows}
\title{Insert new rows of data into a LabKey Server}
\description{
Insert new rows of data into the database.
}
\usage{
labkey.insertRows(baseUrl, folderPath,
    schemaName, queryName, toInsert, na, provenanceParams=NULL)
}
\arguments{
  \item{baseUrl}{a string specifying the \code{baseUrl}for the labkey server}
  \item{folderPath}{a string specifying the \code{folderPath} }
  \item{schemaName}{a string specifying the  \code{schemaName} for the query}
  \item{queryName}{a string specifying the  \code{queryName} }
  \item{toInsert}{a data frame containing rows of data to be inserted}
  \item{na}{(optional) the value to convert NA's to, defaults to NULL}
  \item{provenanceParams}{the provenance parameter object which contains the options to include as part of a provenance recording.
    This is a premium feature and requires the Provenance LabKey module to function correctly, if it is not present this parameter will be ignored.
  }
}
\details{
A single row or multiple rows of data can be inserted.  The \code{toInsert} data frame must contain
values for each column in the dataset and must be created with the \code{stringsAsFactors} option
set to FALSE. The names of the data in the data frame must be the column names from the
LabKey Server.To insert a value of NULL, use an empty string ("") in the data frame (regardless of the database column type). 
Also, when inserting data into a study dataset, the sequence number must be specified..
}

\value{
A list is returned with named categories of \bold{command}, \bold{rowsAffected}, \bold{rows}, \bold{queryName}, \bold{containerPath} and \bold{schemaName}.
The \bold{schemaName}, \bold{queryName} and \bold{containerPath} properties contain the same schema, query 
and folder path used in the request.  The
\bold{rowsAffected} property indicates the number of rows affected by the API action. This will typically be the same
number as passed in the request. The \bold{rows} property contains a list of row objects corresponding to the rows 
inserted.
}
\author{Valerie Obenchain}
\seealso{
\code{\link{labkey.selectRows}}, \code{\link{labkey.executeSql}}, \code{\link{makeFilter}}, 
\code{\link{labkey.importRows}}, \code{\link{labkey.updateRows}}, \cr
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
#	, AttachmentFld = NA    #attachment fields not supported 
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

deleterow <- data.frame(RowId=newRowId, stringsAsFactors=FALSE)
result <- labkey.deleteRows(baseUrl="http://localhost:8080/labkey",
    folderPath="/apisamples", schemaName="lists", queryName="AllTypes",
    toDelete=deleterow)

## Example of creating a provenance run with an initial step with material inputs, a second step
## with provenance mapping to link existing samples with newly inserted samples, and a final step
## with a data output
##
mi <- data.frame(lsid=c("urn:lsid:labkey.com:Sample.251.MySamples:sample1",
        "urn:lsid:labkey.com:Sample.251.MySamples:sample2"))
p <- labkey.provenance.createProvenanceParams(name="step1", description="initial step",
        materialInputs=mi)
ra <- labkey.provenance.startRecording(baseUrl="https://labkey.org/labkey/",
        folderPath = "Provenance", provenanceParams=p)

rows <- fromJSON(txt='[{
        "name" : "sample3",
        "protein" : "p3",
        "prov:objectInputs" : [
            "urn:lsid:labkey.com:Sample.251.MySamples:sample21",
            "urn:lsid:labkey.com:Sample.251.MySamples:sample22"
        ]
    },{
        "name" : "sample4",
        "protein" : "p4",
        "prov:objectInputs" : [
            "urn:lsid:labkey.com:Sample.251.MySamples:sample21",
            "urn:lsid:labkey.com:Sample.251.MySamples:sample22"
        ]
    }
]')

labkey.insertRows(baseUrl="https://labkey.org/labkey/", folderPath = "Provenance",
    schemaName="samples", queryName="MySamples", toInsert=rows,
    provenanceParams=labkey.provenance.createProvenanceParams(name="query step",
        recordingId=ra$recordingId))
labkey.provenance.stopRecording(baseUrl="https://labkey.org/labkey/", folderPath = "Provenance",
    provenanceParams=labkey.provenance.createProvenanceParams(name="final step",
    recordingId=ra$recordingId, dataOutputs=do))
}
}
\keyword{IO}