\name{labkey.provenance.stopRecording}
\alias{labkey.provenance.stopRecording}
\title{Stop a provenance recording}
\description{
Function to end a provenance recording and create and save the provenance run on the server.
Note: this function is in beta and not yet final, changes should be expected so exercise caution when using it.
}
\usage{
labkey.provenance.stopRecording(baseUrl=NULL, folderPath, provenanceParams = NULL)
}
\arguments{
  \item{baseUrl}{a string specifying the \code{baseUrl} for the labkey server}
  \item{folderPath}{a string specifying the \code{folderPath} }
  \item{provenanceParams}{the provenance parameter object which contains the options to include in this recording step}
}
\details{
Function to stop the provenance recording associated with the recording ID, this will create a provenance run
using all the steps (with inputs and outputs) associated with the recording ID. The recording ID that was obtained from
a previous startRecording function call must be passed into the provenanceParams config.
This is a premium feature and requires the Provenance LabKey module to function correctly.
}
\value{
The serialized provenance run that was created.
}

\author{Karl Lum}
\seealso{
\code{\link{labkey.provenance.createProvenanceParams}},
\code{\link{labkey.provenance.startRecording}},
\code{\link{labkey.provenance.addRecordingStep}}
}
\examples{
\donttest{

library(Rlabkey)

## object inputs (from an assay run) and material inputs
##
oi <- labkey.selectRows(baseUrl="https://labkey.org/labkey/", folderPath = "Provenance",
        schemaName="assay.General.titer",
        queryName="Data",
        colSelect= c("LSID"),
        colFilter=makeFilter(c("Run/RowId","EQUAL","253")))
mi <- data.frame(lsid=c("urn:lsid:labkey.com:Sample.251.MySamples:sample1",
        "urn:lsid:labkey.com:Sample.251.MySamples:sample2"))

p <- labkey.provenance.createProvenanceParams(name="step1", description="initial step",
        objectInputs=oi[["LSID"]], materialInputs=mi)
r <- labkey.provenance.startRecording(baseUrl="https://labkey.org/labkey/",
        folderPath = "Provenance", provenanceParams=p)
run <- labkey.provenance.stopRecording(baseUrl="https://labkey.org/labkey/",
        folderPath = "Provenance",
        provenanceParams=labkey.provenance.createProvenanceParams(name="final step",
            recordingId=r$recordingId))
}
}
\keyword{IO}
