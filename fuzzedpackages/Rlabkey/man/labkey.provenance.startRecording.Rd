\name{labkey.provenance.startRecording}
\alias{labkey.provenance.startRecording}
\title{Start a provenance recording}
\description{
Function to start a provenance recording session, if successful a provenance recording ID is returned which
can be used to add additional steps to the provenance run.
Note: this function is in beta and not yet final, changes should be expected so exercise caution when using it.
}
\usage{
labkey.provenance.startRecording(baseUrl=NULL, folderPath, provenanceParams = NULL)
}
\arguments{
  \item{baseUrl}{a string specifying the \code{baseUrl} for the labkey server}
  \item{folderPath}{a string specifying the \code{folderPath} }
  \item{provenanceParams}{the provenance parameter object which contains the options to include in this recording step}
}
\details{
Function to start a provenance recording. A provenance recording can contain an arbitrary number of steps to create a
provenance run, but stopRecording must be called to finish the recording and create the run. If successful
this will return a recording ID which is needed for subsequent steps.
This is a premium feature and requires the Provenance LabKey module to function correctly.
}
\value{
The generated recording ID which can be used in subsequent steps (or queries that support provenance).
}

\author{Karl Lum}
\seealso{
\code{\link{labkey.provenance.createProvenanceParams}},
\code{\link{labkey.provenance.addRecordingStep}},
\code{\link{labkey.provenance.stopRecording}}
}
\examples{
\donttest{

## create provenance params with material inputs and data outputs and start a recording
library(Rlabkey)

mi <- data.frame(lsid=c("urn:lsid:labkey.com:Sample.251.MySamples:sample1",
        "urn:lsid:labkey.com:Sample.251.MySamples:sample2"))
do <- data.frame(
        lsid="urn:lsid:labkey.com:AssayRunTSVData.Folder-251:12c70994-7ce5-1038-82f0-9c1487dbd334")

p <- labkey.provenance.createProvenanceParams(name="step1", description="initial step",
        materialInputs=mi, dataOutputs=do)
labkey.provenance.startRecording(baseUrl="https://labkey.org/labkey/",
        folderPath = "Provenance", provenanceParams=p)
}
}
\keyword{IO}
