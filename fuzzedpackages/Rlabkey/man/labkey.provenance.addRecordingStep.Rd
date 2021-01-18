\name{labkey.provenance.addRecordingStep}
\alias{labkey.provenance.addRecordingStep}
\title{Add a step to a provenance recording}
\description{
Function to add a step to a previously created provenance recording session.
Note: this function is in beta and not yet final, changes should be expected so exercise caution when using it.
}
\usage{
labkey.provenance.addRecordingStep(baseUrl=NULL, folderPath, provenanceParams = NULL)
}
\arguments{
  \item{baseUrl}{a string specifying the \code{baseUrl} for the labkey server}
  \item{folderPath}{a string specifying the \code{folderPath} }
  \item{provenanceParams}{the provenance parameter object which contains the options to include in this recording step}
}
\details{
Function to add a step to a previously created provenance recording. The recording ID that was obtained from
a previous startRecording function call must be passed into the provenanceParams config.
This is a premium feature and requires the Provenance LabKey module to function correctly.
}
\value{
The generated recording ID which can be used in subsequent steps (or queries that support provenance).
}

\author{Karl Lum}
\seealso{
\code{\link{labkey.provenance.createProvenanceParams}},
\code{\link{labkey.provenance.startRecording}},
\code{\link{labkey.provenance.stopRecording}}
}
\examples{
\donttest{

## start a provenance recording and add a recording step
library(Rlabkey)

mi <- data.frame(lsid=c("urn:lsid:labkey.com:Sample.251.MySamples:sample1",
        "urn:lsid:labkey.com:Sample.251.MySamples:sample2"))

p <- labkey.provenance.createProvenanceParams(name="step1", description="initial step",
        materialInputs=mi)
r <- labkey.provenance.startRecording(baseUrl="https://labkey.org/labkey/",
        folderPath = "Provenance", provenanceParams=p)
do <- data.frame(
        lsid="urn:lsid:labkey.com:AssayRunTSVData.Folder-251:12c70994-7ce5-1038-82f0-9c1487dbd334")

labkey.provenance.addRecordingStep(baseUrl="https://labkey.org/labkey/", folderPath = "Provenance",
    provenanceParams=labkey.provenance.createProvenanceParams(name="additional step",
        recordingId=r$recordingId, dataOutputs=do))
}
}
\keyword{IO}
