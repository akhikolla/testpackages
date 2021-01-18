\name{labkey.provenance.createProvenanceParams}
\alias{labkey.provenance.createProvenanceParams}
\title{Create provenance parameter object}
\description{
Helper function to create the data structure that can be used in provenance related APIs.
Note: this function is in beta and not yet final, changes should be expected so exercise caution when using it.
}
\usage{
labkey.provenance.createProvenanceParams(recordingId=NULL, name=NULL, description=NULL,
        materialInputs=NULL, materialOutputs=NULL, dataInputs=NULL, dataOutputs=NULL,
        inputObjectUriProperty=NULL, outputObjectUriProperty=NULL, objectInputs=NULL,
        objectOutputs=NULL, provenanceMap=NULL)
}
\arguments{
  \item{recordingId}{(optional) the recording ID to associate with other steps using the same ID}
  \item{name}{(optional) the name of this provenance step}
  \item{description}{(optional) the description of this provenance step}
  \item{materialInputs}{(optional) the list of materials (samples) to be used as the provenance run input. The data structure should be a dataframe with the column name describing the data type (lsid, id)}
  \item{materialOutputs}{(optional) the list of materials (samples) to be used as the provenance run output. The data structure should be a dataframe with the column name describing the data type (lsid, id)}
  \item{dataInputs}{(optional) the list of data inputs to be used for the run provenance map}
  \item{dataOutputs}{(optional) the list of data outputs to be used for the run provenance map}
  \item{inputObjectUriProperty}{(optional) for incoming data rows, the column name to interpret as the input to the provenance map. Defaults to : 'prov:objectInputs'}
  \item{outputObjectUriProperty}{(optional) for provenance mapping, the column name to interpret as the output to the provenance map. Defaults to : 'lsid'}
  \item{objectInputs}{(optional) the list of object inputs to be used for the run provenance map}
  \item{objectOutputs}{(optional) the list of object outputs to be used for the run provenance map}
  \item{provenanceMap}{(optional) the provenance map to be used directly for the run step. The data structure should be a dataframe with the column names of 'from' and 'to' to indicate which sides of the mapping the identifiers refer to}
}
\details{
This function can be used to generate a provenance parameter object which can then be used as an argument
in the other provenance related functions to assemble provenance runs.
This is a premium feature and requires the Provenance LabKey module to function correctly.
}
\value{
A list containing elements describing the passed in provenance parameters.
}

\author{Karl Lum}
\seealso{
\code{\link{labkey.provenance.startRecording}},
\code{\link{labkey.provenance.addRecordingStep}},
\code{\link{labkey.provenance.stopRecording}}
}
\examples{
\donttest{

## create provenance params with material inputs and data outputs
library(Rlabkey)

mi <- data.frame(lsid=c("urn:lsid:labkey.com:Sample.251.MySamples:sample1",
        "urn:lsid:labkey.com:Sample.251.MySamples:sample2"))
do <- data.frame(
        lsid="urn:lsid:labkey.com:AssayRunTSVData.Folder-251:12c70994-7ce5-1038-82f0-9c1487dbd334")

p <- labkey.provenance.createProvenanceParams(name="step1", description="initial step",
        materialInputs=mi, dataOutputs=do)

## create provenance params with object inputs (from an assay run)
oi <- labkey.selectRows(baseUrl="https://labkey.org/labkey/", folderPath = "Provenance",
        schemaName="assay.General.titer",
        queryName="Data",
        colSelect= c("LSID"),
        colFilter=makeFilter(c("Run/RowId","EQUAL","253")))
mi <- data.frame(lsid=c("urn:lsid:labkey.com:Sample.251.MySamples:sample1",
        "urn:lsid:labkey.com:Sample.251.MySamples:sample2"))

p <- labkey.provenance.createProvenanceParams(name="step1", description="initial step",
        objectInputs=oi[["LSID"]], materialInputs=mi)

}
}
\keyword{IO}
