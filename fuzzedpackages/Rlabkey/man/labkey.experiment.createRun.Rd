\name{labkey.experiment.createRun}
\alias{labkey.experiment.createRun}
\title{Create an experiment run object}
\description{
Create an experiment run object.
}
\usage{
labkey.experiment.createRun(config,
    dataInputs = NULL, dataOutputs = NULL, dataRows = NULL,
    materialInputs = NULL, materialOutputs = NULL, plateMetadata = NULL)
}
\arguments{
  \item{config}{a list of base experiment object properties}
  \item{dataInputs}{(optional) a list of experiment data objects to be used as data inputs to the run}
  \item{dataOutputs}{(optional) a list of experiment data objects to be used as data outputs to the run}
  \item{dataRows}{(optional) a data frame containing data rows to be uploaded to the assay backed run}
  \item{materialInputs}{(optional) a list of experiment material objects to be used as material inputs to the run}
  \item{materialOutputs}{(optional) a list of experiment material objects to be used as material outputs to the run}
  \item{plateMetadata}{(optional) if the assay supports plate templates, the plate metadata object to upload}
}
\details{
Create an experiment run object which can be used in the saveBatch function.
}

\value{
 Returns the object representation of the experiment run object.
}

\author{Karl Lum}
\seealso{
\code{\link{labkey.experiment.saveBatch}},
\code{\link{labkey.experiment.createData}},
\code{\link{labkey.experiment.createMaterial}}
}
\examples{
\donttest{

library(Rlabkey)

## create a non-assay backed run with samples as material inputs and outputs

m1 <- labkey.experiment.createMaterial(
        list(name = "87444063.2604.626"), sampleSetName = "Study Specimens")
m2 <- labkey.experiment.createMaterial(
        list(name = "87444063.2604.625"), sampleSetName = "Study Specimens")
run <- labkey.experiment.createRun(
        list(name="new run"), materialInputs = m1, materialOutputs = m2)
labkey.experiment.saveBatch(baseUrl="http://labkey/", folderPath="home",
        protocolName=labkey.experiment.SAMPLE_DERIVATION_PROTOCOL, runList=run)

## import an assay run which includes plate metadata

df <- data.frame(participantId=c(1:3), visitId = c(10,20,30), welllocation = c("A1", "D11", "F12"))

runConfig <- fromJSON(txt='{"assayId": 310,
    "name" : "api imported run with plate metadata",
    "properties" : {
        "PlateTemplate" : "urn:lsid:labkey.com:PlateTemplate.Folder-6:d8bbec7d-34cd-1038-bd67-b3bd"
    }}
}')

plateMetadata <- fromJSON(txt='{
      "control" : {
        "positive" : {"dilution":  0.005},
        "negative" : {"dilution":  1.0}
      },
      "sample" : {
        "SA01" : {"dilution": 1.0, "ID" : 111, "Barcode" : "BC_111", "Concentration" : 0.0125},
        "SA02" : {"dilution": 2.0, "ID" : 222, "Barcode" : "BC_222"},
        "SA03" : {"dilution": 3.0, "ID" : 333, "Barcode" : "BC_333"},
        "SA04" : {"dilution": 4.0, "ID" : 444, "Barcode" : "BC_444"}
      }
    }')

run <- labkey.experiment.createRun(runConfig, dataRows = df, plateMetadata = plateMetadata)
labkey.experiment.saveBatch(
    baseUrl="http://labkey/", folderPath="home",
    assayConfig=list(assayId = 310), runList=run
)

}
}
\keyword{IO}
