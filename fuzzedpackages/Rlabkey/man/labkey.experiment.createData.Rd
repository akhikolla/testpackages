\name{labkey.experiment.createData}
\alias{labkey.experiment.createData}
\title{Create an experiment data object}
\description{
Create an experiment data object.
}
\usage{
labkey.experiment.createData(config,
    dataClassId = NULL, dataClassName = NULL, dataFileUrl = NULL)
}
\arguments{
  \item{config}{a list of base experiment object properties}
  \item{dataClassId}{(optional) a integer specifying the data class row ID }
  \item{dataClassName}{(optional) a string specifying the name of the data class}
  \item{dataFileUrl}{(optional) a string specifying the local file url of the uploaded file}
}
\details{
Create an experiment data object which can be used as either input or output datas for an experiment run.
}

\value{
 Returns the object representation of the experiment data object.
}

\author{Karl Lum}
\seealso{
\code{\link{labkey.experiment.saveBatch}},
\code{\link{labkey.experiment.createMaterial}},
\code{\link{labkey.experiment.createRun}}
}
\examples{
\donttest{

library(Rlabkey)

## create a non-assay backed run with data classes as data inputs and outputs

d1 <- labkey.experiment.createData(
        list(name = "dc-01"), dataClassId = 400)
d2 <- labkey.experiment.createData(
        list(name = "dc-02"), dataClassId = 402)
run <- labkey.experiment.createRun(
        list(name="new run"), dataInputs = d1, dataOutputs = d2)
labkey.experiment.saveBatch(baseUrl="http://labkey/", folderPath="home",
        protocolName=labkey.experiment.SAMPLE_DERIVATION_PROTOCOL, runList=run)
}
}
\keyword{IO}
