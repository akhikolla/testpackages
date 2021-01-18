\name{labkey.experiment.createMaterial}
\alias{labkey.experiment.createMaterial}
\title{Create an experiment material object}
\description{
Create an experiment material object.
}
\usage{
labkey.experiment.createMaterial(config, sampleSetId = NULL, sampleSetName = NULL)
}
\arguments{
  \item{config}{a list of base experiment object properties}
  \item{sampleSetId}{(optional) a integer specifying the sample set row ID }
  \item{sampleSetName}{(optional) a string specifying the name of the sample set}
}
\details{
Create an experiment material object which can be used as either input or output materials for an experiment run.
}

\value{
 Returns the object representation of the experiment material object.
}

\author{Karl Lum}
\seealso{
\code{\link{labkey.experiment.saveBatch}},
\code{\link{labkey.experiment.createData}},
\code{\link{labkey.experiment.createRun}}
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
}
}
\keyword{IO}
