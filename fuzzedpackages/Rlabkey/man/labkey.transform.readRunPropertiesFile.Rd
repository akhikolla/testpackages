\name{labkey.transform.readRunPropertiesFile}
\alias{labkey.transform.readRunPropertiesFile}
\title{Assay transform script helper function to read a run properties file}
\description{
A function that takes in the full path to the LabKey generated run properties file and returns a data.frame of
the key value pairs for the lines within that file. This helper function would be used as part of an assay transform
script written in R and associated with an assay design.}
\usage{
labkey.transform.readRunPropertiesFile(runInfoPath)
}
\arguments{
  \item{runInfoPath}{the full file system path to the generated run properties file}
}

\details{
The most common scenario is that the assay transform script will get the run properties file path added
into the running script as a replacement variable. To use that replacement variable for this helper function,
you can pass in the runInfoPath parameter as "${runInfo}".\cr
}
\examples{
\donttest{

# library(Rlabkey)

labkey.transform.readRunPropertiesFile("${runInfo}")

}
}
