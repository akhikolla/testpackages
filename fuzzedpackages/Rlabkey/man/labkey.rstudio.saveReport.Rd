\name{labkey.rstudio.saveReport}
\alias{labkey.rstudio.saveReport}
\title{Update RStudio report source back to LabKey}
\description{
LabKey-RStudio integration helper. Not intended for use outside RStudio.
}
\usage{
labkey.rstudio.saveReport(folderPath, reportEntityId, reportFilename,
  useWarning = FALSE)
}
\arguments{
\item{folderPath}{a string specifying the \code{folderPath}}

\item{reportEntityId}{LabKey report's entityId}

\item{reportFilename}{The filename to save}

\item{useWarning}{(TRUE | FALSE) TRUE to prompt user choices to save}
}
\examples{
\donttest{
## RStudio console only
library(Rlabkey)
labkey.rstudio.saveReport(folderPath="home",
    reportEntityId="0123456a-789b-1000-abcd-01234567abcde",
    reportFilename="knitrReport.Rhtml", useWarning=TRUE)

}
}
