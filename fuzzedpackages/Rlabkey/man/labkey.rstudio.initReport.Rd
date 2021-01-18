\name{labkey.rstudio.initReport}
\alias{labkey.rstudio.initReport}
\title{Initialize a RStudio session for LabKey R report source editing}
\description{
LabKey-RStudio integration helper. Not intended for use outside RStudio.
}
\usage{
labkey.rstudio.initReport(apiKey = "", baseUrl = "", folderPath,
  reportEntityId, skipViewer = FALSE, skipEdit = FALSE)
}
\arguments{
\item{apiKey}{session key from your server}

\item{baseUrl}{server location including context path, if any. e.g. https://www.labkey.org/}

\item{folderPath}{a string specifying the \code{folderPath}}

\item{reportEntityId}{LabKey report's entityId}

\item{skipViewer}{(TRUE | FALSE) TRUE to skip setting up LabKey schema viewer in RStudio}

\item{skipEdit}{(TRUE | FALSE) TRUE to open file in editor}
}
\examples{
\donttest{
## RStudio console only
library(Rlabkey)
labkey.rstudio.initReport(apiKey="session|abcdef0123456789abcdef0123456789",
    baseUrl="http://labkey/", folderPath="home",
    reportEntityId="0123456a-789b-1000-abcd-01234567abcde")

}
}
