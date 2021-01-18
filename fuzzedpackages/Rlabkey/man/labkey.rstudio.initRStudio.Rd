\name{labkey.rstudio.initRStudio}
\alias{labkey.rstudio.initRStudio}
\title{Initialize a RStudio session for LabKey integration}
\description{
LabKey-RStudio integration helper. Not intended for use outside RStudio.
}
\usage{
labkey.rstudio.initRStudio(apiKey = "", baseUrl = "", folderPath, skipViewer = FALSE)
}
\arguments{
\item{apiKey}{session key from your server}

\item{baseUrl}{server location including context path, if any. e.g. https://www.labkey.org/}

\item{folderPath}{a string specifying the \code{folderPath}}

\item{skipViewer}{(TRUE | FALSE) TRUE to skip setting up LabKey schema viewer in RStudio}
}
\examples{
\donttest{
## RStudio console only
library(Rlabkey)
labkey.rstudio.initRStudio(apiKey="session|abcdef0123456789abcdef0123456789",
    baseUrl="http://labkey/", folderPath="home")

}
}
