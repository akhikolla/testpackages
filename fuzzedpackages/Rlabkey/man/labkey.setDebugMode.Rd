\name{labkey.setDebugMode}
\alias{labkey.setDebugMode}
\title{Helper function to enable/disable debug mode.}
\description{
When debug mode is enabled, the GET/POST calls with output information about the request being made
and will output a raw string version of the response object.
}
\usage{
labkey.setDebugMode(debug = FALSE)
}
\arguments{
  \item{debug}{a boolean specifying if debug mode is enabled or disabled}
}

\author{Cory Nathe}
\examples{
\donttest{

library(Rlabkey)
labkey.setDebugMode(TRUE)
labkey.executeSql(
    baseUrl="http://localhost:8080/labkey",
    folderPath="/home",
    schemaName="core",
    sql = "select * from containers")

}
}
\keyword{IO}
