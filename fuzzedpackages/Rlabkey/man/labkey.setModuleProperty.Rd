\name{labkey.setModuleProperty}
\alias{labkey.setModuleProperty}
\title{Set module property value}
\description{
Set module property value for a specific folder or as site wide (with folderPath '/')
}
\usage{
labkey.setModuleProperty(baseUrl=NULL, folderPath, moduleName, propName, propValue)
}
\arguments{
\item{baseUrl}{server location including context path, if any. e.g. https://www.labkey.org/}

\item{folderPath}{a string specifying the \code{folderPath}}

\item{moduleName}{name of the module}

\item{propName}{The module property name}

\item{propValue}{The module property value to save}
}
\examples{
\donttest{
library(Rlabkey)
labkey.setModuleProperty(baseUrl="http://labkey/", folderPath="flowProject",
    moduleName="flow", propName="ExportToScriptFormat", propValue="zip")

}
}
