\name{labkey.getModuleProperty}
\alias{labkey.getModuleProperty}
\title{Get effective module property value}
\description{
Get a specific effective module property value for folder
}
\usage{
labkey.getModuleProperty(baseUrl=NULL, folderPath, moduleName, propName)
}
\arguments{
\item{baseUrl}{server location including context path, if any. e.g. https://www.labkey.org/}

\item{folderPath}{a string specifying the \code{folderPath}}

\item{moduleName}{name of the module}

\item{propName}{The module property name}
}
\examples{
\donttest{
library(Rlabkey)
labkey.getModuleProperty(baseUrl="http://labkey/", folderPath="flowProject",
    moduleName="flow", propName="ExportToScriptFormat")

}
}
