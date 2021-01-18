\name{labkey.webdav.mkDirs}
\alias{labkey.webdav.mkDirs}
\title{Create a folder via WebDAV}
\description{
This will create folder(s) under the specified LabKey Server project using WebDAV.
}
\usage{
labkey.webdav.mkDirs(
    baseUrl=NULL,
    folderPath,
    remoteFilePath,
    fileSet='@files'
    )
}
\arguments{
  \item{baseUrl}{a string specifying the \code{baseUrl} for the labkey server}
  \item{folderPath}{a string specifying the \code{folderPath} }
  \item{remoteFilePath}{the folder path to create, relative to the LabKey folder root. }
  \item{fileSet}{(optional) the name of file server fileSet, which is typically "@files" (the default value for this argument). In some cases this might be "@pipeline" or "@fileset". }
}
\details{
Creates a folder on a LabKey Server using WebDAV. If the parent directory or directories no not exist, these will also be created (similar to mkdir -p on linux)
}
\value{
TRUE if the folder was created successfully
}

\author{Ben Bimber, Ph.D.}
\seealso{
\code{\link{labkey.webdav.get}},
\code{\link{labkey.webdav.put}},
\code{\link{labkey.webdav.mkDir}},
\code{\link{labkey.webdav.pathExists}},
\code{\link{labkey.webdav.listDir}},
\code{\link{labkey.webdav.delete}},
\code{\link{labkey.webdav.downloadFolder}}
}
\examples{
\donttest{

library(Rlabkey)

labkey.webdav.mkDirs(
    baseUrl="http://labkey/",
    folderPath="home",
    remoteFilePath="folder1/folder2/toCreate"
)

}
}
\keyword{IO}
