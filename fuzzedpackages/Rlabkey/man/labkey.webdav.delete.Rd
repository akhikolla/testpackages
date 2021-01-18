\name{labkey.webdav.delete}
\alias{labkey.webdav.delete}
\title{Deletes the provided file/folder on a LabKey Server via WebDAV}
\description{
This will delete the supplied file or folder under the specified LabKey Server project using WebDAV.
}
\usage{
labkey.webdav.delete(
    baseUrl=NULL,
    folderPath,
    remoteFilePath,
    fileSet='@files'
    )
}
\arguments{
  \item{baseUrl}{a string specifying the \code{baseUrl} for the labkey server}
  \item{folderPath}{a string specifying the \code{folderPath} }
  \item{remoteFilePath}{the path to delete, relative to the LabKey folder root. }
  \item{fileSet}{(optional) the name of file server fileSet, which is typically "@files" (the default value for this argument). In some cases this might be "@pipeline" or "@fileset". }
}
\details{
This will delete the supplied file or folder under the specified LabKey Server project using WebDAV. Note: if a folder is provided, it will delete that folder and contents.
}
\value{
TRUE if the folder was deleted successfully
}

\author{Ben Bimber, Ph.D.}
\seealso{
\code{\link{labkey.webdav.get}},
\code{\link{labkey.webdav.put}},
\code{\link{labkey.webdav.mkDir}},
\code{\link{labkey.webdav.mkDirs}},
\code{\link{labkey.webdav.listDir}},
\code{\link{labkey.webdav.pathExists}},
\code{\link{labkey.webdav.downloadFolder}}
}
\examples{
\donttest{

library(Rlabkey)

#delete an entire directory and contents
labkey.webdav.delete(baseUrl="http://labkey/", folderPath="home", remoteFilePath="folder1")

#delete single file
labkey.webdav.delete(baseUrl="http://labkey/", folderPath="home", remoteFilePath="folder/file.txt")

}
}
\keyword{IO}
