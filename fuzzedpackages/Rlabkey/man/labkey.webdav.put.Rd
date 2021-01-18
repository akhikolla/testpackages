\name{labkey.webdav.put}
\alias{labkey.webdav.put}
\title{Upload a file via WebDAV}
\description{
This will upload a file to a LabKey Server using WebDAV.
}
\usage{
labkey.webdav.put(
    localFile,
    baseUrl=NULL,
    folderPath,
    remoteFilePath,
    fileSet='@files'
    )
}
\arguments{
  \item{localFile}{the local filepath to upload }
  \item{baseUrl}{a string specifying the \code{baseUrl} for the labkey server}
  \item{folderPath}{a string specifying the \code{folderPath} }
  \item{remoteFilePath}{the destination path of this file on the remote server, relative to the folder root. }
  \item{fileSet}{(optional) the name of file server fileSet, which is typically "@files" (the default value for this argument). In some cases this might be "@pipeline" or "@fileset". }
}
\details{
Upload a single file from the local machine to a LabKey Server using WebDAV.
}
\value{
TRUE if the file was uploaded successfully
}

\author{Ben Bimber, Ph.D.}
\seealso{
\code{\link{labkey.webdav.get}},
\code{\link{labkey.webdav.mkDir}},
\code{\link{labkey.webdav.mkDirs}},
\code{\link{labkey.webdav.pathExists}},
\code{\link{labkey.webdav.listDir}},
\code{\link{labkey.webdav.delete}},
\code{\link{labkey.webdav.downloadFolder}}
}
\examples{
\donttest{
## upload a single file to a LabKey Server
library(Rlabkey)

labkey.webdav.put(
    localFile="myFileToUpload.txt",
    baseUrl="http://labkey/",
    folderPath="home",
    remoteFilePath="myFileToUpload.txt"
)

}
}
\keyword{IO}
