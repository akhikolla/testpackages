\name{labkey.webdav.get}
\alias{labkey.webdav.get}
\title{Download a file via WebDAV}
\description{
This will download a file from a LabKey Server using WebDAV.
}
\usage{
labkey.webdav.get(
    baseUrl=NULL,
    folderPath,
    remoteFilePath,
    localFilePath,
    overwrite=TRUE,
    fileSet='@files'
    )
}
\arguments{
  \item{baseUrl}{a string specifying the \code{baseUrl} for the labkey server}
  \item{folderPath}{a string specifying the \code{folderPath} }
  \item{remoteFilePath}{the path of this file on the remote server, relative to the folder root. }
  \item{localFilePath}{the local filepath where this file will be saved }
  \item{overwrite}{(optional) if true, any pre-existing file at this location will be overwritten.  Defaults to TRUE }
  \item{fileSet}{(optional) the name of file server fileSet, which is typically "@files" (the default value for this argument). In some cases this might be "@pipeline" or "@fileset". }
}
\details{
Download a single file from a LabKey Server to the local machine using WebDAV.
}
\value{
TRUE or FALSE, depending on if this file was downloaded and exists locally.  Will return FALSE if the already file exists and overwrite=F.
}

\author{Ben Bimber, Ph.D.}
\seealso{
\code{\link{labkey.webdav.put}},
\code{\link{labkey.webdav.mkDir}},
\code{\link{labkey.webdav.mkDirs}},
\code{\link{labkey.webdav.pathExists}},
\code{\link{labkey.webdav.listDir}},
\code{\link{labkey.webdav.delete}},
\code{\link{labkey.webdav.downloadFolder}}
}
\examples{
\donttest{
## download a single file from a LabKey Server
library(Rlabkey)

labkey.webdav.get(
    baseUrl="http://labkey/",
    folderPath="home",
    remoteFilePath="folder/myFile.txt",
    localFilePath="myDownloadedFile.txt",
    overwrite=TRUE
)

}
}
\keyword{IO}
