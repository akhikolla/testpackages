\name{labkey.webdav.listDir}
\alias{labkey.webdav.listDir}
\title{List the contents of a LabKey Server folder via WebDAV}
\description{
This will list the contents of a LabKey Server folder using WebDAV.
}
\usage{
labkey.webdav.listDir(
    baseUrl=NULL,
    folderPath,
    remoteFilePath,
    fileSet='@files',
    haltOnError=TRUE
    )
}
\arguments{
  \item{baseUrl}{a string specifying the \code{baseUrl} for the labkey server}
  \item{folderPath}{a string specifying the \code{folderPath} }
  \item{remoteFilePath}{path of the folder on the remote server, relative to the folder root. }
  \item{fileSet}{(optional) the name of file server fileSet, which is typically "@files" (the default value for this argument). In some cases this might be "@pipeline" or "@fileset". }
  \item{haltOnError}{(optional) Specifies whether this request should fail if the requested path does not exist.  Defaults to TRUE }
}
\details{
Lists the contents of a folder on a LabKey Server using WebDAV.
}
\value{
A list with each item under this folder.  Each item (file or directory) is a list with the following attributes:
     \itemize{
       	      \item{ "files": A list of the files, where each has the following attributes: }
              \itemize{
                  \item{ "id": The relative path to this item, not encoded }
                  \item{ "href": The relative URL to this item, HTML encoded }
                  \item{ "text": A dataset in a date based study }
                  \item{ "creationdate": The date this item was created }
                  \item{ "createdby": The user that created this file }
                  \item{ "lastmodified": The last modification time }
                  \item{ "contentlength": The content length }
                  \item{ "size": The file size }
                  \item{ "isdirectory": TRUE or FALSE, depending on whether this item is a directory }
              }
              \item{ "fileCount": If this item is a directory, this property will be present, listing the the total files in this location }
     }
}

\author{Ben Bimber, Ph.D.}
\seealso{
\code{\link{labkey.webdav.get}},
\code{\link{labkey.webdav.put}},
\code{\link{labkey.webdav.mkDir}},
\code{\link{labkey.webdav.mkDirs}},
\code{\link{labkey.webdav.pathExists}},
\code{\link{labkey.webdav.delete}},
\code{\link{labkey.webdav.downloadFolder}}
}
\examples{
\donttest{

library(Rlabkey)

json <- labkey.webdav.listDir(
            baseUrl="http://labkey/",
            folderPath="home",
            remoteFilePath="myFolder"
        )

}
}
\keyword{IO}
