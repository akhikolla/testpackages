\name{labkey.webdav.pathExists}
\alias{labkey.webdav.pathExists}
\title{Tests if a path exists on a LabKey Server via WebDAV}
\description{
This will test if the supplied file/folder exists folder under the specified LabKey Server project using WebDAV.
}
\usage{
labkey.webdav.pathExists(
    baseUrl=NULL,
    folderPath,
    remoteFilePath,
    fileSet='@files'
    )
}
\arguments{
  \item{baseUrl}{a string specifying the \code{baseUrl} for the labkey server}
  \item{folderPath}{a string specifying the \code{folderPath} }
  \item{remoteFilePath}{the path to test, relative to the LabKey folder root. }
  \item{fileSet}{(optional) the name of file server fileSet, which is typically "@files" (the default value for this argument). In some cases this might be "@pipeline" or "@fileset". }
}
\details{
This will test if the supplied file/folder exists folder under the specified LabKey Server project using WebDAV.
}
\value{
TRUE if the folder was created successfully
}

\author{Ben Bimber, Ph.D.}
\seealso{
\code{\link{labkey.webdav.get}},
\code{\link{labkey.webdav.put}},
\code{\link{labkey.webdav.mkDir}},
\code{\link{labkey.webdav.mkDirs}},
\code{\link{labkey.webdav.listDir}},
\code{\link{labkey.webdav.delete}},
\code{\link{labkey.webdav.downloadFolder}}
}
\examples{
\donttest{

library(Rlabkey)

# Test folder
labkey.webdav.pathExists(
    baseUrl="http://labkey/",
    folderPath="home",
    remoteFilePath="pathToTest"
)

# Test file
labkey.webdav.pathExists(
    baseUrl="http://labkey/",
    folderPath="home",
    remoteFilePath="folder/fileToTest.txt"
)

}
}
\keyword{IO}
