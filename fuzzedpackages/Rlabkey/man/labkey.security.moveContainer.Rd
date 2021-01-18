\name{labkey.security.moveContainer}
\alias{labkey.security.moveContainer}
\title{Moves an existing container, which may be a folder or workbook}
\description{
Moves an existing container, which may be a folder or workbook, to be the subfolder of another folder and/or project
on the LabKey server.
}
\usage{
labkey.security.moveContainer(baseUrl=NULL, folderPath,
    destinationParent, addAlias = TRUE)
}
\arguments{
  \item{baseUrl}{A string specifying the \code{baseUrl} for the labkey server.}
  \item{folderPath}{A string specifying the \code{folderPath} to be moved. Additionally, the container entity id is also valid.}
  \item{destinationParent}{The container path of destination parent. Additionally, the destination parent entity id is also valid.}
  \item{addAlias}{Add alias of current container path to container that is being moved (defaults to true).}
}
\details{
This function moves an existing container, which may be a folder or workbook, to be the subfolder of another
folder and/or project on the LabKey server. Projects and the root container can not be moved.
If the target or destination container does not exist or the user does not have permissions, an error message will
be returned.
}
\value{
Returns a success message for the container move action with the new path.
}

\author{Cory Nathe}
\seealso{
\code{\link{labkey.getFolders}},
\code{\link{labkey.security.getContainers}},
\code{\link{labkey.security.createContainer}},
\code{\link{labkey.security.deleteContainer}}
}
\examples{
\donttest{

library(Rlabkey)

labkey.security.moveContainer(baseUrl="http://labkey/", folderPath = "/home/FolderToMove",
    destinationParent = "/OtherProject", addAlias = TRUE
)

}
}
\keyword{IO}
