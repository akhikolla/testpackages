\name{labkey.security.deleteContainer}
\alias{labkey.security.deleteContainer}
\title{Deletes an existing container, which may be a project, folder, or workbook}
\description{
Deletes an existing container, which may be a project, folder, or workbook, and all of its children from the Labeky server.
}
\usage{
labkey.security.deleteContainer(baseUrl=NULL, folderPath)
}
\arguments{
  \item{baseUrl}{A string specifying the \code{baseUrl} for the labkey server.}
  \item{folderPath}{A string specifying the \code{folderPath} to be deleted.}
}
\details{
This function allows for users with proper permissions to delete an existing container, which may be a project, folder,
or workbook, from the LabKey server. This will also remove all subfolders of the container being deleted.
If the container does not exist or the user does not have permissions, an error message will be returned.
}
\value{
Returns a success message for the container deletion action.
}

\author{Cory Nathe}
\seealso{
\code{\link{labkey.getFolders}},
\code{\link{labkey.security.getContainers}},
\code{\link{labkey.security.createContainer}},
\code{\link{labkey.security.moveContainer}}
}
\examples{
\donttest{

library(Rlabkey)

labkey.security.deleteContainer(baseUrl="http://labkey/", folderPath = "/home/FolderToDelete")

}
}
\keyword{IO}
