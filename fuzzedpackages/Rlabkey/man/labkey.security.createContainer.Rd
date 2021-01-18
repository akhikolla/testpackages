\name{labkey.security.createContainer}
\alias{labkey.security.createContainer}
\title{Creates a new container, which may be a project, folder, or workbook, on the server}
\description{
Create a new container, which may be a project, folder, or workbook, on the LabKey server with parameters to control
the containers name, title, description, and folder type.
}
\usage{
labkey.security.createContainer(baseUrl=NULL, parentPath, name = NULL, title = NULL,
    description = NULL, folderType = NULL, isWorkbook = FALSE)
}
\arguments{
  \item{baseUrl}{A string specifying the \code{baseUrl} for the labkey server.}
  \item{parentPath}{A string specifying the \code{parentPath} for the new container.}
  \item{name}{The name of the container, required for projects or folders.}
  \item{title}{The title of the container, used primarily for workbooks.}
  \item{description}{The description of the container, used primarily for workbooks.}
  \item{folderType}{The name of the folder type to be applied (ex. Study or Collaboration).}
  \item{isWorkbook}{Whether this a workbook should be created. Defaults to false.}
}
\details{
This function allows for users with proper permissions to create a new container, which may be a project, folder, or
workbook, on the LabKey server with parameters to control the containers name, title, description, and folder type.
If the container already exists or the user does not have permissions, an error message will be returned.
}
\value{
Returns information about the newly created container.
}

\author{Cory Nathe}
\seealso{
\code{\link{labkey.getFolders}},
\code{\link{labkey.security.getContainers}},
\code{\link{labkey.security.deleteContainer}},
\code{\link{labkey.security.moveContainer}}
}
\examples{
\donttest{

library(Rlabkey)

labkey.security.createContainer(baseUrl="http://labkey/", parentPath = "/home",
    name = "NewFolder", description = "My new folder has this description",
    folderType = "Collaboration"
)

}
}
\keyword{IO}
