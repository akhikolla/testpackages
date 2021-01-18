\name{labkey.security.getContainers}
\alias{labkey.security.getContainers}
\title{Returns information about the specified container}
\description{
Returns information about the specified container, including the user's current permissions within that container.
If the includeSubfolders config option is set to true, it will also return information about all descendants
the user is allowed to see.
}
\usage{
labkey.security.getContainers(baseUrl=NULL, folderPath, includeEffectivePermissions=TRUE,
    includeSubfolders=FALSE, depth=50)
}
\arguments{
  \item{baseUrl}{A string specifying the \code{baseUrl} for the labkey server.}
  \item{folderPath}{A string specifying the \code{folderPath}.}
  \item{includeEffectivePermissions}{If set to false, the effective permissions for this container resource will not be included (defaults to true).}
  \item{includeSubfolders}{If set to true, the entire branch of containers will be returned. If false, only the immediate children of the starting container will be returned (defaults to false).}
  \item{depth}{May be used to control the depth of recursion if includeSubfolders is set to true.}
}
\details{
This function returns information about the specified container, including the user's current permissions within that container.
If the includeSubfolders config option is set to true, it will also return information about all descendants
the user is allowed to see. The depth of the results for the included subfolders can be controlled with the depth parameter.
}
\value{
The data frame containing the container properties for the current folder and subfolders, including name, title, id,
path, type, folderType, and effectivePermissions.
}

\author{Cory Nathe}
\seealso{
\code{\link{labkey.getFolders}},
\code{\link{labkey.security.createContainer}},
\code{\link{labkey.security.deleteContainer}},
\code{\link{labkey.security.moveContainer}}
}
\examples{
\donttest{

library(Rlabkey)

labkey.security.getContainers(
    baseUrl="http://labkey/", folderPath = "home",
    includeEffectivePermissions = FALSE, includeSubfolders = TRUE, depth = 2
)

}
}
\keyword{IO}
