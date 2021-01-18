\name{lsProjects}
\alias{lsProjects}
\title{List the projects available at a given LabKey Server address }
\description{ Lists the projects available. Takes a string URL instead of a session, as it is intended for use before creating a session.}

\usage{
lsProjects(baseUrl) 
}
\arguments{
  \item{baseUrl}{a string specifying the \code{baseUrl}for the LabKey Server, of the form http://<server dns name>/<contextroot> }

}
\details{  List the projects available at a given LabKey Server address.
}
\value{ A character array containing the available projects, relative to the root. These values can be set on a session using curFolder<-
}
\references{https://www.labkey.org/project/home/begin.view}
\author{Peter Hussey}
\seealso{\code{\link{getSession}}, \code{\link{lsFolders}}, \code{\link{lsSchemas}}}
\examples{
\donttest{

## get list of projects on server, connect a session in one project,
## then list the folders in that project
# library(Rlabkey)

lsProjects("https://www.labkey.org")

lkorg <- getSession("https://www.labkey.org", "/home")
lsFolders(lkorg)

lkorg <- getSession("https://www.labkey.org", "/home/Study/ListDemo")
lsSchemas(lkorg)

}
}
\keyword{file}
