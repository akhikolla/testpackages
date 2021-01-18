\name{lsFolders}
\alias{lsFolders}
\title{List the available folder paths }
\description{ Lists the available folder paths relative to the current folder path for a LabKey session }

\usage{
lsFolders(session) 
}
\arguments{
  \item{session}{the session key returned from getSession}
}

\details{  Lists the available folder paths relative to the current folder path for a LabKey session
}

\value{ A character array containing the available folder paths, relative to the project root. These values can be set on a session using curFolder<-
}
\references{https://www.labkey.org/Documentation/wiki-page.view?name=projects}
\author{Peter Hussey}
\seealso{\code{\link{getSession}}, \code{\link{lsProjects}}, \code{\link{lsSchemas}}}
\examples{
\donttest{

##get a list if projects and folders
# library(Rlabkey)

lks<- getSession("https://www.labkey.org", "/home")

#returns values "/home" , "/home/_menus" , ...
lsFolders(lks)

}
}
\keyword{file}
