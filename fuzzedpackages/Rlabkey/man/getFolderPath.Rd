\name{getFolderPath}
\alias{getFolderPath}
\title{Returns the folder path associated with a session}
\description{ Returns the current folder path for a LabKey session }

\usage{
getFolderPath(session) 
}
\arguments{
  \item{session}{the session key returned from getSession}
}

\details{  Returns a string containing the current folder path for the passed in LabKey session
}

\value{ A character array containing the folder path, relative to the root. 
}
\references{https://www.labkey.org/Documentation/wiki-page.view?name=projects}
\author{Peter Hussey}
\seealso{\code{\link{getSession}}  \code{\link{lsFolders}}}
\examples{
\donttest{

# library(Rlabkey)

lks<- getSession("https://www.labkey.org", "/home")
getFolderPath(lks)  #returns "/home"

}
}
\keyword{file}
