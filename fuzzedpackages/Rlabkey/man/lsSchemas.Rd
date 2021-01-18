\name{lsSchemas}
\alias{lsSchemas}
\title{List the available schemas }
\description{ Lists the available schemas given the current folder path for a LabKey session }

\usage{
lsSchemas(session) 
}
\arguments{
  \item{session}{the session key returned from getSession}
}

\details{ Lists the available schemas given the current folder path for a LabKey session
}

\value{ A character array containing the available schema names
}
\author{Peter Hussey}
\seealso{\code{\link{getSession}}, {\code{\link{lsFolders}}, \code{\link{lsProjects}}}
}
\examples{
\donttest{

##  get a list of schemas available in the current session context
# library(Rlabkey)

lks<- getSession(baseUrl="http://localhost:8080/labkey",
    folderPath="/apisamples")

#returns several schema names, e.g. "lists", "core", "MS1", etc.
lsSchemas(lks)

}
}
\keyword{file}
