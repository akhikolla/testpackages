\name{labkey.domain.drop}
\alias{labkey.domain.drop}
\title{Delete a LabKey domain}
\description{
Delete an existing domain.
}
\usage{
labkey.domain.drop(baseUrl=NULL, folderPath, schemaName, queryName)
}
\arguments{
  \item{baseUrl}{a string specifying the \code{baseUrl}for the labkey server}
  \item{folderPath}{a string specifying the \code{folderPath} }
  \item{schemaName}{a string specifying the name of the schema of the domain }
  \item{queryName}{a string specifying the query name }
}
\details{
This function will delete an existing domain along with any data that may have been uploaded to it.
}

\author{Karl Lum}
\seealso{
\code{\link{labkey.domain.get}},
\code{\link{labkey.domain.inferFields}},
\code{\link{labkey.domain.createDesign}},
\code{\link{labkey.domain.createIndices}},
\code{\link{labkey.domain.save}},
\code{\link{labkey.domain.create}}
}
\examples{
\donttest{

## delete an existing domain
library(Rlabkey)

labkey.domain.drop(baseUrl="http://labkey/", folderPath="home",
    schemaName="lists", queryName="test list")
}
}
\keyword{IO}
