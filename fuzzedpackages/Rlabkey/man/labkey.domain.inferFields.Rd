\name{labkey.domain.inferFields}
\alias{labkey.domain.inferFields}
\title{Infer field metadata from a data frame}
\description{
Generate field information from the specified data frame. The resulting list can be used to create or edit a domain using the labkey.domain.create or labkey.domain.save APIs.
}
\usage{
labkey.domain.inferFields(baseUrl = NULL, folderPath, df)
}
\arguments{
  \item{baseUrl}{a string specifying the \code{baseUrl} for the labkey server}
  \item{folderPath}{a string specifying the \code{folderPath} }
  \item{df}{a data frame specifying fields to infer. The data frame must have column names as well as row data to infer the type of the field from.}
}
\details{
Field information can be generated from a data frame by introspecting the data associated with it along with other properties about that column. The data frame
is posted to the server endpoint where the data is analyzed and returned as a list of fields each with it's associated list of properties and values. This list
can be edited and/or used to create a domain on the server.
}
\value{
The inferred metadata will be returned as a list with an element called : "fields" which contains the list of fields inferred from the data frame.
Each field will contain the list of attributes and values for that field definition.
}

\author{Karl Lum}
\seealso{
\code{\link{labkey.domain.get}},
\code{\link{labkey.domain.create}},
\code{\link{labkey.domain.createDesign}},
\code{\link{labkey.domain.createIndices}},
\code{\link{labkey.domain.save}},
\code{\link{labkey.domain.drop}}
}
\examples{
\donttest{

## create a data frame and infer it's fields
library(Rlabkey)

df <- data.frame(ptid=c(1:3), age = c(10,20,30), sex = c("f", "m", "f"))
fields <- labkey.domain.inferFields(baseUrl="http://labkey/", folderPath="home", df=df)

}
}
\keyword{IO}
