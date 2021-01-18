\name{labkey.domain.createDesign}
\alias{labkey.domain.createDesign}
\title{Helper function to create a domain design data structure}
\description{
Create a domain design data structure which can then be used by \code{labkey.domain.create} or \code{labkey.domain.save}
}
\usage{
labkey.domain.createDesign(name, description = NULL, fields, indices = NULL)
}
\arguments{
  \item{name}{a string specifying the name of the domain}
  \item{description}{(optional) a string specifying domain description }
  \item{fields}{a list containing the fields of the domain design, this should be in the same format as returned by \code{labkey.inferFields}.}
  \item{indices}{(optional) a list of indices definitions to be used for this domain design on creation}
}
\details{
This is a function which can be used to create a domain design data structure. Domain designs are used both when creating a new domain or updating an existing domain.
}
\value{
A list containing elements describing the domain design. Any of the APIs which take a domain design parameter can accept this data structure.
}

\author{Karl Lum}
\seealso{
\code{\link{labkey.domain.get}},
\code{\link{labkey.domain.inferFields}},
\code{\link{labkey.domain.createIndices}},
\code{\link{labkey.domain.create}},
\code{\link{labkey.domain.save}},
\code{\link{labkey.domain.drop}},
\code{\link{labkey.domain.createConditionalFormat}},
\code{\link{labkey.domain.createConditionalFormatQueryFilter}},
\code{\link{labkey.domain.FILTER_TYPES}}
}
\examples{
\donttest{

## create a data frame and infer it's fields, then create a domain design from it
library(Rlabkey)

df <- data.frame(ptid=c(1:3), age = c(10,20,30), sex = c("f", "m", "f"))
fields <- labkey.domain.inferFields(baseUrl="http://labkey/", folderPath="home", df=df)
indices = labkey.domain.createIndices(list("ptid", "age"), TRUE)
indices = labkey.domain.createIndices(list("age"), FALSE, indices)
dd <- labkey.domain.createDesign(name="test list", fields=fields, indices=indices)

}
}
\keyword{IO}
