\name{labkey.domain.createIndices}
\alias{labkey.domain.createIndices}
\title{Helper function to create a domain design indices list}
\description{
Create a list of indices definitions which can then be used by \code{labkey.domain.createDesign}
}
\usage{
labkey.domain.createIndices(colNames, asUnique, existingIndices = NULL)
}
\arguments{
  \item{colNames}{a list of string column names for the index}
  \item{asUnique}{a logical TRUE or FALSE value for if a UNIQUE index should be used}
  \item{existingIndices}{a list of previously created indices definitions to append to}
}
\details{
This helper function can be used to construct the list of indices definitions for a domain design structure.
Each call to this function takes in the column names from the domain to use in the index and a parameter
indicating if this should be a UNIQUE index. A third parameter can be used to build up more then one indices
definitions.
}
\value{
The data frame containing the list of indices definitions, concatenated with the existingIndices object if provided.
}

\author{Cory Nathe}
\seealso{
\code{\link{labkey.domain.get}},
\code{\link{labkey.domain.create}},
\code{\link{labkey.domain.createDesign}},
\code{\link{labkey.domain.inferFields}},
\code{\link{labkey.domain.save}},
\code{\link{labkey.domain.drop}}
}
\examples{
\donttest{

## create a list of indices definitions to use for a domain design
library(Rlabkey)

indices = labkey.domain.createIndices(list("intKey", "customInt"), TRUE)
indices = labkey.domain.createIndices(list("customInt"), FALSE, indices)

}
}
\keyword{IO}
