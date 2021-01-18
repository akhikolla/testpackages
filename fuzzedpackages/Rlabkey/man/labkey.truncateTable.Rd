\name{labkey.truncateTable}
\alias{labkey.truncateTable}
\title{Delete all rows from a table}
\description{
Delete all rows from the specified table.
}
\usage{
labkey.truncateTable(baseUrl = NULL, folderPath, schemaName, queryName)
}
\arguments{
  \item{baseUrl}{a string specifying the \code{baseUrl} for the labkey server}
  \item{folderPath}{a string specifying the \code{folderPath} }
  \item{schemaName}{a string specifying the name of the schema of the domain }
  \item{queryName}{a string specifying the query name }
}
\details{
Deletes all rows in the table in a single transaction and will also log a single audit event for the action. Not all tables support truncation, if
a particular table doesn't support the action, an error will be returned. The current list of tables supporting truncation include : lists, datasets,
issues, sample sets, data classes.
}
\value{
Returns the count of the number of rows deleted.
}

\author{Karl Lum}
\seealso{
\code{\link{labkey.deleteRows}}
}
\examples{
\donttest{

## create a data frame and infer it's fields
library(Rlabkey)

labkey.truncateTable(baseUrl="http://labkey/", folderPath="home",
    schemaName="lists", queryName="people")

}
}
\keyword{IO}
