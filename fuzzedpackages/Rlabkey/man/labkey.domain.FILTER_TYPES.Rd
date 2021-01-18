\name{labkey.domain.FILTER_TYPES}
\alias{labkey.domain.FILTER_TYPES}
\title{Provide comparator access}
\description{
A list specifying permitted validator comparators.
}
\usage{
labkey.domain.FILTER_TYPES
}
\details{
This constant contains a list specifying the set of permitted validator operators, using names to map
conventional terms to the expressions used by LabKey filter URL formats. The values are intended to be used in
conjunction with conditional formats or property validators.
}
\value{
A named list of strings.
}

\author{Rosaline Pyktel}
\seealso{
\code{\link{labkey.domain.get}},
\code{\link{labkey.domain.create}},
\code{\link{labkey.domain.createDesign}},
\code{\link{labkey.domain.inferFields}},
\code{\link{labkey.domain.save}},
\code{\link{labkey.domain.drop}},
\code{\link{labkey.domain.createConditionalFormat}},
\code{\link{labkey.domain.createConditionalFormatQueryFilter}}
}
\examples{
\donttest{
library(Rlabkey)

qf <- labkey.domain.FILTER_TYPES

# Example of available comparators
comparator1 <- qf$EQUAL
comparator2 <- qf$GREATER_THAN
comparator3 <- qf$DATE_LESS_THAN_OR_EQUAL
comparator4 <- qf$STARTS_WITH
comparator5 <- qf$CONTAINS_ONE_OF

}
}
\keyword{IO}
