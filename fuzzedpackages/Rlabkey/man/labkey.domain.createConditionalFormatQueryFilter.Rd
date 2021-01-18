\name{labkey.domain.createConditionalFormatQueryFilter}
\alias{labkey.domain.createConditionalFormatQueryFilter}
\title{Create a conditional format query filter}
\description{
Create a conditional format query filter string.
}
\usage{
labkey.domain.createConditionalFormatQueryFilter(filterType, value,
    additionalFilter=NULL, additionalValue=NULL)
}
\arguments{
  \item{filterType}{a string specifying a permitted relational operator}
  \item{value}{a string specifying a comparand}
  \item{additionalFilter}{a string specifying a second relational operator}
  \item{additionalValue}{a string specifying a second comparand}

}
\details{
This function can be used to as a convenience wrapper to construct a query filter string for conditional formats.
Two relational expressions may be formed, one with the first two parameters (for instance, parameter values '50' and 'eq' for value and filter
respectively would create a condition of 'equals 50') and the second with the remaining two optional parameters. If both
conditions are created, they are conjunct with a logical AND, and a value would have to pass both conditions to clear
the filter.
This function can be used in conjunction with \code{labkey.domain.FILTER_TYPES} for easy access to the set of
permitted relational operators.
}
\value{
The string specifying a query filter in LabKey filter URL format.
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
\code{\link{labkey.domain.FILTER_TYPES}}
}
\examples{
\donttest{

library(Rlabkey)

qf <- labkey.domain.FILTER_TYPES

# Filters for values equal to 750
qf1 <- labkey.domain.createConditionalFormatQueryFilter(qf$EQUAL, 750)
# Filters for values greater than 500, but less than 1000
qf2 <- labkey.domain.createConditionalFormatQueryFilter(qf$GREATER_THAN, 500, qf$LESS_THAN, 1000)

}
}
\keyword{IO}
