\name{makeFilter}
\alias{makeFilter}
\title{Builds filters to be used in labkey.selectRows and getRows}
\description{
This function takes inputs of column name, filter value and filter operator and 
returns an array of filters to be used in \code{labkey.selectRows} and \code{getRows}.
}
\usage{
makeFilter(...)
}
\arguments{
\item{...}{Arguments in c("colname","operator","value") form, used to create a filter.}
}
\details{
These filters are applied to the data prior to import into R. The user can specify as many 
filters as desired. The format for specifying a filter is a vector of characters including
the column name, operator and value.
\describe{
\item{colname}{a string specifying the name of the column to be filtered} 
\item{operator}{a string specifying what operator should be used in the filter (see options below)}
\item{value}{an integer or string specifying the value the columns should be filtered on}
}

Operator values:\cr
EQUAL\cr
DATE_EQUAL\cr
NOT_EQUAL\cr
DATE_NOT_EQUAL\cr
NOT_EQUAL_OR_MISSING\cr
GREATER_THAN\cr
DATE_GREATER_THAN\cr
LESS_THAN\cr
DATE_LESS_THAN\cr
GREATER_THAN_OR_EQUAL\cr
DATE_GREATER_THAN_OR_EQUAL\cr
LESS_THAN_OR_EQUAL\cr
DATE_LESS_THAN_OR_EQUAL\cr
STARTS_WITH\cr
DOES_NOT_START_WITH\cr
CONTAINS\cr
DOES_NOT_CONTAIN\cr
CONTAINS_ONE_OF\cr
CONTAINS_NONE_OF\cr
IN\cr
NOT_IN\cr
BETWEEN\cr
NOT_BETWEEN\cr
MEMBER_OF\cr
MISSING\cr
NOT_MISSING\cr
MV_INDICATOR\cr
NO_MV_INDICATOR\cr
Q\cr

When using the MISSING, NOT_MISSING, MV_INDICATOR, or NO_MV_INDICATOR operators, an empty string should be supplied as the value.
See example below.
}


\value{
The function returns either a single string or an array of strings to be use in the
\code{colFilter} argument of the \code{labkey.selectRows} function.
}
\references{http://www.omegahat.net/RCurl/, \cr
https://www.labkey.org/project/home/begin.view}
\author{Valerie Obenchain}
\seealso{\code{\link{labkey.selectRows}}}
\examples{

# library(Rlabkey)

## Two filters, ANDed together
makeFilter(c("TextFld","CONTAINS","h"),
    c("BooleanFld","EQUAL","TRUE"))

## Using "in" operator:
makeFilter(c("RowId","IN","2;3;6"))

## Using "missing" operator:
makeFilter(c("IntFld","MISSING",""))

}
\keyword{file}
