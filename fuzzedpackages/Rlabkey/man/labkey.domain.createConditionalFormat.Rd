\name{labkey.domain.createConditionalFormat}
\alias{labkey.domain.createConditionalFormat}
\title{Create a conditional format data frame}
\description{
Create a conditional format data frame.
}
\usage{
labkey.domain.createConditionalFormat(queryFilter, bold=FALSE, italic=FALSE,
    strikeThrough=FALSE, textColor="", backgroundColor="")
}
\arguments{
  \item{queryFilter}{a string specifying what logical filter should be applied}
  \item{bold}{a boolean for if the text display should be formatted in bold}
  \item{italic}{a boolean for if the text display should be formatted in italic}
  \item{strikeThrough}{a boolean for if the text display should be formatted with a strikethrough}
  \item{textColor}{a string specifying the hex code of the text color for display}
  \item{backgroundColor}{a string specifying the hex code of the text background color for display}

}
\details{
This function can be used to construct a conditional format data frame intended for use within a domain design's
conditionalFormats component while creating or updating a domain. The queryFilter parameter can be used in conjunction
with \code{labkey.domain.createConditionalFormatQueryFilter} for convenient construction of a query filter string. Multiple conditional
formats can be applied to one field, where each format specified constitutes a new row of the field's conditionalFormats
data frame.
If text formatting options are not specified, the default is to display the value as black text on a white background.
}
\value{
The data frame containing values describing a conditional format.
}

\author{Rosaline Pyktel}
\seealso{
\code{\link{labkey.domain.get}},
\code{\link{labkey.domain.create}},
\code{\link{labkey.domain.createDesign}},
\code{\link{labkey.domain.inferFields}},
\code{\link{labkey.domain.save}},
\code{\link{labkey.domain.drop}},
\code{\link{labkey.domain.createConditionalFormatQueryFilter}},
\code{\link{labkey.domain.FILTER_TYPES}}
}
\examples{
\donttest{

library(Rlabkey)

domain <- labkey.domain.get(baseUrl="http://labkey/", folderPath="home",
    schemaName="lists", queryName="test list")

## update the third field to use two conditional formats
qf <- labkey.domain.FILTER_TYPES
cf1 = labkey.domain.createConditionalFormat(labkey.domain.createConditionalFormatQueryFilter(qf$GT,
      100), bold=TRUE, text_color="D33115", background_color="333333")
cf2 = labkey.domain.createConditionalFormat(labkey.domain.createConditionalFormatQueryFilter(
      qf$LESS_THAN, 400), italic=TRUE, text_color="68BC00")
domain$fields$conditionalFormats[[3]] = rbind(cf1,cf2)

labkey.domain.save(baseUrl="http://labkey/", folderPath="home",
    schemaName="lists", queryName="test list", domainDesign=domain)

}
}
\keyword{IO}
