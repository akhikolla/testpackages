\name{labkey.saveBatch}
\alias{labkey.saveBatch}
\title{Save an assay batch object to a labkey database}
\description{
Save an assay batch object to a labkey database
}
\usage{
labkey.saveBatch(baseUrl, folderPath, assayName, resultDataFrame,
    batchPropertyList=NULL, runPropertyList=NULL)
}
\arguments{
  \item{baseUrl}{a string specifying the \code{baseUrl}for the labkey server}
  \item{folderPath}{a string specifying the \code{folderPath} }
  \item{assayName}{a string specifying the  name of the assay instance}
  \item{resultDataFrame}{a data frame containing rows of data to be inserted}
  \item{batchPropertyList}{a list of batch Properties }
  \item{runPropertyList}{a list of run Properties }
  
}
\details{
This function has been deprecated and will be removed in a future release, please use labkey.experiment.saveBatch instead as
it supports the newer options for saving batch objects.

To save an R data.frame an assay results sets, you must create a named assay using the "General" assay provider.  Detailed 
instructions are available in the Rlabkey Users Guide, accessible by entering RlabkeyUsersGuide() at the R command prompt.
Note that saveBatch currently supports only a single run with one result set per batch.
}

\value{
 Returns the object representation of the Assay batch.
}
\references{https://www.labkey.org/Documentation/wiki-page.view?name=createDatasetViaAssay}

\author{Peter Hussey}
\seealso{
\code{\link{labkey.selectRows}}, \code{\link{labkey.executeSql}}, \code{\link{makeFilter}}, 
\code{\link{labkey.updateRows}}, \cr
\code{\link{labkey.deleteRows}},
\code{\link{labkey.experiment.saveBatch}}
}
\examples{
\donttest{

## Very simple example of an analysis flow:  query some data, calculate
## some stats, then save the calculations as an assay result set in
## LabKey Server
##  Note this example expects to find an assay named "SimpleMeans" in
##  the apisamples project
# library(Rlabkey)

simpledf <- labkey.selectRows(
	baseUrl="http://localhost:8080/labkey",
	folderPath="/apisamples",
	schemaName="lists", 
	queryName="AllTypes")

## some dummy calculations to produce and example analysis result
testtable <- simpledf[,3:4]
colnames(testtable) <- c("IntFld", "DoubleFld")
row <- c(list("Measure"="colMeans"), colMeans(testtable, na.rm=TRUE))
results <- data.frame(row, row.names=NULL, stringsAsFactors=FALSE)
row <- c(list("Measure"="colSums"), colSums(testtable, na.rm=TRUE))
results <- rbind(results, as.vector(row))

bprops <- list(LabNotes="this is a simple demo")
bpl <- list(name=paste("Batch ", as.character(date())),properties=bprops) 
rpl <- list(name=paste("Assay Run ", as.character(date())))

assayInfo<- labkey.saveBatch(
	baseUrl="http://localhost:8080/labkey",
	folderPath="/apisamples",
	"SimpleMeans", 
	results, 
	batchPropertyList=bpl,
	runPropertyList=rpl
)

}
}
\keyword{IO}
