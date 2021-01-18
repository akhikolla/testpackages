\name{saveResults}
\alias{saveResults}
\title{ Returns an object representing a LabKey schema }
\description{
A wrapper function to labkey.saveBatch which uses a session object and provides defaults for the Batch/Run names.
 }

\usage{
saveResults(session, assayName, resultDataFrame, 
  batchPropertyList= list(name=paste("Batch ", as.character(date()))),
  runPropertyList= list(name=paste("Assay Run ", as.character(date()))))
}
\arguments{
  \item{session}{the session key returned from getSession}
  \item{assayName}{a string specifying the  name of the assay instance}
  \item{resultDataFrame}{a data frame containing rows of data to be inserted}
  \item{batchPropertyList}{a list of batch Properties }
  \item{runPropertyList}{a list of run Properties }

}
\details{  saveResults is a wrapper function to labkey.saveBatch with two changes: First, it uses a session object in place
of the separate baseUrl and folderPath arguments. Second, it provides defaults for generating Batch and Run names based on a current
timestamp.

To see the save result on LabKey server, click on the "SimpleMeans" assay in the Assay List web part.
}

\value{an object representing the assay.
}
\references{
https://www.labkey.org/project/home/begin.view}
\author{Peter Hussey}
\seealso{\code{\link{getSession}}, \code{\link{getSchema}}, \code{\link{getLookups}} \code{\link{getRows}}}
\examples{
\donttest{

## Very simple example of an analysis flow:  query some data,
## calculate some stats, then save the calculations as an assay
## result set in LabKey Server
# library(Rlabkey)

s<- getSession(baseUrl="http://localhost:8080/labkey",
    folderPath="/apisamples")
scobj <- getSchema(s, "lists")
simpledf <- getRows(s, scobj$AllTypes)

## some dummy calculations to produce and example analysis result
testtable <- simpledf[,3:4]
colnames(testtable) <- c("IntFld", "DoubleFld")
row <- c(list("Measure"="colMeans"), colMeans(testtable, na.rm=TRUE))
results <- data.frame(row, row.names=NULL, stringsAsFactors=FALSE)
row <- c(list("Measure"="colSums"), colSums(testtable, na.rm=TRUE))
results <- rbind(results, as.vector(row))

bprops <- list(LabNotes="this is a simple demo")
bpl<- list(name=paste("Batch ", as.character(date())),properties=bprops) 
rpl<- list(name=paste("Assay Run ", as.character(date())))

assayInfo<- saveResults(s, "SimpleMeans", results,
    batchPropertyList=bpl, runPropertyList=rpl)

}
}
\keyword{file}
