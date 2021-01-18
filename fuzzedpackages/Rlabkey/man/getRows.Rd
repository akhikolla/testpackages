\name{getRows}
\alias{getRows}
\title{ Retrieve data from LabKey Server}
\description{ Retrive rows from a LabKey Server given a session and query object }

\usage{
getRows(session, query, maxRows=NULL, colNameOpt='fieldname', ...) 
}
\arguments{
  \item{session}{the session key returned from getSession}
  \item{query}{an object representing a query on LabKey Server, a child object of the object returned by \code{getSchema()}}
  \item{maxRows}{(optional) an integer specifying how many rows of data to return. If no value is specified, all rows are returned.}
 \item{colNameOpt}{(optional) controls the name source for the columns of the output dataframe, with valid values of 'caption', 'fieldname', and 'rname'}
 \item{...}{Any of the remaining options to \code{link{labkey.selectRows}}}
}

\details{  This function works as a shortcut wrapper to \code{\link{labkey.selectRows}}.  All of the arguments are the same as documented in
labkey.selectRows. 

See \code{\link{labkey.selectRows}} for a discussion of the valid options and defaults for \code{colNameOpt}.  Note in particular 
that with getRows the default is 'fieldname' instead of 'caption'.  
}

\value{ A data frame containing the query results corresponding to the default view of the specified query.
}
\author{Peter Hussey}
\seealso{\code{\link{getSession}}, \code{\link{getSchema}}, \code{\link{getLookups}}, \code{\link{saveResults}}
\code{\link{labkey.selectRows}}
}
\examples{
\donttest{

## simple example of getting data using schema objects
# library(Rlabkey)

s<-getSession(baseUrl="http://localhost:8080/labkey",
    folderPath="/apisamples")
s  # shows schemas

scobj <- getSchema(s, "lists")
scobj   # shows available queries

scobj$AllTypes  ## this is the query object

getRows(s, scobj$AllTypes)
	
}
}
\keyword{file}
