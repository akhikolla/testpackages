\name{getSession}
\alias{getSession}
\title{Creates and returns a LabKey Server session}
\description{
The session object holds server address and folder context for a user working with LabKey Server.  The session-based model 
supports more efficient interactive use of LabKey Server from R.
}
\usage{
getSession(baseUrl, folderPath="/home",
    curlOptions=NULL, lkOptions=NULL)
}
\arguments{
  \item{baseUrl}{a string specifying the address of the LabKey Server, including the context root}
  \item{folderPath}{a string specifying the hierarchy of folders to the current folder (container) for the operation, starting with the project folder }
  \item{curlOptions}{(optional) a list of curlOptions to be set on connections to the LabKey Server, see details}
  \item{lkOptions}{(optional) a list of settings for default behavior on naming of objects, see details}
}
\details{
Creates a session key that is passed to all the other session-based functions.  Associated with the key are a baseUrl and a folderPath which
determine the security context.  

\code{curlOptions}

The curlOptions parameter gives a mechanism to pass control options down to the RCurl library used by Rlabkey.  
This can be very useful for debugging problems or setting proxy server properties.  See example for debugging.

\code{lkOptions}  

The lkOptions parameter gives a mechanism to change default behavior in the objects returned by Rlabkey.  
Currently the only available options are \code{colNameOpt}, which affects the names of columns in the data frames returned
by getRows(), and \code{maxRows}, which sets a default value for this parameter when calling getRows() 
}

\value{
getSession returns a session object that represents a specific user within a specific project folder within the LabKey Server identified by the baseUrl.  
The combination of user, server and project/folder determines the security context of the client program.  See the Rlabkey Users Guide for more
discussion of how the user identity is established.
}
\author{Peter Hussey}
\seealso{\code{\link{getRows}}, \code{\link{getSchema}}, \code{\link{getLookups}} \code{\link{saveResults}}\cr
The Rlabkey Users Guide is available by typing RlabkeyUsersGuide().
}
\examples{
\donttest{

# library(Rlabkey)

s <- getSession("https://www.labkey.org", "/home")
s   #shows schemas

## using the curlOptions for generating debug tracesof network traffic
d<- debugGatherer()
copt <- curlOptions(debugfunction=d$update, verbose=TRUE,
    cookiefile='/cooks.txt')
sdbg<- getSession(baseUrl="http://localhost:8080/labkey",
    folderPath="/apisamples", curlOptions=copt)
getRows(sdbg, scobj$AllTypes)
strwrap(d$value(), 100)

}
}
\keyword{file}
