\name{labkey.getRequestOptions}
\alias{labkey.getRequestOptions}
\title{Helper function to get the HTTP request options for a specific method type.}
\description{
The internal functions for labkey.get() and labkey.post() use this labkey.getRequestOptions()
helper to build up the HTTP request options for things like CSRF, CURL options, and authentication
properties. This function is also exposed for general use if you would like to make your own
HTTP request but need to use those request options as set in your session context.
}
\usage{
labkey.getRequestOptions(method = 'GET', encoding = NULL)
}
\arguments{
  \item{method}{a string specifying the HTTP method for the request options you want to get}
  \item{encoding}{a string specifying the type of encoding to add to the header properties, defaults to UTF-8 when NULL}
}

\author{Cory Nathe}
\examples{
\donttest{

library(Rlabkey)
labkey.getRequestOptions()

}
}
\keyword{IO}
