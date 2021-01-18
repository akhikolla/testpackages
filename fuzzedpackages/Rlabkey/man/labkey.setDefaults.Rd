\name{labkey.setDefaults}
\alias{labkey.setDefaults}
\title{Set the default parameters used for all http or https requests}
\description{
Use this function to set the default baseUrl and authentication parameters as package environment variables to be used
for all http or https requests. You can also use labkey.setDefaults() without any parameters to reset/clear
these settings.}
\usage{
labkey.setDefaults(apiKey="", baseUrl="", email="", password="")
}
\arguments{
  \item{apiKey}{session key from your server}
  \item{baseUrl}{server location including context path, if any. e.g. https://www.labkey.org/}
  \item{email}{user email address}
  \item{password}{user password}
}

\details{
Note: Support for API keys was added in LabKey Server release 16.2; they are not supported in 16.1 or earlier.

An API key can be used to authorize Rlabkey functions that access protected content on LabKey Server. Using an API key
avoids copying and storing credentials on the client machine. Also, all Rlabkey access is tied to the current browser
session, which means the code runs in the same context as the browser (e.g. same user, same authorizations, same
declared terms of use and PHI level, same impersonation state, etc.).\cr

A site administrator must first enable the use of API keys on that LabKey Server. Once enabled, any logged in user can
generate an API key by clicking their display name (upper right) and selecting "API Keys". The API Key page creates
and displays keys that can be copied and pasted into a labkey.setDefaults() statement to tie an Rlabkey session to the
authorization and session information already set in the browser.\cr

If an API key is not provided, you can also use this function for basic authentication via email and password. Note
that both email and password must be set via a labkey.setDefaults() call. If an apiKey is also set, that will be given
preference and the email/password will not be used for authentication. Once authenticated via email/password, you
can make multiple labkey.get or labkey.post API calls using that same connection.\cr
}
\examples{

## Example of setting and clearing an API key and/or email/password.
# library(Rlabkey)

labkey.setDefaults(apiKey="session|abcdef0123456789abcdef0123456789")
labkey.setDefaults(email="testing@localhost.test", password="password")

## Functions invoked at this point share authorization
## and session information with the browser session

labkey.setDefaults() # called without any parameters will reset/clear the environment variables

}
