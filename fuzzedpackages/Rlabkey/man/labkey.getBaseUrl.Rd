\name{labkey.getBaseUrl}
\alias{labkey.getBaseUrl}
\title{Get the default baseUrl parameter used for all http or https requests}
\description{
Use this function to get "baseUrl" package environment variables to be used for all http or https requests.
}
\usage{
labkey.getBaseUrl(baseUrl=NULL)
}
\arguments{
  \item{baseUrl}{server location including context path, if any. e.g. https://www.labkey.org/}
}

\details{
The function takes an optional baseUrl parameter. When non empty parameter is passed in and if baseUrl
has not been previously set, the function will remember the baseUrl value in package environment variables and
return the formatted baseUrl. Skip baseUrl parameter to get previously set baseUrl.
}
\examples{
\donttest{
## Example of getting previously set baseUrl
library(Rlabkey)
labkey.setDefaults(apiKey="session|abcdef0123456789abcdef0123456789",
    baseUrl="http://labkey/")
labkey.getBaseUrl()

}
}
