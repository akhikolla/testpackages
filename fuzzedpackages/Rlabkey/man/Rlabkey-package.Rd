\name{Rlabkey-package}
\alias{Rlabkey-package}
\alias{Rlabkey}
\docType{package}
\title{
Exchange data between LabKey Server and R
}
\description{
This package allows the transfer of data between a LabKey Server and an R session. Data can be retrieved from 
LabKey into a data frame in R by specifying the query schema information (\code{labkey.selectRows} and \code{getRows}) 
or by using sql commands (\code{labkey.executeSql}). From an R session, existing data can be updated
 (\code{labkey.updateRows}), new data can be inserted (\code{labkey.insertRows} and \code{labkey.importRows}) or
data can be deleted from the LabKey database (\code{labkey.deleteRows}). Interactive R users can discover available data via
schema objects (\code{labkey.getSchema}).

}
\details{
\tabular{ll}{
Package: \tab Rlabkey\cr
Type: \tab Package\cr
Version: \tab 2.5.5\cr
Date: \tab 2020-12-11\cr
License: \tab Apache License 2.0\cr
LazyLoad: \tab yes\cr
}
The user must have appropriate authorization on the LabKey Server in order to access or modify data using these functions.
All access to secure content on LabKey Server requires authentication via an api key (see \code{labkey.setDefaults} for
more details) or a properly configured netrc file that includes the user's login information.

The netrc file is a standard mechanism for conveying configuration and autologin information to the File Transfer Protocol
client (ftp) and other programs such as CURL. On a Linux or Mac system this file should be named .netrc (dot netrc) and on
Windows it should be named _netrc (underscore netrc). The file should be located in the user's home directory and the
permissions on the file should be unreadable for everybody except the owner.

To create the _netrc on a Windows machine, first create an environment variable called 'HOME' set to your home directory
(e.g., c:/Users/<User-Name> on recent versions of Windows) or any directory you want to use. In that directory, create a
text file named _netrc (note that it's underscore netrc, not dot netrc like it is on Linux/Mac).

The following three lines must be included in the .netrc or _netrc file either separated by white space
(spaces, tabs, or newlines) or commas.

machine <remote-machine-name>\cr
login <user-email>\cr
password <user-password>


One example would be:\cr
machine localhost\cr
login peter@labkey.com\cr
password mypassword\cr

Another example would be:\cr
machine atlas.scharp.org login vobencha@fhcrc.org password mypassword\cr
\cr
Multiple such blocks can exist in one file.
}

\author{Valerie Obenchain
}
\references{http://www.omegahat.net/RCurl/,\cr
https://www.labkey.org/project/home/begin.view}
\keyword{ package }
\seealso{
\code{\link{labkey.selectRows}}, \code{\link{labkey.executeSql}}, \code{\link{makeFilter}}, 
\code{\link{labkey.insertRows}}, \code{\link{labkey.importRows}}, \cr
\code{\link{labkey.updateRows}}, \code{\link{labkey.deleteRows}}
\cr
The Rlabkey Users Guide is available by typing RlabkeyUsersGuide().
}

