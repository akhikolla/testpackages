\name{saveOldFile}
\alias{saveOldFile}
\title{
Save a copy of an existing file.
}
\description{
This checks if a file of given name already exists in the current directory, and if so saves a copy of it under an automatically generated name (see below).  
}
\usage{
saveOldFile(filename)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filename}{
Name of file to be saved.
}
}
\details{
This function copies the file named \dQuote{\emph{first names}}{\emph{.ext}} under a name created by inserting a string of the form \code{.old_}\eqn{n} between \dQuote{\emph{first names}} and \emph{.ext}, where \eqn{n} is one more than the highest value for any file, matching the first names and extension, already in the current directory, and 0 if no file matches. For example, if \code{filename} is \code{my.beautiful.pdf}, it is copied as \code{my.beautiful.old_0.pdf} if no \code{my.beautiful.old_}\eqn{n}\code{.pdf} file exists, and is is copied as \code{my.beautiful.old_4.pdf} if \code{my.beautiful.old_3.pdf} (and any lower \eqn{n}) file exists.
}
\value{
Returns \code{""} if no file with given name was present on disk, FALSE if it failed to copy an existing old file, 
and the name of the copy if it successfully copied such a file.
}
\examples{
\dontrun{
saveOldFile("same.story")
}
}