\name{writeFinalInfo}
\alias{writeFinalInfo}
\title{
Pretty output, and management of output files
}
\description{
Final code of the R script written by the Migraine software (\url{http://kimura.univ-montp2.fr/~rousset/Migraine.htm}; see main documentation page for the package, for the context). This prints some information, close output files, and beeps to warn that a possibly long computation is finished.
}
\usage{
writeFinalInfo(cleanResu = "")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cleanResu}{
  A connection, or a character string naming a file for some nicely formated output. If \code{""} (the default), print to the standard output connection.
}
}
\value{
returns NULL invisibly.
}
