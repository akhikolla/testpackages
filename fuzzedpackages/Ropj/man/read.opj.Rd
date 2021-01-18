% vi:et:ts=2
\name{read.opj}
\alias{read.opj}
\docType{package}
\title{Parse Origin(R) project file into a list of objects}
\description{
  This function parses an OPJ file into a list of objects it consists
  of. Items understood by \code{read.opj} include spreadsheets, matrices,
  and notes.
}
\usage{
  read.opj(file, encoding = 'latin1', tree = FALSE, ...)
}
\arguments{
  \item{file}{
    Path to the OPJ file to parse.  Only file paths are supported,
    not R connections.
  }
  \item{encoding}{
    Encoding of the strings inside the file being opened. This should
    correspond to the ANSI code page of Windows installation used to
    produce the file. The default of \code{'latin1'} is usually safe to
    use. See \code{\link{iconvlist}()} for the full list of encodings
    supported by your R installation.
  }
  \item{tree}{
    Control the structure of the returned list. When \code{FALSE}
    (default), returned list is flat and its unique names correspond
    to the short names of the objects, while the \code{\link{comment}}
    attribute contains their long names and comments. When \code{TRUE},
    the list itself becomes a recursive data structure containing the
    tree of the objects, making it possible to access objects by their
    paths in that tree.
  }
  \item{\dots}{
    The rest of the arguments is passed to \code{\link{as.data.frame}}
    when converting spreadsheets from lists columns, making it possible
    to set \code{stringsAsFactors} and other parameters as needed.
  }
}
\value{
  A named \code{\link{list}} containing objects stored in the file.

  \itemize{
    \item Spreadsheets are presented as \code{\link{data.frame}}s,
    with additional \code{\link{attributes}}:
    \itemize{
      \item \code{\link{comment}} contains the fields Long name, Units,
      and Comment, joined by \code{\\r\\n}, if present.  Due to a
      possible bug in \code{liborigin}, these values are sometimes
      followed by \code{@} and some text that wasn't present in
      the original document. (In versions prior to v0.2-2 it was called
      \code{comments}, which should be still supported until v1.0.)
      \item \code{commands} contains the formula that was used to
      create the values of the column (e.g. \code{col(A) * 2 + 1}).
    }
    \item Multi-sheet spreadsheets are stored as named lists of
    \code{\link{data.frame}}s described above.
    \item Matrices are presented as \code{\link{list}}s
    of \code{\link{matrix}} objects containing numeric
    data. \code{\link{dimnames}} are also assigned.  The list also
    has attributes:
    \itemize{
      \item \code{commands} contains the formula that was used to
      compute the values in the matrix.
    }
    \item Notes are stored as plain strings.
  }

  When \code{tree = FALSE}, the list is flat, its names are short names
  of the objects, and the \code{\link{comment}} attribute of the list
  contains the long names of the objects
  stored in the file.

  When \code{tree = TRUE}, the list names are long names (if available)
  and the list itself represents the folder structure of the project.
}
\note{
  While Origin(R) and its scripting language seem to rely on the
  \emph{short} names of the objects being unique across a project,
  neither \emph{long} names nor \emph{folder} names are guaranteed to
  avoid collisions. Tree-like lists returned by \code{read.opj(..., tree =
  TRUE)} might be easier to navigate interactively but present problems
  if someone gives one or many folders or objects the same long name.
}
\examples{
  x <- read.opj(system.file('test.opj', package = 'Ropj'))
  head(x$Book2, 7)
  comment(x$Book2)
  attr(x$Book2, 'commands')
  with(x$Book1, head(Sheet2 - Sheet1))
  x$MBook1$MSheet1[1:4,1:4]
  x$Note1

  \dontrun{
    # encoding names aren't guaranteed to be supported across all platforms
    x <- read.opj(system.file('test.opj', package = 'Ropj'), 'WINDOWS-1251')
    x$cyrillic
  }

  str(read.opj(system.file('tree.opj', package = 'Ropj'), tree = TRUE))
}
\keyword{IO}
\keyword{file}
