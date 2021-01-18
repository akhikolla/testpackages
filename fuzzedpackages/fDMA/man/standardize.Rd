\name{standardize}
\alias{standardize}

\title{Standardizes a Numeric Matrix by Columns.}

\description{
Sometimes it is desirable to have all variables to have mean 0 and standard deviation 1. This function rescales the values in such a way.

If the argument is not a \code{\link[base]{matrix}}, the function tries to convert the object into a \code{\link[base]{matrix}}. For example, it works smoothly for \code{\link[xts]{xts}} objects.
}

\usage{
standardize(data)
}

\arguments{
\item{data}{\code{\link[base]{matrix}}, observations are put in rows, and variables are grouped by columns}
}

\value{\code{\link[base]{matrix}}}

\examples{
standardize(crudeoil)
}

\seealso{\code{\link{normalize}}}
