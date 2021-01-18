\name{descstat}
\alias{descstat}

\title{Computes Basic Descriptive Statistics.}

\description{
This function computes descriptive statistics which are most useful for Dynamic Model Averaging. 

It is a wrapper of \code{\link[psych]{describe}}.

If the argument is not a \code{\link[base]{matrix}}, the function tries to convert the object into a \code{\link[base]{matrix}}. For example, it works smoothly for \code{\link[xts]{xts}} objects.
}

\details{
See \code{\link[psych]{describe}}.
}

\usage{
descstat(data)
}

\arguments{
\item{data}{\code{\link[base]{matrix}}, observations are put in rows, and variables are grouped by columns}
}

\value{\code{\link[base]{matrix}}}

\examples{
descstat(crudeoil)
}
