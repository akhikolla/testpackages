\name{icons}
\alias{icons}
\alias{icons_matrix}
\alias{oneill}
\docType{data}
\title{Dataset on climate change due to O'Neill}
\description{
  Object \code{icons_matrix} is 
  a matrix of nine rows and six columns, one column for each of six
  icons relevant to climate change.  The matrix entries show the number
  of respondents who indicated which icon they found most concerning.
  The nine rows show different classes of respondents who were exposed
  to different subsets (of size four) of the six icons.

  The columns correspond to the different stimulus icons used, detailed
  below.

  Object \code{icons} is the corresponding likelihood function, which
  can be created with \code{saffy(icons_matrix)}.
}
\usage{data(oneill)}
\details{
  The six icons were used in this study were:
  \describe{
    \item{PB}{polar bears, which face extinction through loss of ice
    floe hunting grounds}
    \item{NB}{The Norfolk Broads, which flood due to intense rainfall
    events}
    \item{LF}{London flooding, as a result of sea level rise}
    \item{THC}{The Thermo-haline circulation, which may slow or stop as
      a result of anthropogenic modification of the hydrological cycle}
    \item{OA}{Oceanic acidification as a result of anthropogenic emissions
    of carbon dioxide}
    \item{WAIS}{The West Antarctic Ice Sheet, which is calving into the sea
    as a result of climate change}
  }
}
\source{
  Data kindly supplied by Saffron O'Neill of the University of East Anglia
}
\references{
  \itemize{
    \item
    S. O'Neill 2007. \emph{An Iconic Approach to Communicating
      Climate Change}, University of East Anglia,  School of Environmental
    Science (in prep)
    \item I. Lorenzoni and N. Pidgeon 2005. \emph{Defining Dangers of
      Climate Change and Individual Behaviour: Closing the Gap}. 
    In  \emph{Avoiding Dangerous Climate Change} (conference
    proceedings), UK Met Office, Exeter, 1-3 February
    \item R. K. S. Hankin 20101. \dQuote{A generalization of the Dirichlet
      distribution}.   \emph{Journal of Statistical software}, 33:11
  }
}
\examples{
data(oneill)
maxp(icons)

icons == saffy(icons_matrix)   # should be TRUE
}
\keyword{datasets}
