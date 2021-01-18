\name{los.data}
\docType{data}
\alias{los.data}
\title{Length of hospital stay}

\description{ The \code{los.data} data frame has 756 rows, one row for
  each patient, and 7 columns.
}

\usage{data(los.data)}

\format{A data frame with the following columns:
  \describe{
  \item{adm.id}{ admision id of the patient}
  \item{j.01}{ observed time for jump from 0 (initial state) to 1
    (intermediate state)}
  \item{j.02}{ observed time for jump from 0 to 2 (discharge)}
  \item{j.03}{ observed time for jump from 0 to 3 (death)}
  \item{j.12}{ observed time for jump from 1 to 2}
  \item{j.13}{ observed time for jump from 1 to 3}
  \item{cens}{ censoring time (either in initial or intermediate state) }
  }
}

\examples{
data(los.data)
my.data <- prepare.los.data(los.data)
}

\keyword{datasets}
