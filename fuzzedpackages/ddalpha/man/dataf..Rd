\name{dataf.*}
\alias{dataf.*}
\docType{data}
\title{
Functional Data Sets
}
\description{
The functions generate data sets of functional two-dimensional data of two or more classes.
}
\usage{
# dataf.[name]()

# load the data set by name
# data(list = "name")

# load the data set by name to a variable
# getdata("name")
}
\format{
  The functional data as a data structure.
  \describe{
    \item{\code{dataf}}{
      The functional data as a list of objects. Each object is characterized by two coordinates 
      \describe{
        \item{\code{args}}{The arguments vector}
        \item{\code{vals}}{The values vector}
       }
    }
    \item{\code{labels}}{The classes of the objects}
  }
}
\details{
More details about the datasets in the topics:

\code{\link{dataf.geneexp}}

\code{\link{dataf.growth}}

\code{\link{dataf.medflies}}

\code{\link{dataf.population}}

\code{\link{dataf.population2010}}

\code{\link{dataf.tecator}}

\code{\link{dataf.tecator}}

The following datasets provide simulated data:

\code{\link{dataf.sim.1.CFF07}}

\code{\link{dataf.sim.2.CFF07}}

}
\seealso{
\code{\link{plot.functional}} for building plots of functional data
}
\examples{
## load the Growth dataset
dataf = dataf.growth()

## view the classes
unique(dataf$labels)

## access the 5th point of the 2nd object
dataf$dataf[[2]]$args[5]
dataf$dataf[[2]]$vals[5]

\dontrun{plot.functional(dataf)}
}
\keyword{datasets}
\keyword{functional}
