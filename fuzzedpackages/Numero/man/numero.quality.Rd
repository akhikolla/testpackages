\name{numero.quality}
\alias{numero.quality}

\title{Self-organizing map statistics}

\description{Assign new data to map districts and
calculate quality measures}

\usage{
numero.quality(model, data = NULL)
}

\arguments{
  \item{model}{
A list object that contains a self-organizing map and a data layout.
  }

  \item{data}{
A matrix or a data frame.
  }
}

\details{
The input argument \code{model} must be in the the output format
as returned by \code{\link{numero.create}()}.
}

\value{
A list with named elements: \code{som} contains the self-organizing map;
\code{layout} contains the district assignments for data points;
\code{planes} contains smoothed district averages of quality measures, see
\code{\link{nroAggregate}()} and \code{\link{nroMatch}()}; the element
\code{ranges} contains the reference ranges to be used in
\code{\link{nroColorize}()}; the element \code{palette} is the name of
the colormap to be used for colorings; and \code{statistics} contains
the output from \code{\link{nroPermute}()}.
}

\examples{
# Import data.
fname <- system.file("extdata", "finndiane.txt", package = "Numero")
dataset <- read.delim(file = fname)

# Set identities and manage missing data.
dataset <- numero.clean(dataset, identity = "INDEX")

# Prepare training variables.
trvars <- c("CHOL", "HDL2C", "TG", "CREAT", "uALB")
trdata <- numero.prepare(data = dataset, variables = trvars)

# Create a self-organizing map.
modl <- numero.create(data = trdata)

# Analyze map quality.
qc <- numero.quality(model = modl)
}