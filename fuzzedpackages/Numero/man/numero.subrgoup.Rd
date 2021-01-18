\name{numero.subgroup}
\alias{numero.subgroup}

\title{Interactive subgroup assignment}

\description{Plot self-organizing map colorings and let the user choose
multi-district regions as subgroups}

\usage{
numero.subgroup(results, variables, topology = NULL, reference = NULL,
                gain = 1, detach = FALSE, capacity = 9)
}

\arguments{
\item{results}{
A list object that contains the self-organizing map and its
statistical colorings.
}
\item{variables}{
A string vector that contains names of variables to show on screen.
}
\item{topology}{
A SOM topology or the output from a previous subgrouping session.
}
\item{reference}{
Reference color ranges and scales.
}
\item{gain}{
Modifier for overall color intensity.
}
\item{detach}{
Use a detached window.
}
\item{capacity}{
Maximum number of subplots to show on screen.
}
}

\details{
The input \code{results} must contain the output from
code{\link{numero.evaluate}()} or similar.

The input argument \code{topology} can be the structure of a SOM or with
additional columns as in the output from \code{\link{nroPlot}()}.

The input argument \code{reference} follows the output format from
\code{\link{numero.evaluate}()}.

Setting detach to FALSE will also clear all devices whenever the figure is
refreshed. This may be inconvenient when using R from the terminal,
for example; please see the help page of \code{\link{numero.plot}()} for
using detached window device instead.

If any districts are left unmarked, they are automatically collected
into a subgroup of their own.
}

\value{
A data frame similar to the format returned by \code{\link{nroPlot}()}.
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

# Evaluate map statistics for all variables.
stats <- numero.evaluate(model = modl, data = dataset)

# Define subgroups, uncomment to launch interactive window.
#elem <- numero.subgroup(results = stats, variables = trvars)
}