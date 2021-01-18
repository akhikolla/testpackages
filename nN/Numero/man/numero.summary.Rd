\name{numero.summary}
\alias{numero.summary}

\title{Summarize subgroup statistics}

\description{Estimates subgroup statistics after self-organizing map
analysis}

\usage{
numero.summary(results, topology, data = NULL, capacity = 10)
}

\arguments{

\item{results}{
A list object that contains the self-organizing map and its
statistical colorings.
}

\item{topology}{
A SOM topology with additional labels that indicate selected regions.
}

\item{data}{
A matrix or a data frame.
}

\item{capacity}{
Maximum number of subgroups to compare.
}
}

\details{
The input \code{results} must contain the output from
\code{\link{numero.evaluate}()} or similar.

The input argument \code{topology} must be a definition of a SOM with
additional columns as in the output from \code{\link{numero.subgroup}()}.

The function first looks for row names in \code{data} that are also included
in \code{results}. The rows are then divided into subgroups according to the
district assignments in \code{results} and the region labels in
\code{topology}. 
}

\value{
A data frame of summary statistics, see \code{\link{nroSummary}()}
for details. The data frame also contains additional information on which
variables were used for the training of the SOM.

The attribute 'layout' is added to the output. It indicates the location
on the map and the subgroup name and label for each data row that were
included in the analysis.
}

\author{
Ville-Petteri Makinen
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

# Define subgroups.
x <- stats$planes[,"uALB"]
tops <- which(x >= quantile(x, 0.75, na.rm=TRUE))
bottoms <- which(x <= quantile(x, 0.25, na.rm=TRUE))
elem <- data.frame(stats$map$topology, stringsAsFactors = FALSE)
elem$REGION <- "MiddleAlb"
elem$REGION[tops] <- "HighAlb"
elem$REGION[bottoms] <- "LowAlb"
elem$REGION.label <- "M"
elem$REGION.label[tops] <- "H"
elem$REGION.label[bottoms] <- "L"

# Compare subgroups.
cmp <- numero.summary(results = stats, topology = elem, data = dataset)
}