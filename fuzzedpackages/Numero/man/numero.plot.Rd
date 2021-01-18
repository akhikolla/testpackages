\name{numero.plot}
\alias{numero.plot}

\title{Plot results from SOM analysis}

\description{Plot map colorings and save them as vector graphics}

\usage{
numero.plot(results, variables = NULL, topology = NULL, folder = NULL,
            prefix = "figure", reference = NULL, subplot = NULL,
	    gain = 1, detach = FALSE, capacity = 500, font = NULL)
}

\arguments{
\item{results}{
A list object that contains the self-organizing map and its
statistical colorings.
}
\item{variables}{
A string vector that contains names of variables to show.
}
\item{topology}{
The topology of a SOM with subgroup labels.
}
\item{folder}{
Folder path for saving figures.
}
\item{prefix}{
Prefix for each figure file (if saving enabled).
}
\item{reference}{
Reference color ranges and scales.
}
\item{gain}{
Modifier for overall color intensity.
}
\item{subplot}{
A two-element vector that sets out the number of rows and
columns for subplots per figure.
}
\item{detach}{
Use detached windows for figures.
}
\item{capacity}{
Maximum number of subplots to show on screen.
}
  \item{font}{
Multiplier to adjust font size for SVG and HTML output.
  }
}

\details{
The input \code{results} must contain the output from
code{\link{numero.evaluate}()} or similar.

The input argument \code{topology} can be the topology of a SOM or with
additional columns as in the output from \code{\link{numero.subgroup}()}.

The input argument \code{reference} follows the output format from
\code{\link{numero.evaluate}()}.

Possible values for \code{detach} include "X11", "aqua", TRUE or FALSE.
Using multiple figures may result in different behaviour in terminal vs.
RStudio instances. The default behaviour is to create detached windows for each
figure when the X11 display server is available (e.g. in Linux). To use
detached windows in Mac, use the value "aqua". Setting \code{detach = TRUE}
will use a more general approach, however, some systems may behave
unpredictably. To create multiple figures that remain docked within the
RStudio work window, set \code{detach = FALSE}.

If a destination folder is provided, all plots are saved in files
without plotting them on screen.
}

\value{
The number of figures that were created.
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
trdata <- numero.prepare(data = dataset, variables = trvars,
    batch = "MALE", confounders = c("AGE", "T1D_DURAT"))

# Create a self-organizing map.
modl <- numero.create(data = trdata)

# Evaluate map statistics for all variables.
stats <- numero.evaluate(model = modl, data = dataset)

# Plot map colorings.
numero.plot(results = stats)
}