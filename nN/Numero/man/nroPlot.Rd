\name{nroPlot}
\alias{nroPlot}
\alias{nroPlot.save}

\title{Plot a self-organizing map}

\description{
Create a graphical interface for selecting subgroups from multiple
map colorings simultaneously.
}

\usage{
nroPlot(topology, colors, labels = NULL, subplot = NULL,
        interactive = FALSE, clear = NULL)

nroPlot.save(file, topology, colors, labels = NULL,
             subplot = NULL, font = 1.0)
}

\arguments{
  \item{topology}{
A data frame with K rows and six or more columns that contain the district
positions of a self-organizing map and optional region assignments.
  }
  \item{colors}{
A character vector with K elements or a K x N matrix of hexadecimal
color codes as strings.
  }
  \item{labels}{
A character vector with K elements or a K x N matrix of district labels.
  }
  \item{subplot}{
A two-element vector that sets out the number of rows and columns for
a grid layout of multiple colorings.
  }
  \item{clear}{
If TRUE, all graphics devices are cleared when the plot is refreshed.
}
  \item{interactive}{
If TRUE, an interactive version of the plot is launched.
  }
  \item{file}{
If non-empty, the figure is saved as an SVG or HTML file
instead of plotting on graphics device.
  }
  \item{font}{
Multiplier to adjust font size for SVG and HTML output.
  }
}

\details{
The input \code{topology} must follow the format from
\code{\link{nroKohonen}()}, but may also contain the columns REGION,
and REGION.label that specify the names for subsets of districts and the
single character labels to be shown on top of those districts. The input
can also be the list object as returned by \code{\link{nroKohonen}()}.

The color input can include the attribute 'contrast' that contains a
binary vector or a matrix of equal size. If an element is set, it means
a dark label or highlight will have better contrast with the background.

The label input can include the attribute 'visible' that contains a
binary vector or a matrix of equal size. If an element is set, it means
a label is visible, otherwise it will not be shown on the map.

Some non-alphanumeric characters are not supported and will be
automatically converted to "_". Too long labels or column names will be
truncated.

The default value for \code{clear} is TRUE to prevent multiple plot windows
from accumulating within the RStudio. If you are running R from the terminal
or using detached devices, setting \code{clear} to FALSE will retain the
current window when refreshing.

If the file name ends with ".html", an interactive HTML document is
produced, otherwise an SVG document is created. We recommend opening the
HTML file with a web browser to select regions on large maps (i.e. when the
R plot window becomes too slow to use). The HTML page allows you to assign
subgroups and to save the results as tab-delimited text.
}

\value{
The main function returns a data frame with K rows that contains the
topology and subgrouping information. The \code{.save} subroutine returns
the number of bytes written in the output file.
}

\examples{
# Import data.
fname <- system.file("extdata", "finndiane.txt", package = "Numero")
dataset <- read.delim(file = fname)

# Detect binary columns.
dataset <- nroPreprocess(dataset, method = "")

# Prepare training data.
trvars <- c("CHOL", "HDL2C", "TG", "CREAT", "uALB")
trdata <- scale.default(dataset[,trvars])

# K-means clustering.
km <- nroKmeans(data = trdata)

# Self-organizing map.
sm <- nroKohonen(seeds = km)
sm <- nroTrain(map = sm, data = trdata)

# Assign data points into districts.
matches <- nroMatch(centroids = sm, data = trdata)

# Select a subset of variables and detect binary data.
vars <- c("AGE", "MALE", "uALB", "CHOL", "DIAB_KIDNEY", "DECEASED")
selected <- nroPreprocess(dataset[,vars], method = "")

# Calculate district averages for seleted variables.
vars <- c("AGE", "MALE", "uALB", "CHOL", "DIAB_KIDNEY", "DECEASED")
planes <- nroAggregate(topology = sm, districts = matches, data = selected)

# Estimate statistics.
stats <- nroPermute(map = sm, districts = matches, data = selected)

# Set visuals.
colrs <- nroColorize(values = planes, amplitudes = stats)
labls <- nroLabel(topology = sm, values = planes)

# Add subgrouping information.
topo <- sm$topology
topo$REGION <- ""
topo$REGION[1:8] <- "Center"
topo$REGION[9:21] <- "Perimeter"

# Add subgroup labels.
topo$REGION.label <- ""
topo$REGION.label[1:8] <- "C"
topo$REGION.label[9:21] <- "P"

# Add subgroup colors.
topo$REGION.color <- ""
topo$REGION.color[1:8] <- "#00f00060"
topo$REGION.color[9:21] <- "#f000f060"

# Plot colorings on screen.
nroPlot(topology = topo, colors = colrs, labels = labls)

# Save colorings in file.
#fn <- "colorings.html"
#n <- nroPlot.save(file = fn, topology = topo,
#    colors = colrs, labels = labls)
#cat(n, " bytes saved in '", fn, "'\n", sep="")
}