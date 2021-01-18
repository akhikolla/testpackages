\name{prediction.plot}
\alias{prediction.plot}
\title{
Plot observed versus predicted values for categorical data.
}
\description{
Plot observed minus predicted responses from a cognitive model for categorical data fit with MPTinR. Values above 0 indicate that there are to many responses in that category compared to the predictions, values below 0 indicate that there are to little responses compared to the predictions.
}
\usage{
prediction.plot(results, model.filename, 
dataset = 1, 
absolute = TRUE, 
spacing = 2, 
axis.labels = NULL, 
ylim, model.type = c("easy", "eqn", "eqn2"), 
args.plot = list(), args.rect = list(), args.box = list(), args.points = list(), 
args.labels = list(), numbers = c("individual", "continuous"), 
pos.numbers = c("plot", "axis"), args.numbers = list(), args.abline = list(), abline)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{results}{
list. Results from \code{\link{fit.mpt}}, \code{\link{fit.model}}, or \code{\link{fit.mptinr}}.
  }
  \item{model.filename}{
Same as in \code{\link{fit.mpt}}.
}
  \item{dataset}{
integer scalar or \code{"aggregated"} defining which dataset to plot.
}
  \item{absolute}{
logical. Should absolute deviations (the default) or G^2 deviations be plotted. See Details.
}
  \item{spacing}{
The spacing between two trees in x-axis ticks.
}
  \item{axis.labels}{
The labels on the x-axis. Default is \code{Tree 1} to code{Tree n}.
}
  \item{ylim}{
the \code{ylim} argument to \code{plot}. If missing taken from data.
}
  \item{model.type}{
Same as in \code{\link{fit.mpt}}.
}
  \item{args.plot}{
\code{list}. Further arguments to \code{\link{plot}}, see details.
}
  \item{args.rect}{
\code{list}. Further arguments to \code{\link{rect}}, see details.
}
  \item{args.box}{
\code{list}. Further arguments to \code{\link{box}}, see details.
}
  \item{args.points}{
\code{list}. Further arguments to \code{\link{points}}, see details.
}
  \item{args.labels}{
\code{list}. Further arguments to \code{\link{axis}}, see details.
}
  \item{numbers}{
character vector (using partial matching) or \code{NULL} indicating where/if to plot numbers. Possible values are \code{"individual"}, \code{"continuous"}, or \code{NULL}. \code{"individual"} will start with \code{1} for the first response category in each tree/item type. \code{"continuous"} will use consecutive numbering matching the column numbers/ position of the data, \code{NULL} will plot no numbers. The default plots \code{"individual"} numbers.
}
  \item{pos.numbers}{
character vector, indicating where to plot the points. Possible values are \code{"plot"} or \code{"axis"} (using partial matching).
}
  \item{args.numbers}{
\code{list}. Further arguments to either \code{\link{text}} (if \code{pos.numbers = "plot"}) or \code{\link{axis}} (if \code{pos.numbers = "axis"}), see details.
}
  \item{args.abline}{
\code{list}. Further arguments to \code{\link{abline}}, see details.
}
  \item{abline}{
logical. Whether to print vertical line at the positions of each point. If missing is set to \code{TRUE} if \code{pos.numbers = "axis"}.
}
}
\details{
This function uses base graphics to produce the plots and calls the following functions in the order given to do so: \code{plot} (produces an empty plot with axes), \code{rect} (produces the shaded area for each tree/ item type), \code{box} (produces another box around the plot), possibly \code{abline} (produces the vertical lines for each point), \code{points} (adds the data points), and depending on the value of \code{numbers} and \code{pos.numbers} either \code{text} (adds the numbers in the plot) or \code{axis} (adds the numbers below the plot).

For all of those functions default values are set but can be changed using the corresponding argument. These argument must be a named list containing arguments to that function (see Examples). Default arguments are:

\itemize{
 \item \code{args.plot}:  \code{list(xlab = "", ylab = "", main = "")}
  \item \code{args.rect}: \code{list(col = "grey", border = "transparent", density = 30, angle = 45)}
  \item \code{args.box}: nothing
  \item{args.points}: \code{list(pch = 1, cex = 2.25)}
  \item{args.labels}: either \code{list(line = -1)} or if plotting \code{numbers} at \code{"axis"}, \code{list(line = 1)}
  \item{args.numbers} (if \code{pos.numbers = "plot"}): \code{list(labels = as.character(numbers.pch), cex = 0.7)} (If you change \code{"labels"} you can define what to plot instead of numbers).
  \item{args.numbers} (if \code{pos.numbers = "axis"}): \code{list(labels = numbers.pch, cex.axis = 0.6, mgp = c(0,0.3,0))} (If you change \code{"labels"} you can define what to plot instead of numbers)
  \item{args.abline}: \code{list(col = "darkgrey")}.
}

If \code{absolute = TRUE} (the default) absolute deviations are plotted (i.e., observed - predicted from the model). If \code{absolute = FALSE} G^2 values are plotted which are computed for all predictions where data is non 0 with: \deqn{2 \times \mbox{data} \times (log(\mbox{data}) - log(\mbox{predictions}))}

}
\value{
Invoked for its side effects, but invisibly returns a \code{list} with the x and y positions for each point.
}
\author{
Henrik Singmann. Thanks to David Kellen for discussion and ideas.
}
\note{
Please report all problems.
}

\seealso{
\code{\link{fit.mpt}}
}
\examples{
\dontrun{
#### using the model and data from Broeder & Schuetz:
data(d.broeder, package = "MPTinR")
m.2htm <- system.file("extdata", "5points.2htm.model", package = "MPTinR")
m.sdt <- "pkg/MPTinR/inst/extdata/broeder.sdt.model"

m.sdt <- system.file("extdata", "broeder.sdt.model", package = "MPTinR")

# fit the 2HTM
br.2htm <- fit.mpt(d.broeder, m.2htm)

# graphical parameters
par(mfrow = c(2,2))
prediction.plot(br.2htm, m.2htm, 4)
prediction.plot(br.2htm, m.2htm, 4, ylim = c(-4, 4), numbers = NULL, 
	args.points = list(pch = 16, cex = 1.5))
prediction.plot(br.2htm, m.2htm, 4, ylim = c(-4, 4), args.plot = list(main = "Dataset 4 - A"),
	abline = TRUE, numbers = "continuous")
prediction.plot(br.2htm, m.2htm, 4, ylim = c(-4, 4), args.plot = list(main = "Dataset 4 - B"),
	pos.numbers = "axis", abline = TRUE, 
	args.numbers = list(mgp = c(3, 0.2, 0), cex.axis = 0.35),
	args.points = list(pch = 4, cex = 1.5))
dev.off()

prediction.plot(br.2htm, m.2htm, "aggregated", axis.labels = unlist(lapply(c(10, 25, 50, 75, 90),
	paste, c("o.o", "o.n"), sep = "")))

# fit the SDT
br.sdt <- fit.model(d.broeder, m.sdt, lower.bound = c(rep(-Inf, 5), 0, 1), upper.bound = Inf)

axis.labels <- unlist(lapply(c(10, 25, 50, 75, 90), paste, c("o.o", "o.n"), sep = ""))
# compare predictions for aggregated data:
par(mfrow = c(2,2))
prediction.plot(br.2htm, m.2htm, "aggregated", ylim = c(-30, 30), 
	args.plot = list(main = "MPT model - absolute"), axis.labels = axis.labels)
prediction.plot(br.sdt, m.2htm, "aggregated", ylim = c(-30, 30),
	args.plot = list(main = "SDT model - absolute"), axis.labels = axis.labels)
prediction.plot(br.2htm, m.2htm, "aggregated", ylim = c(-60, 60), 
	args.plot = list(main = "MPT model - G.squared"), absolute = FALSE, 
	axis.labels = axis.labels, pos.numbers = "axis", args.points = list(pch = 8, cex = 1))
prediction.plot(br.sdt, m.2htm, "aggregated", ylim = c(-60, 60),
	args.plot = list(main = "SDT model - G.squared"), absolute = FALSE, 
	axis.labels = axis.labels, pos.numbers = "axis", args.points = list(pch = 8, cex = 1))

# comparing absoulte and G-sqaured plot with zero counts in cell 2:
par(mfrow = c(2,2))
prediction.plot(br.2htm, m.2htm, 2, ylim = c(-1, 1),
	args.plot = list(main = "MPT model - absolute"))
prediction.plot(br.sdt, m.2htm, 2, ylim = c(-1, 1),
	args.plot = list(main = "SDT model - absolute"))
prediction.plot(br.2htm, m.2htm, 2, ylim = c(-2, 2),
	args.plot = list(main = "MPT model - G.squared"), absolute = FALSE)
prediction.plot(br.sdt, m.2htm, 2, ylim = c(-2, 2),
	args.plot = list(main = "SDT model - G.squared"), absolute = FALSE)
}

}
