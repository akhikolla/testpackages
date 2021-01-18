#' Plot phenotypic value for a population.
#'
#' Plot the phenotypic value for a population over the course of a
#' prior simulation run.
#'
#' The plot is a line graph depicting the mean, minimum and maximum
#' phenotypic value in the population across generations. This method
#' can only be used if the population has been run via the simulator.
#'
#' @param x an object of class \code{'Population'} which has been run
#'   in the simulator
#' @param ... additional parameters (ignored)
#'
#' @return A plot of the population's simulation run is displayed.
#' @author Dion Detterer, Paul Kwan, Cedric Gondro
#' @export
#'
#' @examples
#' # Build a population
#' pop <- Population(
#'   popSize = 100, map = map100snp, QTL = 20,
#'   alleleFrequencies = runif(100), broadH2 = 0.9,
#'   narrowh2 = 0.5, traitVar = 40
#' )
#' pop <- addEffects(pop)
#' pop <- attachEpiNet(pop)
#' 
#' # Run population in simulation
#' pop <- runSim(pop)
#' 
#' # Plot population's run
#' plot(pop)
#' @seealso \code{\link{Population}}, \code{\link{runSim}},
#'   \code{\link{addEffects}}, \code{\link{attachEpiNet}}
plot.Population <- function(x, ...) {
  testPop(x)

  if (!x$hasRun) {
    stop("Population must be run prior to plotting")
  }

  # Plot minimum, maximum and mean phenotypic values
  summ <- x$summaryData

  summ <- data.frame(Minimum = summ[, 1], Mean = summ[, 4], Maximum = summ
  [
    ,
    6
  ])
  summ$Generation <- 1:nrow(summ)
  summ <- reshape2::melt(summ, id = "Generation")
  names(summ)[2] <- "Range"
  names(summ)[3] <- "Phenotype"
  titstr <- paste("Simulation run across", nrow(x$summaryData), "generations")
  Generation <- NULL
  Phenotype <- NULL
  Range <- NULL
  ggplot2::ggplot(data = summ, ggplot2::aes(
    x = Generation, y = Phenotype,
    color = Range
  )) + ggplot2::geom_line() + ggplot2::ggtitle(titstr)
}


#' Plot epistatic network.
#'
#' Plot an epistatic network between a set of QTLs.
#'
#' An object of class \code{EpiNet} is typically first retrieved from
#' a \code{Population} object (using \code{\link{getEpiNet}}) before
#' being plotted using \code{plot.EpiNet()}.
#'
#' @param x an object of class \code{'EpiNet'}.
#' @param ... additional parameters (ignored)
#'
#' @return A plot of the epistatic network is displayed.
#' @author Dion Detterer, Paul Kwan, Cedric Gondro
#' @export
#'
#' @examples
#' # Build a population with an epistatic network attached
#' pop <- Population(
#'   popSize = 100, map = map100snp, QTL = 20,
#'   alleleFrequencies = runif(100), broadH2 = 0.9,
#'   narrowh2 = 0, traitVar = 40
#' )
#' pop <- attachEpiNet(pop)
#' 
#' # Retrieve and plot the epistatic network
#' epinet <- getEpiNet(pop)
#' plot(epinet)
#' @seealso \code{\link{Population}}, \code{\link{attachEpiNet}},
#'   \code{\link{getEpiNet}}
plot.EpiNet <- function(x, ...) {
  if (!is(x, "EpiNet")) {
    stop("Object must be of EpiNet class")
  }

  net <- splitMatrix(x$Incidence)
  n <- nrow(net[[1]])
  xx <- numeric(0)
  count <- n

  for (i in 1:length(net)) {
    intorder <- sum(net[[i]][, 1])

    if (intorder == 2) {
      for (j in 1:ncol(net[[i]])) xx <- c(xx, which(net[[i]][, j] ==
          1))
    } else {
      for (j in 1:ncol(net[[i]])) {
        count <- count + 1
        for (k in which(net[[i]][, j] == 1)) xx <- c(xx, k, count)
      }
    }
  }

  gg <- igraph::make_graph(xx, directed = FALSE, count)
  igraph::V(gg)$size <- log(igraph::degree(gg) + 1.5) * 500 * sqrt(n / 50) / n / (sum(net[[length(net)]]
  [
    ,
    1
  ]) - 1)
  igraph::V(gg)$color <- 2
  if (count > n) {
    igraph::V(gg)$size[(n + 1):count] <- 0
    igraph::V(gg)$color[(n + 1):count] <- "black"
  }

  edgecount <- 0
  for (i in 1:length(net)) {
    intorder <- sum(net[[i]][, 1])
    edges <- ncol(net[[i]])

    if (intorder > 2) {
      edges <- edges * intorder
    }

    if (i == 1) {
      igraph::E(gg)$color <- intorder
    } else {
      igraph::E(gg)$color[(edgecount + 1):(edgecount + edges)] <- intorder
    }

    edgecount <- edgecount + edges
  }

  igraph::plot.igraph(gg, vertex.label = NA)
  title(main = paste(
    "Epistatic interations between", nrow(x$Incidence),
    "QTLs"
  ))
}
