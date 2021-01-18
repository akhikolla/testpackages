#'@title Create dosimetric cluster system
#'
#'@description In order to allow interaction of an spatial a correlation clusters
#'in RLumCarlo, first a dosimetric system needs to be created in a three-dimensional space,
#'which is the purpose of this function.
#'
#'@detail
#'
#'**The creation of the dosimetric space**
#'To date, this function supports only a random distribution of clusters within
#'an arbitrary cube with dimensions running from 0-1 for `x`,`y`, and `z`. Positions
#'of clusters are assigned by sampling from a uniform distribution ([stats:runif]).
#'
#'**The grouping of clusters**
#'Clusters are grouped according their euclidean distance calculated with [stats::dist].
#'The grouping is done by [stats::hclust] and the clusters are further cut
#'using [stats::cutree]. `plot = TRUE` displays the created system.
#'The cluster creation is somewhat arbitrary and it may change in future.
#'To that end, for the moment, there is no deeper scientific connection between
#'the parameters used to cut the cluster tree and the physics attempted to
#'be simulated.
#'
#'@param n [numeric] (*with default*): number of clusters to be created
#'in an arbitrary 3-dimensional cube. x, y, z  distances range between 0 and 1.
#'
#'@param h [numeric] (*with default*): numeric scalar the cut the cluster tree
#'using [stats::cutree]. The number must range between 0 and 1.
#'
#'@param plot [logical] (*with default*): enables/disables plot output
#'
#'@param ... further arguments to be passed to the plot output
#'
#'@return The function returns a [list] of class `RLumCarlo_clusters` consisting
#'of [numeric] vector of cluster groups and a [matrix] of the cluster positions
#'in the arbitrary space. If `plot = TRUE` the system is displayed using
#'[scatterplot3d::scatterplot3d]
#'
#'@seealso [stats::dist], [stats::hclust], [stats::cutree]
#'
#'@section Function version: 0.1.0
#'
#'@author Sebastian Kreutzer, Geography & Earth Sciences, Aberystwyth University (United Kingdom)
#'
#'@examples
#'create_ClusterSystem(n = 10, plot = TRUE)
#'
#'@keywords models data
#'@encoding UTF-8
#'@md
#'@export
create_ClusterSystem <- function(
  n = 100,
  h = 0.5,
  plot = FALSE,
  ...){

  ## assign coordinates
  m <- matrix(data = stats::runif(n[1] * 3), ncol = 3)
  colnames(m) <- c("x", "y", "z")

  ## calculate euclidean distance between all points
  m_dist <- stats::dist(m, diag = TRUE)

  ## group clusters ... for this we have to find out points closely
  ## connected and cut the clusters in five groups
  cl_groups <- stats::cutree(stats::hclust(m_dist), h = h[1])

# Plot output -------------------------------------------------------------
 if(plot){
    ## set plot settings
    plot_settings <- modifyList(
      x = list(
      xlim = c(0, 1),
      ylim = c(0, 1),
      zlim = c(0, 1),
      xlab = "Distance [a.u.]",
      ylab = "Distance [a.u.]",
      zlab = "Distance [a.u.]",
      main = paste0("Cluster system (n = ",n,")"),
      color = khroma::color("smooth rainbow")(34)[cl_groups],
      col.grid="grey",
      pch = 16,
      x.ticklabs = NULL,
      y.ticklabs = NULL,
      z.ticklabs = NULL,
      mtext = paste0("h = ",h," | n_groups = ", max(unique(cl_groups)))
      ),
    val = list(...))

    ## create plot
    scatterplot3d::scatterplot3d(
      x = m,
      color = plot_settings$color,
      pch = plot_settings$pch,
      xlim = plot_settings$xlim,
      ylim = plot_settings$ylim,
      zlim = plot_settings$zlim,
      xlab = plot_settings$xlab,
      ylab = plot_settings$ylab,
      zlab = plot_settings$zlab,
      main = plot_settings$main,
      x.ticklabs = plot_settings$x.ticklabs,
      y.ticklabs = plot_settings$y.ticklabs,
      z.ticklabs = plot_settings$z.ticklabs,
      col.grid = plot_settings$col.grid
    )
    mtext(side = 3, text = plot_settings$mtext)
 }

# Return ------------------------------------------------------------------
  output <- list(cl_groups = cl_groups, m = m)
  class(output) <- "RLumCarlo_ClusterSystem"
  return(output)
}
