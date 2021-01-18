#####################################################
# Name: heatmap.R                                   #
# Author: Poppy Miller <p.miller@lancaster.ac.uk> & #
# Chris Jewell <c.jewell@lancaster.ac.uk>           #
# Created: 20161206                                 #
# Copyright: Poppy Miller & Chris Jewell 2016       #
# Purpose: Draws a clustered heatmap                #
#####################################################
are_colours <- function(object) {
  sapply(object, function(x) {
    tryCatch(
      is.matrix(grDevices::col2rgb(x)),
      error = function(e)
        FALSE
    )
  })
}

clusterHeatMap <- function(object, cols, xnames = 1:length(object), hclust_method) {

  # Check colours
  if (length(cols) != 2 |
      !mode(cols) %in% c("character") | !all(are_colours(cols))) {
    message(
      "The argument cols contain colours that are not valid. The defaults will be used instead."
    )
    cols <- c("blue", "white")
  }

  # compute dissimilarity matrix for the type effect clusters
  disim_clust_g <- cluster::daisy(object)
  clu <-
    stats::hclust(disim_clust_g, hclust_method) # default method is complete
  dend <- stats::as.dendrogram(clu)

  # OPTIONAL: change the colour of the heatmap. The lighter the colour
  # (when using the default white blue colour scheme),
  # the higher the dissimilarity between the 2 types (i.e. the less
  # often two type effects are assigned to the same group in the mcmc)
  hmcols <- grDevices::colorRampPalette(cols)(299)

  heatmap_data <- as.matrix(disim_clust_g)

  rownames(heatmap_data) <- colnames(heatmap_data) <- xnames

  gplots::heatmap.2(
    heatmap_data,
    density.info = "none",
    # turns off density plot in the legend
    trace = "none",
    # turns off trace lines in the heat map
    col = hmcols,
    # use color palette defined earlier
    dendrogram = "col",
    # only draw a row dendrogram
    Colv = dend,
    Rowv = dend,
    symm = TRUE,
    key = F
  )
}
