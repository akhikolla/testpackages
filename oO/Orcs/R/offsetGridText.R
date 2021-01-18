#' Insert offset text annotation into 'trellis' plot
#' 
#' @description
#' This is a wrapper function around \code{Orcs:::calcOffsetGridText} and 
#' \strong{grid}-based text drawing functions (currently including 
#' \code{\link{grid.text}} and \code{\link{grid.stext}}) that automatically adds
#' offset text annotations to a 'trellis' plot.
#' 
#' @param x A \code{numeric} vector containing x coordinates, or a 2-column
#' \code{matrix} containing x and y coordinates.
#' @param y A \code{numeric} vector containing y coordinates, or \code{NULL} 
#' if 'x' is a two-column \code{matrix}.
#' @param labels The text to be written as \code{character}.
#' @param xlim,ylim X and Y-axis limits (\code{c(min, max)}) of the current plot. 
#' If not supplied, limits are automatically calculated from supplied x and y
#' coordinates.
#' @param pos Text position specifier(s) as \code{integer} used by 
#' \code{\link{text}}. If not supplied, optimal text positions will be 
#' determined with respect to neighboring locations using 
#' \code{\link[plotrix]{thigmophobe}}. 
#' @param stext \code{logical}, defaults to \code{FALSE}. If \code{TRUE}, shadow 
#' text will be drawn around 'labels'.
#' @param offset A \code{numeric} offset in normalized parent coordinates
#' ("npc", see \code{\link[grid]{unit}}).
#' @param ... Additional arguments passed to the respective \strong{grid} text 
#' drawing function (depends on 'stext'). 
#' 
#' @author
#' Florian Detsch
#' 
#' @seealso
#' \code{\link[grid]{grid.text}}, \code{\link{grid.stext}}, 
#' \code{\link[plotrix]{thigmophobe}}, \code{Orcs::calcOffsetGridText}.
#' 
#' @examples
#' stopifnot(
#'   require(sf)
#'   , require(latticeExtra)
#'   , require(grid)
#' )
#' 
#' # kilimanjaro peaks
#' peaks = data.frame(Peak = c("Kibo", "Mawenzi", "Shira")
#'                    , Lon = c(37.359031, 37.455061, 37.210408)
#'                    , Lat = c(-3.065053, -3.095436, -3.038222))
#' 
#' coordinates(peaks) = ~ Lon + Lat
#' proj4string(peaks) = "+init=epsg:4326"
#' 
#' # visualization
#' xlim_kili <- c(37.15, 37.55)
#' ylim_kili <- c(-3.25, -2.9)
#' 
#' p = spplot(KiLi[[1]], col.regions = "transparent", colorkey = FALSE, 
#'            xlim = xlim_kili, ylim = ylim_kili,
#'            scales = list(draw = TRUE, y = list(rot = 90)), 
#'            sp.layout = rgb2spLayout(KiLi, quantiles = c(0, 1), alpha = .8)) + 
#'   layer(sp.points(peaks, cex = 1.5, pch = 20, col = "black"))
#' 
#' print(p)
#' 
#' downViewport(trellis.vpname(name = "figure"))
#' offsetGridText(x = coordinates(peaks), labels = peaks$Peak,  
#'                xlim = xlim_kili, ylim = ylim_kili, stext = TRUE, offset = .02,
#'                gp = gpar(fontsize = 16))
#'                                
#' @export offsetGridText
#' @name offsetGridText
offsetGridText <- function(x, y = NULL, labels, xlim = NULL, ylim = NULL, 
                           pos = NULL, stext = FALSE, offset = .02, ...) {

  if (is.matrix(x)) {
    y <- x[, 2]
    x <- x[, 1]
  }
  
  # best label locations (if 'pos' is not supplied)
  int_loc_lbl <- if (is.null(pos)) plotrix::thigmophobe(x, y) else pos
  ch_loc_lbl <- pos2just(int_loc_lbl)
  
  # calculate offset point coordinates
  mat_crd_rel_off <- calcOffsetGridText(x = x, y = y, xlim = xlim, ylim = ylim, 
                                        pos = pos, offset = offset)
  
  for (tmp_cnt in 1:nrow(mat_crd_rel_off)) {
    if (stext) {
      grid.stext(labels[tmp_cnt], 
                 x = grid::unit(mat_crd_rel_off[tmp_cnt, 1], "npc"), 
                 y = grid::unit(mat_crd_rel_off[tmp_cnt, 2], "npc"), 
                 just = ch_loc_lbl[tmp_cnt], ...)
    } else {
      grid::grid.text(labels[tmp_cnt], 
                      x = grid::unit(mat_crd_rel_off[tmp_cnt, 1], "npc"), 
                      y = grid::unit(mat_crd_rel_off[tmp_cnt, 2], "npc"),
                      just = ch_loc_lbl[tmp_cnt], ...)
    }
  }

  return(invisible())
}