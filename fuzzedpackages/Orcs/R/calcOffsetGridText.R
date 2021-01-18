# Calculate offset coordinates for grid-based text annotations
# 
# @description
# Calculate offset coordinates for (\strong{grid}-based) text drawing 
# functions, e.g. \code{\link{grid.text}}.
# 
# @param x A \code{numeric} vector containing x coordinates, or a 2-column
# \code{matrix} containing x and y coordinates.
# @param y A \code{numeric} vector containing y coordinates, or \code{NULL} 
# if 'x' is a two-column \code{matrix}.
# @param offset A \code{numeric} offset in normalized parent coordinates
# ("npc", see \code{\link[grid]{unit}}).
# @param pos Text position specifier(s) as \code{integer} used by 
# \code{\link{text}}. If not supplied, optimal text positions will be 
# determined with respect to neighboring locations using 
# \code{\link[plotrix]{thigmophobe}}. 
# @param xlim,ylim X and Y-axis limits (\code{c(min, max)}) of the current plot. 
# If not supplied, limits are automatically calculated from supplied x and y
# coordinates.
# @param ... Currently not used. 
# 
# @return
# A \code{numeric matrix} containing offset coordinates.
# 
# @author
# Florian Detsch
# 
# @seealso
# \code{\link[grid]{grid.text}}, \code{\link{text}}, 
# \code{\link[plotrix]{thigmophobe}}
# 
# @examples
# \dontrun{
# stopifnot(
#   require(mapview)
#   , require(lattice)
#   , require(grid)
# )
# 
# ## calculate offsets for breweries with more than 3 different types of beer
# brw = as(subset(breweries, number.of.types > 3), "Spatial")
# loc = calcOffsetGridText(coordinates(brw), offset = .025)
# 
# ## create plot
# p = spplot(brw, zcol = "number.of.types", auto.key = FALSE)
# 
# plot.new()
# print(p, newpage = FALSE)
# 
# ## add text labels
# downViewport(trellis.vpname(name = "figure"))
# for (i in 1:length(brw)) {
#   grid.text(label = brw$number.of.types[i], x = loc[i, 1], y = loc[i, 2])
# }
# }
#               
# @export calcOffsetGridText
# @aliases calcOffsetGridText
calcOffsetGridText <- function(x, y = NULL, offset = 0.02, pos = NULL, 
                               xlim = NULL, ylim = NULL, ...) {
  
  if (is.matrix(x)) {
    y <- x[, 2]
    x <- x[, 1]
  }
  
  # relative ("npc") pointcoordinates
  num_xmin <- if (is.null(xlim)) min(x) - .04 * (max(x) - min(x)) else xlim[1]
  num_xmax <- if (is.null(xlim)) max(x) + .04 * (max(x) - min(x)) else xlim[2]
  num_xrng <- num_xmax - num_xmin
  num_x_rel <- (x-num_xmin) / num_xrng
  
  num_ymin <- if (is.null(ylim)) min(y) - .04 * (max(y) - min(y)) else ylim[1]
  num_ymax <- if (is.null(ylim)) max(y) + .04 * (max(y) - min(y)) else ylim[2]
  num_yrng <- num_ymax - num_ymin
  num_y_rel <- (y-num_ymin) / num_yrng
    
  # best label locations (if 'pos' is not supplied)
  int_loc_lbl <- if (is.null(pos)) plotrix::thigmophobe(num_x_rel, num_y_rel) else pos
  ch_loc_lbl <- pos2just(int_loc_lbl)
  
  # apply offset to point coordinates
  ls_off <- lapply(1:length(num_x_rel), function(tmp_cnt) {

    tmp_x <- num_x_rel[tmp_cnt]
    tmp_y <- num_y_rel[tmp_cnt]
    
    ch_jst <- ch_loc_lbl[tmp_cnt]
    
    if (ch_jst %in% c("left", "right")) {
      if (ch_jst == "left") {tmp_x <- tmp_x+offset} else {tmp_x <- tmp_x-offset}
    } else {
      if (ch_jst == "top") {tmp_y <- tmp_y-offset} else {tmp_y <- tmp_y+offset*1.5}
    }
    
    tmp_mat <- matrix(c(tmp_x, tmp_y), byrow = TRUE, ncol = 2)

    return(tmp_mat)
  })

  mat_off <- do.call("rbind", ls_off)
  return(mat_off)
}


### pos2just(): convert integer text position to string ----

# Convert integer text position specifier to string
# 
# @description
# Convert integer position specifiers as supported by
# \code{\link{text}} to character position specifiers as supported by
# \strong{grid}-based text drawing functions (e.g. \code{\link{grid.text}}).
# 
# @param pos Integer. A position specifier for text annotations as used by
# \code{\link{text}}.
# @param ... Currently not in use.
# 
# @return
# A character vector used as input for text justification in \strong{grid}-based text
# drawing functions.
# 
# @author
# Florian Detsch
# 
# @seealso
# \code{\link{text}}, \code{\link{grid.text}}
# 
# @examples
# set.seed(100)
# pos <- sample(1:4, 5, replace = TRUE)
# 
# pos2just(pos)
# 
# @export pos2just
# @aliases pos2just
pos2just <- function(pos, ...) {
  
  sapply(pos, function(x) {
    if (x == 1) {
      return("top")
    } else if (x == 2) {
      return("right")
    } else if (x == 3) {
      return("bottom")
    } else if (x == 4) {
      return("left")
    } else {
      stop("Invalid position specifier supplied: ", x)
    }
  })
  
}