#' Create a new polygon from the coordinates of mapped points
#' @title Create polygon from mapped points
#' @param r integer specifying the region of current interest
#' @param my_end n x 2  list of mapped points, i.e. the points to which the polygon should extend
#' @param poly_name character string to name the new polygon (defaults to "unspecified")
#' @param loop_r boolean indicating whether the points are going in a loop
#' @return \code{SpatialPolygons} object created from the mapped points
#' @importFrom rgeos gIntersects gIsValid gDifference
#' @importFrom raster aggregate
#' @importFrom sp SpatialPoints SpatialPolygons Polygon Polygons
#' @importFrom maptools spRbind
#' @examples
#' new_poly <- make_polygons(r = 5, my_end = mappedPoints, loop_r = FALSE)
#' plot(new_poly)
make_polygons <- function(r, my_end, poly_name = "unspecified", loop_r) {
  #current values
  start_line_r <- reg_info$start_lines[[r]]
  start_lines_coords_r <- reg_info$start_lines_coords[[r]]

  #interpolate along end line
  pts_aug <- interp_new_pts(r, new_pts = my_end, reg_info)
  if (loop_r) {
    new_poly <- indiv_poly(pts = pts_aug, r, loop_r, reg_info,
                           poly_name = "loop_poly")
  } else {
    #interpolate along start line
    pts_aug <- interp_new_pts(r, new_pts = pts_aug, reg_info, end = FALSE)

    # Add points where line connecting point sequence crosses start_line
    pts_aug <- inter_start_line(r, pts = pts_aug, reg_info)

    #identify indices of "transition" points when polygons should be started/ended
    trans <- find_trans(pts_aug, r, reg_info)
    n_trans <- length(trans)

    #Make and add invidual polygons
    first <- TRUE
    for (i in 1:(n_trans - 1)) {
      t1 <- trans[i]; t2 <- trans[i + 1]
      poly_to_add <- indiv_poly(pts = pts_aug, r, loop_r, reg_info, t1, t2,
                                sprintf("poly_%i", i))
      if (!is.null(poly_to_add)) {
        if (first) {
          new_poly <- poly_to_add
          first <- FALSE
        } else {
          new_poly <- aggregate(spRbind(new_poly, poly_to_add))
        }
      }
    }

    #no polygons ever added
    if (first) {
      new_poly <- NULL
    }
  }


  #for reg  8, add back the one pixel that connects with reg at
  #its self-intersection, if the new_poly touches it

  if (r == 8) {
    sq_pts <- rbind(c(-1000, -750), c(-1000, -725), c(-1025, -725),
                    c(-1025, -750))
    square <- SpatialPolygons(list(Polygons(list(Polygon(sq_pts)),
                                            "square")))
    if (gIntersects(square, new_poly)) {
      new_poly <- aggregate(spRbind(new_poly, square))
    }
  }


  #keep only polygon within the region
  if (!is.null(new_poly)) {
    new_poly <- keep_poly(gIntersection(new_poly, reg_info$regions[[r]]))
    if (!is.null(new_poly)) { #can become NULL in previous line
      new_poly@polygons[[1]]@ID <- poly_name
    }
  }


  return(new_poly)
}


#' Interpolate contour points that are very close or on the region boundaries.
#' @title Interpolate along region boundaries
#' @param r integer indicating for which region the contours are being generated
#' @param new_pts coordinates of the contour
#' @param reg_info  a \code{reg_info} list (see documentation for \code{reg_info})
#' @param end indicator determining if the points are being interpolated on
#' the ending coordinates or the starting coordinates. Defaults to \code{TRUE}.
#' @param close how close a point must be to the line to count
#' as being on it, defaults to 12.5
#'
interp_new_pts <- function(r, new_pts, reg_info, end = TRUE, close = 12.5) {
  if (end) {
    bd_r <- reg_info$end_coords[[r]]
  } else {
    bd_r <- reg_info$start_lines_coords[[r]]
  }
  loop_r <- reg_info$loop[[r]]

  #Find indices of new pts within 'close' of a region boundary (to interpolate)
  n_pts <- nrow(new_pts)
  on_bd <- rep(FALSE, n_pts)
  for (i in 1:n_pts) {
    test <- min(apply(bd_r, 1, function(x){get_dist(x, new_pts[i,])}))
    if (test <= close) {
      on_bd[i] <- TRUE
    }
  }

  #no points intersecting with bound line, just return orignal line
  if (!any(on_bd)) {
    rownames(new_pts) <- NULL
    return(new_pts)
  }

  new_pts_interp <- matrix(nrow = 0, ncol = 2)
  #points before start of sequence
  if (on_bd[1]) {
    new_pts_interp <- rbind(new_pts_interp,
                            sec_to_interp(p1 = new_pts[1, ], bd_r = bd_r))
  }

  #typical points
  for (s in 1:(n_pts - 1)) {
    if (!(on_bd[s] & on_bd[s + 1])) {
      new_pts_interp <- rbind(new_pts_interp, new_pts[s, ])
    } else {
      new_pts_interp <- rbind(new_pts_interp,
                              sec_to_interp(p1 = new_pts[s,], p2 = new_pts[s + 1,],
                                            bd_r))
    }
  }

  #points at end of sequence
  if (loop_r) {
    if (!(on_bd[n_pts] & on_bd[1])) {
      new_pts_interp <- rbind(new_pts_interp, new_pts[n_pts,] )
    } else {  #interpolate over loop
      new_pts_interp <- rbind(new_pts_interp,
                              sec_to_interp(new_pts[n_pts, ], new_pts[1, ], bd_r, loop_r))
    }
  } else {
    if (!on_bd[n_pts]) {
      new_pts_interp <- rbind(new_pts_interp, new_pts[n_pts,])
    } else {
      new_pts_interp <- rbind(new_pts_interp,
                              sec_to_interp(p2 = new_pts[n_pts,], bd_r = bd_r))
    }
  }

  rownames(new_pts_interp) <- NULL
  return(new_pts_interp)
}


#' Check if a point crosses a line segment
#' @param pt_to_test numeric vector of length two giving the point to test
#' @param fixed_pts matrix of dimension 2 by 2 giving the line segment to test
pt_line_inter <- function(pt_to_test, fixed_pts) {
  if (any(apply(fixed_pts, 1, function(x){all(x == pt_to_test)}))) {
    #if point touches the point boundary, return FALSE looking for intersections
    return(FALSE)
  } else {
    pt_to_test <- SpatialPoints(matrix(round(pt_to_test), ncol = 2))
    fixed_line <- SpatialLines(list(Lines(Line(fixed_pts), ID = "fixed")))
    return(gIntersects(pt_to_test, fixed_line))
  }
}


#' Interpolate a section of line
#' @param p1 vector of length two giving the coordinates of the first point
#' @param p2 vector length two giving the coordinates of the second point
#' @param bd_r matrix with two columns giving the fixed line on which to interpolate
#' @param loop_r boolean indicating whether the points are going in a loop
#' @details If only p1 is given the point is assumed to be the first in the
#' sequence. If only p2 is given the point is assumed to be the last point
#' in the sequence
sec_to_interp <- function(p1 = NULL, p2 = NULL, bd_r, loop_r = FALSE) {
  stopifnot(!is.null(p1) || !is.null(p2))
  n_pts <- nrow(bd_r)

  #find nearest pts on bd_r for p1 and p2
  if (!is.null(p1)) {
    fix1 <- which.min(((bd_r[, 1] - p1[1])^2 +(bd_r[, 2] - p1[2])^2))
  }
  if (!is.null(p2)) {
    fix2 <- which.min((bd_r[, 1] - p2[1])^2 + (bd_r[, 2] - p2[2])^2)
  }

  #check if closest point on bd_r is before or after p1 and adjust if needed
  if (!is.null(p1)) {
    if (fix1 != n_pts) {
      if (pt_line_inter(p1, bd_r[fix1:(fix1 + 1),])) {
          fix1 <- fix1 + 1
      }
    }
  }

  #check if closest point on bd_r is before or after p2 and adjust if needed
  if (!is.null(p2)) {
    if (fix2 != 1) {
      if (pt_line_inter(p2, bd_r[(fix2 - 1):fix2,])) {
        fix2 <- fix2 - 1
      }
    }
  }

  if (!is.null(p1) & !is.null(p2)) {
    if (fix1 == fix2) { #p1 and p2 are in same 1/2 of line segment of bd_r, no need to match up to points on bd_r
      pts_ret <- p1
    } else if (!(loop_r)) {
      if (fix1 > fix2) {
        pts_ret <- rbind(p1, bd_r[fix2, ])
      } else {
        pts_ret <- rbind(p1, bd_r[fix1:fix2,])
      }

    } else {
      if (fix2 != 1) {
        pts_ret <- rbind(p1, bd_r[c(fix1:n_pts, 1:fix2),])
      } else {
        pts_ret <- rbind(p1, bd_r[c(fix1:n_pts, fix2),])
      }
    }

  } else if (!is.null(p1)) { #p1 only
    if (fix1 > 1) {
      pts_ret <- bd_r[1:fix1,]
    } else {
      pts_ret <- matrix(nrow = 0, ncol = 2)
    }

  } else { #p2 only
    pts_ret <- rbind(bd_r[fix2:n_pts, ])
  }

  #remove first point if matches second and last point if matches
  #second-to-last (occurs, for examplem, when p1 = fix1 or p2 = fix2)
  pts_ret <- matrix(pts_ret, ncol = 2)
  n_ret <- nrow(pts_ret)
  if (n_ret >= 2) {
    if (all(pts_ret[1,] == pts_ret[2,])) {
      pts_ret <- pts_ret[2:n_ret,]
      n_ret <- n_ret - 1
    }
    if (n_ret >= 2) {
      if (all(pts_ret[n_ret - 1, ] == pts_ret[n_ret, ])) {
        pts_ret <- pts_ret[1:(n_ret - 1),]
      }
    }
  }

  rownames(pts_ret) <- NULL
  return(pts_ret)
}

#' Add points where line connecting point sequence crosses start_line
#' @param r region number
#' @param pts matrix with two columns giving coordinates of the points
#' @param reg_info  a \code{reg_info} list (see documentation for \code{reg_info})
inter_start_line <- function(r, pts, reg_info) {
  n_pts <- nrow(pts)
  start_line_r <- reg_info$start_lines[[r]]
  pts_new <- matrix(nrow = 0, ncol = 2)

  for (s in 1:(n_pts - 1)) {
    pts_s <- SpatialPoints(pts[s:(s + 1),])
    on_start_line <- gIntersects(pts_s, start_line_r, byid = T)
    if (all(on_start_line)) { #both already on start line, nothing to do
      pts_new <- rbind(pts_new, pts[s, ])
    } else {
      line_seg <- SpatialLines(list(Lines(Line(pts[s:(s+1),]), ID = "line_seg")))
      inter <- gIntersection(line_seg, start_line_r)
      if (!is.null(inter)) {
        #line segment crosses start line, add the intersection point
        new_inter <- get_coords(inter)
        d_to_inter1 <-   get_dist(new_inter[1,], pts[s,])
        d_to_interN <- get_dist(new_inter[nrow(new_inter),], pts[s,])
        if (d_to_interN > d_to_inter1) {
          to_add <- rbind(pts[s, ], new_inter)
        } else {
          to_add <- rbind(pts[s, ], new_inter[nrow(new_inter):1, ])
        }
        n_add <- nrow(to_add)
        if (all(to_add[1,] == to_add[2,])) { #don't need first point twice
          to_add <- matrix(to_add[2:n_add, ], ncol = 2)
          n_add <- n_add - 1
        }
        if (all(to_add[n_add,] == pts[s + 1,])) { #don't need last point twice
          to_add <- matrix(to_add[1:(n_add - 1),], ncol = 2)
        }
        pts_new <- rbind(pts_new, to_add)
      } else {
        #line segement does not cross start line, nothing to add
        pts_new <- rbind(pts_new, pts[s, ])
      }
    }
  }

  #add back the last point
  pts_new <- rbind(pts_new, pts[n_pts,])

  rownames(pts_new) <- NULL
  return(pts_new)
}


#' Generate an individual polygon from a set of points
#' @param pts matrix with two columns giving the coordinates of the
#' generated points
#' @param r integer specifying the region
#' @param loop_r boolean indicating whether the points are going in a loop
#' @param reg_info  a \code{reg_info} list (see documentation for \code{reg_info})
#' @param t1 index of first transition point under consideration in
#' \code{pts} matrix, NULL if \code{loop_r == TRUE}
#' @param t2  index of second transition point under consideration in
#' \code{pts} matrix, NULL if \code{loop_r == TRUE}
#' @param poly_name string giving name of polygon
#' @importFrom rgeos gContains
indiv_poly <- function(pts, r, loop_r, reg_info, t1 = NULL, t2 = NULL,
                       poly_name = "unspecified") {

  stopifnot(loop_r == is.null(t1))
  stopifnot(loop_r == is.null(t2))

  if (loop_r) {
    poly_pts <- pts
  } else {
    #transition points
    p1 <- pts[t1, ]
    p2 <- pts[t2, ]

    #pts are a line seg from p1 to p2 on start line, no polygon to generate
    p1_to_p2 <- SpatialLines(list(Lines(Line(rbind(p1, p2)), ID = "p1_to_p2")))
    if (gContains(reg_info$start_lines[[r]], p1_to_p2)) {
      return(NULL)
    }

    #points on start line between transition points
    on_start <- sec_to_interp(p1, p2, bd_r = reg_info$start_lines_coords[[r]])

    #combine on_start points with generated points; re-ordered generated
    #points, so points form a closed loop
    poly_pts <- matrix(rbind(on_start, pts[t2:t1,]), ncol = 2)

  }

  #need more than two points to make a polygon
  n_poly_pts <- nrow(poly_pts)
  if (nrow(poly_pts) > 3) {
    poly <- SpatialPolygons(list(Polygons(list(Polygon(poly_pts)),
                                                ID = poly_name)))
  } else {
    return(NULL)
  }

  #identify any self-intersections and overlap with land
  if (suppressWarnings(gIsValid(poly))) {
    poly <- gDifference(poly, land)
    if (is.null(poly)) {
      return(NULL)
    }
  }
  if (!suppressWarnings(gIsValid(poly))) {
    poly <- untwist(poly)
    if (is.null(poly)) {
      return(NULL)
    }
    poly <- gDifference(poly, land)
    if (is.null(poly)) {
      return(NULL)
    }
  }

  poly@polygons[[1]]@ID <- poly_name
  return(poly)
}

#' Find transition points (points at the start/end of polygons)
#' @param pts matrix with two columns giving coordinates of the points
#' @param r integer specifying the region
#' @param reg_info  a \code{reg_info} list (see documentation for \code{reg_info})
#' @param close how close a point must be to the line to count
#' as being on it, defaults to 12.5
find_trans <- function(pts, r, reg_info, close = 12.5) {
  n_pts <- nrow(pts)
  trans <- rep(FALSE, n_pts)
  #identify points 'close' to start line
  for (i in 2:(n_pts - 1)) {
    dist_to_start <- min(apply(reg_info$start_lines_coords[[r]], 1,
                               function(x){get_dist(x, pts[i,])}))
    if (dist_to_start <= close) {
      trans[i] <- TRUE
    }
  }
  trans <- which(trans == TRUE)

  #add in first and last point always
  trans <- unique(c(1, trans, n_pts))

  return(trans)
}
