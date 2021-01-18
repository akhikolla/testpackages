# GCSM between vectors in pure R
gcsm_r <- function(x, y, rescale = FALSE,
                   xmin = NULL, xmax = NULL,
                   ymin = NULL, ymax = NULL,
                   comp = "si") {
  if (anyNA(x) || anyNA(y)) {
    x[is.na(y)] <- NA
    y[is.na(x)] <- NA
    if (all(is.na(x)))  return(NA)
  }
  if (is.null(xmin)) xmin <- min(x, na.rm = T)
  if (is.null(xmax)) xmax <- max(x, na.rm = T)
  if (is.null(ymin)) ymin <- min(y, na.rm = T)
  if (is.null(ymax)) ymax <- max(y, na.rm = T)
  if (xmin > xmax) stop("xmin > xmax, please reset them!");
  if (ymin > ymax) stop("ymin > ymax, please reset them!");
  if (xmax < min(x, na.rm = T) || xmin > max(x, na.rm = T)) stop("[xmin, xmax] is beyond the range of x!");
  if (ymax < min(y, na.rm = T) || ymin > max(y, na.rm = T)) stop("[ymin, ymax] is beyond the range of x!");
  xymin <- min(xmin, ymin)
  xymax <- max(xmax, ymax)
  if (xymin == xymax) stop("global min equals to global max!");
  maxmmin <- xymax - xymin
  maxpmin <- xymax + xymin
  if (rescale) {
    if (xmin == xmax) x[!is.na(x)] <- 1 else x <- (x - xmin) / (xmax - xmin)
    if (ymin == ymax) y[!is.na(y)] <- 1 else y <- (y - ymin) / (ymax - ymin)
    maxmmin <- maxpmin <- 1
  }

  sdx <- stats::sd(x, na.rm = T)
  sdy <- stats::sd(y, na.rm = T)
  s3 <- if (sdx == 0 && sdy == 0) {
    1
  } else if (sdx == 0 || sdy == 0) {
    0
  } else {
    stats::cor(x, y, use = "na.or.complete")
  }
  d2 <- abs(sdx - sdy) / (maxmmin / 2)
  if (d2 > 1) d2 <- 1
  mx <- mean(x, na.rm = T)
  my <- mean(y, na.rm = T)
  d1 <- if (s3 >= 0) abs(mx - my) / maxmmin else abs(maxpmin - mx - my) / maxmmin
  if (d1 > 1) d1 <- 1
  s2 <- 1 - d2
  s1 <- 1 - d1
  if (comp == "si") return(s1 * s2 * s3)
  if (comp == "s1") return(s1)
  if (comp == "s2") return(s2)
  if (comp == "s3") return(s3)
  stop("comp should be 'si' or 's1', 's2', 's3'!")
}
