#' `Q4` class for storing rotation data as quaternions
#'
#' Creates or tests for objects of class "Q4".
#'
#' Construct a single or sample of rotations in 3-dimensions in quaternion form.
#' Several possible inputs for \code{x} are possible and they are differentiated
#' based on their class and dimension.
#'
#' For \code{x} an n-by-3 matrix or a vector of length 3, the angle-axis
#' representation of rotations is utilized.  More specifically, each quaternion
#' can be interpreted as a rotation of some reference frame about the axis
#' \eqn{U} (of unit length) through the angle \eqn{\theta}.  For each axis and
#' angle the quaternion is formed through
#' \deqn{q=[cos(\theta/2),sin(\theta/2)U]^\top.}{q=[cos(theta/2),sin(theta/2)U]'.}
#' The object \code{x} is treated as if it has rows \eqn{U} and \code{theta} is
#' a vector or angles. If no angle is supplied then the length of each axis is
#' taken to be the angle of rotation theta.
#'
#' For \code{x} an n-by-9 matrix of rotation matrices or an object of class
#' \code{"SO3"}, this function will return the quaternion equivalent of
#' \code{x}.  See \code{\link{SO3}} or the vignette "rotations-intro" for more
#' details on rotation matrices.
#'
#' For \code{x} an n-by-4 matrix, rows are treated as quaternions; rows that
#' aren't of unit length are made unit length while the rest are returned
#' untouched.  A message is printed if any of the rows are not quaternions.
#'
#' For \code{x} a \code{"data.frame"}, it is translated into a matrix of the
#' same dimension and the dimensionality of \code{x} is used to determine the
#' data type: angle-axis, quaternion or rotation (see above). As demonstrated
#' below, \code{is.Q4} may return \code{TRUE} for a data frame, but the
#' functions defined for objects of class \code{'Q4'} will not be called until
#' \code{as.Q4} has been used.
#'
#' @name Q4
#'
#' @param x object to be coerced or tested
#' @param theta vector or single rotation angle; if \code{length(theta)==1}, the
#'   same theta is used for all axes
#' @param ... additional arguments.
#'
#' @format \code{id.Q4} is the identity rotation given by the matrix
#'   \eqn{[1,0,0,0]^\top}{[1,0,0,0]'}.
#'
#' @return
#'   \item{as.Q4}{coerces its object into a Q4 type}
#'   \item{is.Q4}{returns \code{TRUE} or \code{FALSE} depending on whether its
#'     argument satisfies the conditions to be an quaternion; namely it must be
#'     four-dimensional and of unit length}
#'
#' @examples
#' # Pull off subject 1's wrist measurements
#' Subj1Wrist <- subset(drill, Subject == '1' & Joint == 'Wrist')
#'
#'                                ## The measurements are in columns 5:8
#' all(is.Q4(Subj1Wrist[,5:8]))   #TRUE, even though Qs is a data.frame, the rows satisfy the
#'                                #conditions necessary to be quaternions BUT,
#'                                #S3 methods (e.g. 'mean' or 'plot') for objects of class
#'                                #'Q4' will not work until 'as.Q4' is used
#'
#' Qs <- as.Q4(Subj1Wrist[,5:8])  #Coerce measurements into 'Q4' type using as.Q4.data.frame
#' all(is.Q4(Qs))                 #TRUE
#' mean(Qs)                       #Estimate central orientation for subject 1's wrist, see ?mean.Q4
#' Rs <- as.SO3(Qs)               #Coerce a 'Q4' object into rotation matrix format, see ?as.SO3
#'
#' #Visualize the measurements, see ?plot.Q4 for more
#' \donttest{
#'   plot(Qs, col = c(1, 2, 3))
#' }
NULL

#' @rdname Q4
#' @export
as.Q4 <- function(x, ...) {
  UseMethod("as.Q4")
}

#' @rdname Q4
#' @export
as.Q4.default <- function(x, theta = NULL, ...) {
  p <- ncol(x)
  n <- nrow(x)

  if (is.null(p)) p <- length(x)

  if (p == 3) {
    # If input is length 3, x is assumed to be vectors in R^3
    U <- x
    U <- matrix(U, ncol = 3)
    ulen <- sqrt(rowSums(U^2))

    if (is.null(theta)) theta <- ulen %% pi

    ntheta <- length(theta)

    if (nrow(U) != ntheta) {
      if (ntheta == 1)
        theta <- rep(theta, n)
      else
        stop("Number of angles must match number of axes.")
    }

    nonZ <- which(ulen != 0)
    U[nonZ, ] <- U[nonZ, ] / ulen[nonZ]
    q <- cbind(cos(theta / 2), sin(theta / 2) * U)
  } else if (p == 4) {
    # If input has length divisible by 4, data are normalized and made into
    # class "Q4"
    q <- matrix(x, ncol = 4)
    rowLens <- (rowSums(q^2))^0.5
    notRots <- which(abs(rowLens - 1) > 1.0e-7)

    if (length(notRots) > 0) {
      txt <- "Row(s)"
      for (i in notRots)
        txt <- paste(txt, i)
      txt <- paste(txt, "was(were) not quaternions.")
      message(txt)
    }

    q <- q / rowLens
  } else if (p == 9) {
    # If input has 9 columns, q is assumed to be rotations
    q <- as.Q4.SO3(x)
  } else
    stop("Unknown data type. Pease see ?as.Q4 for more details.")

  class(q) <- "Q4"
  q
}

#' @rdname Q4
#' @export
as.Q4.SO3 <- function(x, ...) {
  R <- x
  R <- formatSO3(R)
  theta <- mis.angle(R)
  u <- mis.axis(R)
  x <- as.Q4.default(u, theta)
  x
}

#' @rdname Q4
#' @export
as.Q4.Q4 <- function(x, ...) {
  x
}

#' @rdname Q4
#' @export
as.Q4.data.frame <- function(x, ...) {
  q <- data.matrix(x)
  q <- matrix(q, ncol = 4)
  as.Q4.default(q)
}

#' @rdname Q4
#' @export
is.Q4 <- function(x) {
  if (length(x) == 9)
    FALSE
  else
    apply(x, 1, function(q) {
      sum(q^2) - 1 < 1.0e-9 & length(q) == 4
    })
}

#' @rdname Q4
#' @export
id.Q4 <- as.Q4(c(1, 0, 0, 0))
