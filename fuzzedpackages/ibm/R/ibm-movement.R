# Particle diffusion ------------------------------------------------------

#' @title Brownian diffusion of a set of particles
#' @description This funtions performs a brownian difussion over a set of particles.
#' The dimension is automatically calculated from the number of columns of the object.
#' @param object The positions of the particles, dimension is taken from the number
#' of columns or assumed to be 1 is no columns.
#' @param sd Standard deviation for the gaussian jump, for dynamics models should be
#' set proportional to \code{sqrt(dt)}.
#' @param \dots Additional arguments for different methods.
#' @details This functions apply a brownian diffusion to a set of point coordinates.
#' @export
diffusion = function(object, sd, ...) {
  UseMethod("diffusion")
}

#' @export
diffusion.default = brownian1D

#' @export
diffusion.matrix = function(object, sd, N=NULL, ...) {
  if(length(N)!=1) stop("Only one N must be provided.")
  if(N==0) return(object)
  n = ncol(object)
  if(n<1 | n>3) stop("Dimension for brownian diffusion must be 1, 2 or 3.")
  out = switch (n,
                "1" = brownian1D(object=object, sd=sd, N=N),
                "2" = brownian2D(object=object, sd=sd, N=N),
                "3" = brownian3D(object=object, sd=sd, N=N)
  )
  return(out)
}

#' @export
diffusion.array = function(object, sd, N=NULL, ...) {
  out  = apply(X=object, MARGIN=seq_along(dim(object))[-c(1,2)], FUN=diffusion, sd=sd, N=NULL)
  dim(out) = dim(object)
  return(out)
  
}

#' @export
#' @method diffusion list
diffusion.list = function(object, sd, N=-1, ...) {
  out = mapply(object=object, FUN=diffusion, sd=sd, N=N, SIMPLIFY = FALSE)
  return(out)
}


# Boundaries --------------------------------------------------------------

#' @title Spatial boundary restrictions
#' @description Set spatial restrictions to the domain.
#' @param x The positions of the particles.
#' @param \dots Additional arguments for different methods.
#' @details Boundaries is a generic and methods can be written. The default
#' applies simmetric boundaries (dynamics over a torus) or reflexive barriers.
#' @export
boundaries = function(x, ...) {
  UseMethod("boundaries")
}

#' @export
#' @method boundaries default
boundaries.default = function(x, periodic=TRUE, ...) {
  if(all(is.na(x))) return(x)
  out = if(isTRUE(periodic)) torus(x) else mirror(x)
  return(out)
}

#' @export
#' @method boundaries list
boundaries.list = function(x, periodic=TRUE, ...) {
  out = lapply(X=x, FUN=boundaries, periodic=periodic)
  return(out)
}

torus = function(x) {
  x = ifelse(x<0, x+1, x)
  x = ifelse(x>1, x-1, x)
  return(x)
}

mirror = function(x) {
  x = ifelse(x<0, -x, x)
  x = ifelse(x>1, 2-x, x)
  return(x)
}

