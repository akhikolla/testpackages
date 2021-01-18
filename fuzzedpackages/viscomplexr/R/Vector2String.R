# R package viscomplexr - phase portraits of functions in the
# complex number plane
# Copyright (C) 2020 Peter Biber
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>



#' Convert a vector into a comma-separated string
#'
#' A simple utility function that transforms any vector into a single character
#' string, where the former vector elements are separated by commas. This is can
#' be useful, in some circumstances, for feeding a series of constant numeric
#' values to \code{\link{phasePortrait}} (see examples). For most applications
#' we recommend, however, to use \code{\link{phasePortrait}}'s parameter
#' \code{moreArgs} instead.
#'
#' @param vec The (usually real or complex valued) vector to be converted.
#'
#' @return A string, where the former vector elements are separated by commas,
#'   enclosed between "c(" and ")".
#'
#' @family helpers
#'
#' @export
#'
#' @examples
#' # Make a vector of 77 complex random numbers inside the unit circle
#' n <- 77
#' a <- complex(n, modulus = runif(n), argument = 2*pi*runif(n))
#' a <- vector2String(a)
#' print(a)
#'
#'
#' # Use this for portraying a Blaschke product
#' \donttest{
#' # x11(width = 9.45, height = 6.30) # Screen device commented out
#'                                    # due to CRAN test requirements.
#'                                    # Use it when trying this example
#' op <- par(mar = c(1, 1, 1, 1), bg = "black")
#' n <- 77
#' a <- complex(n, modulus = runif(n), argument = 2*pi*runif(n))
#' a <- vector2String(a)
#' FUN <- paste("vapply(z, function(z, a){
#'                     return(prod(abs(a)/a * (a-z)/(1-Conj(a)*z)))
#'                    }, a =", a,
#'              ", FUN.VALUE = complex(1))", sep = "")
#' phasePortrait(FUN, pType = "p", axes = FALSE,
#'               xlim = c(-3, 3), ylim = c(-2.0, 2.0),
#'               nCores = 2) # Max. two cores allowed on CRAN
#'                           # not a limit for your own use
#' par(op)
#'   \dontshow{
#'   # R CMD check: make sure any open connections are closed afterward
#'   foreach::registerDoSEQ()
#'   doParallel::stopImplicitCluster()
#'   }
#' }
#'
#'
vector2String <- function(vec) {

  n    <- length(vec)
  rVec <- paste(vec, collapse = ", ")
  rVec <- paste("c(", rVec, ")", sep = "")
  return(rVec)

}
