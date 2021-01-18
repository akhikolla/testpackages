# Part of the landsepi R package.
# Copyright (C) 2017 Loup Rimbaud <loup.rimbaud@inrae.fr>
#                    Julien Papaix <julien.papaix@inrae.fr>
#                    Jean-François Rey <jean-francois.rey@inra.fr>
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software Foundation, Inc.,i
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
#


#' @title Periodic covariance function
#' @name periodic_cov
#' @description Periodic function used to compute the variance-covariance matrix of the fields of the landscape.
#' @param d a numeric object containing pairwise distances between the centroids of the fields
#' @param range range (half-period of oscillations)
#' @param phi amplitude of the oscillations
#' @details The periodic covariance is defined by exp(-2 * sin(d*pi/(2*range))^2 / phi^2). It is used to generate
#' highly fragmented or highly aggregated landscapes.
#' @return An object of the same type as d.
#' @examples
#' periodic_cov(10, range = 5)
#' @export
periodic_cov <- function(d, range, phi = 1) {
  return(exp(-2 * sin(d * pi / (2 * range))^2 / phi^2))
}


#' @title Inverse logit function
#' @name invlogit
#' @description Given a numeric object, return the invlogit of the values. Missing values (NAs) are allowed.
#' @param x a numeric object
#' @details The invlogit is defined by exp(x) / (1+exp(x)). Values in x of -Inf or Inf return invlogits
#' of 0 or 1 respectively.
#' Any NAs in the input will also be NAs in the output.
#' @return An object of the same type as x containing the invlogits of the input values.
#' @examples
#' invlogit(10)
#' @export
invlogit <- function(x) {
  exp(x) / (1 + exp(x))
}

#' @title Logit function
#' @name logit
#' @description Given a numeric object, return the logit of the values. Missing values (NAs) are allowed.
#' @param x a numeric object containing values between 0 and 1
#' @details The logit is defined by log(x/(1-x)). Values in x of 0 or 1 return logits of -Inf or Inf respectively.
#' Any NAs in the input will also be NAs in the output.
#' @return An object of the same type as x containing the logits of the input values.
#' @examples
#' logit(0.5)
#' @export
logit <- function(x) {
  log(x / (1 - x))
}

#' @title Allocation of cultivars
#' @name multiN
#' @description Algorithm based on latent Gaussian fields to allocate two different types of crops across
#' a landscape.
#' @param d a symmetric matrix of the pairwise distances between the centroids of the fields of the landscape.
#' @param area vector containing field areas.
#' @param prop proportion of landscape surface covered by the second type of crop.
#' @param range range of spatial autocorrelation between fields (must be greater or equal 0). The greater
#' the value of range, the higher the degree of spatial aggregation (roughly, range between 0 and 0.1 for
#' fragmented landscapes, between 0.1 and 0.5 for balanced
#' landscapes, between 0.5 and 3 for aggregated landscapes, and above 3 for highly aggregated landscapes).
#' @param algo the algorithm used for the computation of the variance-covariance matrix of the multivariate
#' normal distribution: "exp" for exponential function, "periodic" for periodic function,
#' "random" for random draw (see details). If algo="random", the parameter range is ignored.
#' @details This algorithm allows the control of the proportions of each type of crop in terms of surface
#' coverage, and their level of spatial aggregation. A random vector of values is drawn from a multivariate
#' normal distribution with expectation 0 and a variance-covariance matrix which depends on the pairwise
#' distances between the centroids of the fields. Two different functions allow the computation of the
#' variance-covariance matrix to allocate crops with more or less spatial aggregation (depending on the
#' value of the range parameter). The exponential function codes for an exponential decay of the spatial
#' autocorrelation as distance between fields increases. The periodic function codes for a periodic fluctuation
#' of the spatial autocorrelation as distance between fields increases. Alternatively, a normal distribution
#' can be used for a random allocation of the types of crops.
#' Next, the two types of crops are allocated to different fields depending on whether
#' the value drawn from the multivariate normal distribution is above or below a threshold. The proportion
#' of each type of crop in the landscape is controlled by the value of this threshold (parameter prop).
#' @return A dataframe containing the index of each field (column 1) and the index (0 or 1) of the type
#' of crop grown on these fields (column 2).
#' @examples
#' \dontrun{
#' d <- matrix(rpois(100, 100), nrow = 10)
#' d <- d + t(d) ## ensures that d is symmetric
#' area <- data.frame(id = 1:10, area = 10)
#' multiN(d, area, prop = 0.5, range = 0.5, algo = "periodic")
#' }
#' @export
multiN <- function(d, area, prop, range = 0, algo = "random") {
  nPoly <- nrow(area)
  d <- d / max(d) ## normalization between 0 and 1
  if (range == 0) { ## security
    range <- 1E-6
  }

  ## Multivariate normal distribution with covariance matrix computed via several algorithms
  if (algo == "exp") {
    covMat <- Exponential(d, range)
    s <- as.vector(mvtnorm::rmvnorm(1, mean = rep(0, nPoly), sigma = as.matrix(covMat), method = "chol")[1, ])
  } else if (algo == "periodic") {
    covMat_tmp <- periodic_cov(d, range)
    ## conversion in a positive-definite matrix
    covMatPD <- nearPD(covMat_tmp, keepDiag = TRUE, ensureSymmetry = TRUE, maxit = 1000)
    if (covMatPD$converge == TRUE) {
      covMat <- covMatPD$mat
    } else {
      print(paste("WARNING: covariance cannot be converted in a positive-definite matrix -- range =", range))
      print("Exponential kernel used instead")
      covMat <- Exponential(d, range)
    }
    s <- as.vector(mvtnorm::rmvnorm(1, mean = rep(0, nPoly), sigma = as.matrix(covMat), method = "chol")[1, ])
    
  } else { ## i.e. algo=="random"
    s <- rnorm(nPoly, mean = 0, sd = 1) ## random allocation
  }

  area$s <- s
  ## Ordered gaussian beam
  ord_S <- order(s)
  area.ord_S <- area[ord_S, ]
  area.ord_S$cumsum <- cumsum(area.ord_S$area)

  ## Threshold for allocation of the second type of crop
  prop.areaTot <- prop * area.ord_S$cumsum[nPoly]
  area.ord_S$dif <- abs(area.ord_S$cumsum - prop.areaTot)
  th_index <- which.min(area.ord_S$dif)
  area.ord_S$type <- as.numeric((1:nPoly) < th_index)

  ## Final landscape
  alloc <- area.ord_S[order(area.ord_S$id), c("id", "type")]

  return(alloc)
  # plot(centroid, col=grey((s-min(s))/diff(range(s))), pch=16, cex=2)
  # plot(centroid, col=(s<0)+1, pch=16, cex=2)
}
