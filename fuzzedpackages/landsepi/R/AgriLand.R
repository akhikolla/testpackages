# Part of the landsepi R package.
# Copyright (C) 2017 Loup Rimbaud <loup.rimbaud@inrae.fr>
#                    Julien Papaix <julien.papaix@inrae.fr>
#                    Jean-François Rey <jean-francois.rey@inrae.fr>
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

#' @title Landscape allocation
#' @name AgriLand
#' @description Generates a landscape composed of fields where croptypes
#' are allocated with controlled proportions and spatio-temporal aggregation.
#' @param landscape a spatialpolygon object containing field coordinates.
#' @param Nyears an integer giving the number of simulated cropping seasons.
#' @param rotation_period number of years before rotation of the landscape. There is no rotation if rotation_period=0 or rotation_period=Nyears.
#' @param rotation_sequence a list, each element of the list contains indices of croptypes that are cultivated during a period given by "rotation_period".
#' There is no change in cultivated croptypes if the list contains only one element (e.g. only one vector c(0,1,2), indicating cultivation of
#' croptypes 0, 1 and 2).
#' @param rotation_realloc a logical indicating if a new random allocation of croptypes is performed when the landscape is rotated
#' (FALSE=static allocation, TRUE=dynamic allocation). Note that if rotation_realloc=FALSE, all elements of the list "rotation_sequence"
#' must have the same length, and only the first element of the lists "prop" and "aggreg" will be used.
#' @param prop a list of the same size as "rotation_sequence", each element of the list contains a vector of the proportions (in surface)
#' associated with the croptypes in "rotation_sequence". A single vector can be given instead of a list if all elements of
#' "rotation_sequence" are associated with the same proportions.
#' @param aggreg a list of the same size as "rotation_sequence", each element of the list is a single double indicating the degree of
#' aggregation of the landscape. This double must greater or equal 0; the greater its value, the higher the degree of
#' spatial aggregation (roughly, aggreg between 0 and 0.1 for fragmented landscapes, between 0.1 and 0.5 for balanced
#' landscapes, between 0.5 and 3 for aggregated landscapes, and above 3 for highly aggregated landscapes).
#' A single double can be given instead of a list if all elements of
#' "rotation_sequence" are associated with the same level of aggregation.
#' @param algo the algorithm used for the computation of the variance-covariance matrix of the multivariate
#' normal distribution: "exp" for exponential function, "periodic" for periodic function,
#' "random" for random draw (see details of function multiN). If algo="random", the parameter aggreg is not used.
#' @param croptype_names a vector of croptype names (for legend in graphic).
#' @param graphic a logical indicating if a graphic of the landscape must be generated (TRUE) or not (FALSE).
#' @param outputDir a directory to save graphic
#' @details An algorithm based on latent Gaussian fields is used to allocate two different croptypes 
#' across the simulated landscapes (e.g. a susceptible and a resistant cultivar, denoted as 
#' SC and RC, respectively). This algorithm allows the control of the proportions of each croptype 
#' in terms of surface coverage, and their level of spatial aggregation. 
#' A random vector of values is drawn from a multivariate normal distribution with expectation 0 
#' and a variance-covariance matrix which depends on the pairwise distances between
#' the centroids of the fields. Next, the croptypes are allocated to different fields 
#' depending on whether each value drawn from the multivariate normal distribution is above 
#' or below a threshold. The proportion of each cultivar in the landscape is controlled by the value
#' of this threshold. To allocate more than two croptypes, \code{AgriLand} uses sequentially 
#' this algorithm. For instance, the allocation of three croptypes (e.g. SC, RC1 and RC2) 
#' is performed as follows:
#' \enumerate{
#' \item the allocation algorithm is run once to segregate the fields where the susceptible 
#' cultivar is grown, and
#' \item the two resistant cultivars (RC1 and RC2) are assigned to the remaining candidate 
#' fields by re-running the allocation algorithm.
#' }
#' @return a gpkg (shapefile) containing the landscape structure (i.e. coordinates of field boundaries),
#' the area and composition (i.e. croptypes) in time (i.e. each year) for each field. A png graphic can be generated if graphic=TRUE.
#' @importFrom sf st_as_sf
#' @importFrom sf st_write
#' @include Math-Functions.R graphics.R
#' @references Rimbaud L., Papaïx J., Rey J.-F., Barrett L. G. and Thrall P. H. (2018). Assessing the durability and efficiency of landscape-based strategies to deploy plant resistance to pathogens. \emph{PLoS Computational Biology} 14(4):e1006067.
#' @seealso \link{multiN}, \link{periodic_cov}, \link{allocateLandscapeCroptypes}
#' @examples
#' \dontrun{
#' data(landscapeTEST)
#' landscape <- get("landscapeTEST1")
#' set.seed(12345)
#' ## Generate a mosaic of three croptypes in balanced proportions
#' ## and high level of spatial aggregation
#' AgriLand(landscape,
#'   Nyears = 10,
#'   rotation_sequence = c(0, 1, 2), prop = rep(1 / 3, 3),
#'   aggreg = rep(10, 3), algo = "periodic",
#'   graphic = TRUE, outputDir = getwd()
#' )
#'
#' ## Generate a dynamic mosaic of two croptypes in unbalanced proportions
#' ## and low level of spatial aggregation,
#' ## the second croptype being replaced every 5 years without changing field allocation
#' AgriLand(landscape,
#'   Nyears = 20, rotation_period = 5, rotation_sequence = list(c(0, 1), c(0, 2)),
#'   prop = c(1 / 3, 2 / 3), aggreg = c(0.07, 0.07), algo = "periodic", graphic = TRUE,
#'   outputDir = getwd()
#' )
#'
#' ## Generate a dynamic mosaic of four croptypes in balanced proportions
#' ## and medium level of spatial aggregation,
#' ## with field allocation changing every year
#' AgriLand(landscape,
#'   Nyears = 5, rotation_period = 1, rotation_realloc = TRUE,
#'   rotation_sequence = c(0, 1, 2, 3),
#'   prop = rep(1 / 4, 4), aggreg = 0.25, algo = "exp", graphic = TRUE, outputDir = getwd()
#' )
#' }
#' @export
AgriLand <- function(landscape, Nyears, rotation_period = 0, rotation_sequence = list(c(0, 1, 2)),
                     rotation_realloc = FALSE, prop = list(c(1 / 3, 1 / 3, 1 / 3)), aggreg = list(1),
                     algo = "periodic", croptype_names = c(), graphic = FALSE, outputDir = "./") {
  ## Checks
  if (!is.list(rotation_sequence)) {
    rotation_sequence <- list(rotation_sequence)
  }
  if (!is.list(prop)) {
    prop <- list(prop)
  }
  if (!is.list(aggreg)) {
    aggreg <- list(aggreg)
  }

  ## a "rotation" refers to a set of croptypes to allocate during a given period to the fields of the landscape.
  ## a "cycle" is the period corresponding to the cultivation of all rotations
  ## --> the total number of landscapes is K = Nrotations * Ncycles
  Nrotations <- length(rotation_sequence)
  if (rotation_period == 0 | rotation_period == Nyears) {
    rotation_period <- Nyears + 1
    Nrotations <- 1
  }
  ## need to have Nyears+1 because of C algorithm for host plantation
  Ncycles <- ceiling((Nyears + 1) / (Nrotations * rotation_period))
  K <- Nrotations * Ncycles

  ## Other checks
  if (length(prop) != Nrotations) {
    if (length(prop) == 1) {
      prop <- rep(prop, Nrotations)
    } else {
      stop("Error: rotation_sequence and prop must have the same length")
    }
  }
  if (length(aggreg) != Nrotations) {
    if (length(aggreg) == 1) {
      aggreg <- rep(aggreg, Nrotations)
    } else {
      stop("Error: rotation_sequence and aggreg must have the same length")
    }
  }
  if (!prod(sapply(prop, length) == sapply(rotation_sequence, length))) {
    stop("Error: elements of rotation_sequence must have the same length as elements of prop")
  }
  sameNbElts <- prod(sapply(rotation_sequence, function(x) {
    length(x) == length(rotation_sequence[[1]])
  }))
  if (rotation_realloc == FALSE & sameNbElts == FALSE) {
    stop("Error: all elements of rotation_sequence must have the same number of croptypes since 
         there is no re-allocation")
  }
  sumProp <- prod(sapply(prop, function(x) {
    sum(x) == 1
  }))
  if (sumProp == FALSE) {
    stop("Error: proportions of the croptypes (i.e. elements of list prop) must have a sum at 1")
  }

  Npoly.tmp <- length(landscape)

  ## Centroid of the paddocks
  centroid <- NULL
  area <- NULL
  for (i in 1:Npoly.tmp) {
    for (j in 1:length(landscape@polygons[[i]])) {
      centroid <- rbind(centroid, apply(landscape@polygons[[i]]@Polygons[[j]]@coords, 2, mean))
      area <- c(area, landscape@polygons[[i]]@Polygons[[j]]@area)
    }
  }
  Npoly <- nrow(centroid)
  d <- as.matrix(dist(centroid)) ## 2-by-2 distance between centroid of each paddock
  area.df <- data.frame(id = 1:Npoly, area)

  ## Generation of K landscapes
  croptype_alloc <- matrix(NA, nrow = Npoly, ncol = K * rotation_period)
  rotationID <- rep(1:Nrotations, Ncycles)
  for (k in 1:K) {
    rotation_sequence_k <- rotation_sequence[[rotationID[k]]]

    if (k == 1 | rotation_realloc == TRUE) {
      prop_k <- prop[[rotationID[k]]]
      aggreg_k <- aggreg[[rotationID[k]]]
      nAlloc <- length(rotation_sequence_k) - 1

      ## Transform prop and aggreg into nestedProp and nestedAggreg (for sequencial allocation in algorithm)
      nestedProp <- numeric()
      for (i in 1:nAlloc) {
        nestedProp[i] <- sum(prop_k[(i + 1):length(prop_k)]) / sum(prop_k[i:length(prop_k)])
      }
      nestedAggreg <- rep(aggreg_k, nAlloc)

      ## allocation algorithm
      croptype_tmp <- rep(0, Npoly)
      i <- 0
      fieldsToAllocate <- 1:Npoly
      while (sum(fieldsToAllocate) > 0 & i < nAlloc) {
        ## update area and d for allocation of the cultivar
        area.tmp <- area.df[croptype_tmp == i, ]
        d.tmp <- d[croptype_tmp == i, croptype_tmp == i]
        i <- i + 1
        ## Compute the multivariate normal distribution
        alloc.tmp <- multiN(d.tmp, area.tmp, nestedProp[i], nestedAggreg[i], algo)

        ## update croptype_tmp by allocating i
        fieldsToAllocate <- alloc.tmp$type == 1
        croptype_tmp[alloc.tmp$id[fieldsToAllocate]] <- i
      }
    } ## if realloc
    year <- (k - 1) * rotation_period + 1
    croptype_alloc[, year:(year + rotation_period - 1)] <- rotation_sequence_k[croptype_tmp + 1]
  } ## for k

  croptype.df <- data.frame(croptype_alloc[, 1:(Nyears + 1)])
  colnames(croptype.df) <- paste("year_", 1:(Nyears + 1), sep = "")

  landscape_croptypes <- SpatialPolygonsDataFrame(landscape, croptype.df, match.ID = T)

  results <- st_as_sf(landscape_croptypes)

  ## Graphic representing the landscape
  if (graphic) {
    yearsToPlot <- seq(1, Nyears, by = rotation_period)[1:Nrotations]
    if (rotation_realloc) {
      yearsToPlot <- seq(1, Nyears, by = rotation_period)
    }

    for (year in yearsToPlot) {
      plot_allocation(landscape_croptypes, year,
        croptype_names = croptype_names,
        title = paste("Simulated landscape")
        # , subtitle = paste("Cropping ratios = (", paste(round(prop[[1]], 2), collapse = " ; "), ")",
        #                    "   Aggregation =", paste(names(range)[aggreg[[rotationID[k]]]]))
        , subtitle = paste("Year", year),
        filename = paste(outputDir, "/landscape_year", year, ".png", sep = "")
      )
    }
  }
  return(results)
}
