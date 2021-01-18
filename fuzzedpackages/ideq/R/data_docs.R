#' Simulated data from a standard IDE model
#'
#' A matrix containing simulated data points on a 20 X 20 grid at 20 time points.
#' The rows correspond to the 400 spatial locations,
#' and the columns, to the 20 time points.
#'
#' @seealso [ide_locations] for a matrix containing the spatial locations
#' corresponding to the rows of `ide_standard` and
#' [ide_spatially_varying] for a similar data set generated from
#' a spatially varying ide model
#' 
#' @format A numeric matrix with 400 rows and 20 columns
#' @source Randomly generated as part of Easton Huch's MS project at BYU
"ide_standard"

#' Simulated data from a spatially varying IDE model
#'
#' A matrix containing simulated data points on a 20 X 20 grid at 20 time points.
#' The rows correspond to the 400 spatial locations,
#' and the columns, to the 20 time points.
#'
#' @seealso [ide_locations] for a matrix containing the spatial locations
#' corresponding to the rows of `ide_spatially_varying` and
#' [ide_standard] for a similar data set generated from
#' a standard ide model
#' 
#' @format A numeric matrix with 400 rows and 20 columns
#' @source Randomly generated as part of Easton Huch's MS project at BYU
"ide_spatially_varying"

#' Spatial locations for IDE data sets
#'
#' A matrix containing the 400 two-dimensional spatial locations of the 
#' following data sets: `ide_standard` and `ide_spatially_varying`.
#' The rows of these two data sets correspond with the rows of `ide_locations`.
#'
#' @seealso [ide_standard]/[ide_spatially_varying], the data sets that
#' correspond to these spatial locations
#' 
#' @format A numeric matrix with 400 rows and 2 columns
#' @source Generated as part of Easton Huch's MS project at BYU
"ide_locations"