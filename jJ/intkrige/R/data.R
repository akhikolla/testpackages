#' An interval-valued design ground snow load dataset for Utah.
#'
#' A dataset containing the interval-valued data used in the analysis
#' of Bean et. at (2019). The 415 measurement locations
#' included in the dataset are taken from Bean et. al. (2018).
#'
#' @format A data frame with 415 rows and 8 variables:
#'
#'  \describe{
#'  \item{STATION}{The global historical climatological network (GHCN)
#'  station identifier}
#'  \item{STATION_NAME}{The GHCN station name}
#'  \item{LATITUDE}{Latitude coordinate position}
#'  \item{LONGITUDE}{Longitude coordinate position}
#'  \item{ELEVATION}{Elevation of the measurement location (meters)}
#'  \item{minDL}{The lower endpoint of the interval-valued design snow load as measured in kilopascals (kpa)}
#'  \item{maxDL}{The upper endpoint of the design snow load interval (kpa)}
#'  \item{pointDL}{The original point-valued design snow load from the 2018 Utah Snow Load Study (kpa)}
#'
#'  }
#'
#' @details The interval valued kriging analysis described in Bean et. al. (2019) analyzes this dataset on a
#'  log scale after removing the effect of elevation. An example of predictions using this workflow are
#'  found in the README associated with this package.
#'
#' @references
#' \insertRef{Bean2019-int}{intkrige}
#'
#' \insertRef{Bean2018-report}{intkrige}
"utsnow"

#' An interval-valued design ground snow load dataset for Utah that only
#'   considers depth to load conversions.
#'
#' A dataset containing the interval-valued data used in the analysis
#' of Bean et. at (2019). The 415 measurement locations
#' included in the dataset are taken from Bean et. al. (2018).
#'
#' @format A data frame with 415 rows and 8 variables:
#'
#'  \describe{
#'  \item{STATION}{The global historical climatological network (GHCN)
#'  station identifier}
#'  \item{STATION_NAME}{The GHCN station name}
#'  \item{LATITUDE}{Latitude coordinate position}
#'  \item{LONGITUDE}{Longitude coordinate position}
#'  \item{ELEVATION}{Elevation of the measurement location (meters)}
#'  \item{minDL}{The lower endpoint of the interval-valued design snow load as measured in kilopascals (kpa)}
#'  \item{maxDL}{The upper endpoint of the design snow load interval (kpa)}
#'  \item{pointDL}{The original point-valued design snow load from the 2018 Utah Snow Load Study (kpa)}
#'
#'  }
#'
#' @details The interval valued kriging analysis described in Bean et. al. (2019) analyzes this dataset on a
#'  log scale after removing the effect of elevation. An example of predictions using this workflow are
#'  found in the README associated with this package. Note that this dataset differs from utsnow in that
#'  intervals only consider differences in depth to load conversions.
#'
#' @references
#' \insertRef{Bean2019-int}{intkrige}
#'
#' \insertRef{Bean2018-report}{intkrige}
"utsnow_dtl"

#' An interval-valued design ground snow load dataset for Utah that only
#'   considers depth to load conversions.
#'
#' A dataset containing the interval-valued data used in the analysis
#' of Bean et. at (2019). The 415 measurement locations
#' included in the dataset are taken from Bean et. al. (2018).
#'
#' @format A data frame with 415 rows and 8 variables:
#'
#'  \describe{
#'  \item{STATION}{The global historical climatological network (GHCN)
#'  station identifier}
#'  \item{STATION_NAME}{The GHCN station name}
#'  \item{LATITUDE}{Latitude coordinate position}
#'  \item{LONGITUDE}{Longitude coordinate position}
#'  \item{ELEVATION}{Elevation of the measurement location (meters)}
#'  \item{minDL}{The lower endpoint of the interval-valued design snow load as measured in kilopascals (kpa)}
#'  \item{maxDL}{The upper endpoint of the design snow load interval (kpa)}
#'  \item{pointDL}{The original point-valued design snow load from the 2018 Utah Snow Load Study (kpa)}
#'
#'  }
#'
#' @details The interval valued kriging analysis described in Bean et. al. (2019) analyzes this dataset on a
#'  log scale after removing the effect of elevation. An example of predictions using this workflow are
#'  found in the README associated with this package. Note that this dataset differs from utsnow in that
#'  intervals only consider differences in depth to load conversions. This dataset differs from utsnow_dtl
#'  in that intervals are only calculated at the final step of the analysis: after finding 50 year events
#'  using all 8 depth to load conversion techniques. utsnow_dtl rather created annual intervals, only fitting
#'  distributions to two sets of maximums (the annual lower and upper bounds), rather fitting 8 sets of
#'  maximums on all the depth-to-load conversion types.
#'
#' @references
#' \insertRef{Bean2019-int}{intkrige}
#'
#' \insertRef{Bean2018-report}{intkrige}
"utsnow_dtl2"

#' 30 year trimmed mean daily maximum and minimum temperatures for the
#' Ohio river basin.
#'
#' Intervals are defined by the mean daily maximum and minimum temperatures for
#' the Ohio river basin from January 1, 1988 - December 31, 2018. The 116
#' observations in this dataset all had at least 300 daily observations of
#' temperature in at least 30 of the 31 considered years. The mean was
#' calculated after trimming 10% of the data in the tails to remove the
#' influence of potential outliers.
#'
#' @format A data frame with 161 rows and 8 variables:
#'
#'  \describe{
#'  \item{ID}{The global historical climatological network (GHCN)
#'  station identifier}
#'  \item{NAME}{The GHCN station name}
#'  \item{STATE}{The two-digit designation for the state in which each station resides}
#'  \item{LATITUDE}{Latitude coordinate position}
#'  \item{LONGITUDE}{Longitude coordinate position}
#'  \item{ELEVATION}{Elevation of the measurement location (meters)}
#'  \item{minm}{The 30 year mean daily minimum temperature (tenths of degrees Celsius)}
#'  \item{maxm}{The 30 year mean daily maximum temperature (tenths of degrees Celsius)}
#'
#'  }
"ohtemp"

#' Ohio River Basin map
#'
#' @format A SpatialPolygons objects with the boundaries of
#'  the Ohio river basin. Shapefile obtained from the U.S. Geological Survey
#'  National Watershed Boundary Dataset.
"ohMap"
