# This file is part of SMITIDvisu package.
# Copyright (C) 2018-2019 Jean-François Rey <jean-francois.rey@inra.fr>
#
# SMITIDvisu is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# SMITIDvisu is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with SMITIDvisu. If not, see <https://www.gnu.org/licenses/>.
#


#===============================================================================
# * Functions to convert a data frame to a GeoJSON string
#===============================================================================

#-------------------------------------------------------------------------------
# * Main functions
#-------------------------------------------------------------------------------

#' df2geojson
#' @description Transform a data frame into a string formatted in GeoJSON
#' @param df Data frame to convert in GeoJSON.
#' It must contain at least columns 'id', 'time', 'X' and 'Y'.
#' Additionnal columns will be added as features' properties.
#' @param multipleValuesByTime Vector of strings indicating the df columns names which can contain several values by time.
#' @return a geojson string
#' @examples
#' library(SMITIDvisu)
#' data(transmissiontree)
#' geojson <- df2geojson(tt.events, multipleValuesByTime = c('infectedby', 'probabilities'))
#' @export
df2geojson <- function(df, multipleValuesByTime = c()) {
  # cat('===== DF2GEOJSON ======\n')
  # Check if the df contains the necessary columns & these columns are not empty
  if (!("id" %in% names(df) && "time" %in% names(df) &&
        "X" %in% names(df) && "Y" %in% names(df))) {
    # cat('> Empty df\n')
    geojson <- formatGeoJsonCollection('')
  } else {
    # Clean & sort data
    df <- prepareData(df)
    # Gather all the features as a unique string
    features <- getFeatures(df, multipleValuesByTime)
    # Get the final GeoJSON string
    geojson <- formatGeoJsonCollection(features)
  }

  geojson
}

# Clean and sort data by id and time
prepareData <- function(df) {
  if (!("id" %in% names(df) && "time" %in% names(df))) {
    return(df)
  }

  # Converts time column into timestamps
  df$time <- formatTimes(df$time)

  # Order rows by id and time
  df <- df[order(df$id, df$time),]

  # Replace NA coordinates (X,Y) with the previous one according to the time attribute
  df$X <- replaceNA(df$X)
  df$Y <- replaceNA(df$Y)

  df
}

# Convert each feature described inside the given data frame, then return all
# of the features as one GeoJSON string.
getFeatures <- function(df, multipleValuesByTime = c()) {
  # Get all the features strings formatted in GeoJSON
  features <- sapply(unique(df$id),
                     FUN = getFeature,
                     df = df,
                     multipleValuesByTime = multipleValuesByTime)
  # Concatenate GeoJSON features into one string
  features <- paste(features, sep='', collapse=',')
}

# Determine a feature's data and convert it to a GeoJSON string, then return it
getFeature <- function(df, featureId, multipleValuesByTime = c()) {
  # cat(sprintf('-------- Feature: %s --------\n', featureId))
  # Get feature's data
  featureData <- getFeatureDataById(df, featureId)
  # Convert each feature's attribute to a string
  featureAttributes <- getFeatureAttributes(featureData, multipleValuesByTime)
  # Surround feature's id if it's a string
  featureId <- wrapFeatureId(featureId)

  # Format feature attributes in GeoJSON
  properties <- formatFeatureProperties(featureId, featureAttributes)
  geometry <- formatFeatureGeometry(featureAttributes)
  feature <- formatGeoJsonFeature(properties, geometry)

  feature
}

# Return feature's data for the given id, without the column 'id'
getFeatureDataById <- function(df, featureId) {
  if (!is.element('id', names(df))) { return(data.frame()) }
  df[which(df$id == featureId), names(df) != 'id']
}

# Return a feature's attributes as a named list of string
getFeatureAttributes <- function(df, multipleValuesByTime = c()) {
  # Generate the feature's attributes for each time
  featureAttributes <- sapply(unique(df$time),
                              FUN = getFeatureAttributesAtTime,
                              df = df,
                              multipleValuesByTime = multipleValuesByTime)
  featureAttributes <- as.data.frame(t(featureAttributes))

  # Concatenate each column strings separately, surrounded with brackets '[]'
  # and separated by commas ','
  featureAttributes <- lapply(featureAttributes, FUN = wrapProperty)
  featureAttributes
}

# Convert a data frame's columns to a list of strings, one string for each
# column. The data frame must have at least columns 'time', 'X' and 'Y'.
getFeatureAttributesAtTime <- function(df, time, multipleValuesByTime = c()) {
  # cat(sprintf('--- Time: %s ---\n', time))
  if (nrow(df) == 0) {
    return(data.frame())
  }

  # Keep the rows matching the given time, without the 'time' column
  featureAtTime <- df[which(df$time == time),
                           names(df) != 'time']

  # Convert the necessary columns into strings
  coordinates <- formatCoordinates(featureAtTime)

  # Determine and format the additionnal columns
  necessaryColumns <- c('time', 'X', 'Y')
  additionnalColumns <- setdiff(names(featureAtTime), necessaryColumns)
  additionnalProperties <- sapply(additionnalColumns,
                                  FUN = handleAdditionnalProperty,
                                  featureAtTime = featureAtTime,
                                  multipleValuesByTime = multipleValuesByTime)

  # Set the attributes strings
  attributes <- cbind(time, coordinates)
  attributes <- c(attributes, unname(additionnalProperties))
  names(attributes) <- c('time', 'coordinates', additionnalColumns)

  attributes
}

# Convert the first non-NA value (only the first) of a df column into a string.
# If there is no value, return an empty string ''.
# If it is a character, return it surrounded by quotation marks '"%s"'.
# Otherwise, return the value as a string.
handleUniqueValueByTimeProperty <- function(propertyName, featureAtTime) {
  if (!is.element(propertyName, names(featureAtTime))) {
    property <- ''
  } else {
    property <- featureAtTime[propertyName]
    property <- unlist(property, use.names=FALSE)
    property <- property[!is.na(property)] # Remove NA values

    if (is.null(property) || length(property) == 0) {
      # print('Unique-value property: null!!!')
      property <- ''
    } else if (is.character(property) || is.factor(property)) {
      # print('Unique-Value property: character!!!')
      property <- sprintf('"%s"', property[1])
    } else {
      # print('Unique-Value property: else!!!')
      property <- sprintf('%s', property[1])
    }
  }

  property
}

# Does as the 'handleUniqueValueByTimeProperty' function, but takes all the
# column values, separated by commas and surrounded with brackets '[]'.
# If there is no value, return empty brackets '[]'
handleMultipleValuesByTimeProperty <- function(propertyName, featureAtTime) {
  if (!is.element(propertyName, names(featureAtTime))) {
    property <- '[]'
  } else {
    property <- featureAtTime[propertyName]
    property <- unlist(property, use.names=FALSE)
    property <- property[!is.na(property)] # Remove NA values

    if (is.null(property) || length(property) == 0) {
      # print('Multiple-Values property: null!!!')
      property <- '[]'
    } else if (is.character(property) || is.factor(property)) {
      # print('Multiple-Values property: character!!!')
      property <- wrapProperty(property, wrap = '["%s"]', sep = '","')
    } else {
      # print('Multiple-Values property: else!!!')
      property <- wrapProperty(property, wrap = '[%s]', sep = ',')
    }
  }

  property
}

# Convert a Data Frame column into a string
handleAdditionnalProperty <- function(propertyName, featureAtTime, multipleValuesByTime = c()) {
  isUniqueValue <- (is.null(multipleValuesByTime) ||
                    !is.element(propertyName, multipleValuesByTime))

  if (isUniqueValue) {
    property <- handleUniqueValueByTimeProperty(propertyName, featureAtTime)
  } else {
    property <- handleMultipleValuesByTimeProperty(propertyName, featureAtTime)
  }

  property
}

#-------------------------------------------------------------------------------
# * Format strings as GeoJSON
#-------------------------------------------------------------------------------

# Return a collection of GeoJSON features
formatGeoJsonCollection <- function(features) {
  sprintf('{"type":"FeatureCollection","features":[%s]}', paste(features))
}

# Return a GeoJSON feature as a string
formatGeoJsonFeature <- function(properties, geometry) {
  sprintf('{"type":"Feature","properties":{%s},"geometry":{%s}}', properties, geometry)
}

# Return a GeoJSON feature's geometry as a string
formatFeatureGeometry <- function(featureAttributes) {
  sprintf('"type":"LineString","coordinates":%s', featureAttributes$coordinates)
}

# Return a GeoJSON feature's properties as a string
formatFeatureProperties <- function(id, featureProperties) {
  minimalProperties <- sprintf('"id":%s,"linestringTimestamps":%s',
                               id, featureProperties$time)
  additionnalProperties <- formatAdditionnalProperties(featureProperties)
  properties <- sprintf('%s%s', minimalProperties, additionnalProperties)
  properties
}

# Concatenate and return the additionnal properties (if any) as a string,
# starting with a comma ','.
# If there are none, return an empty string ''.
formatAdditionnalProperties <- function(featureAtTime) {
  necessaryColumns <- c('time', 'coordinates')
  additionnalColumns <- setdiff(names(featureAtTime), necessaryColumns)

  if (length(additionnalColumns) == 0) {
    additionnalProperties <- ''
  } else {
    additionnalProperties <- tapply(featureAtTime[additionnalColumns], additionnalColumns, FUN = function(property) {
      return(sprintf('"%s":%s', names(property), property[[1]]))
    })
    additionnalProperties <- paste(additionnalProperties, collapse = ',', sep = '')
    additionnalProperties <- sprintf(',%s', additionnalProperties)
  }

  additionnalProperties
}

#-------------------------------------------------------------------------------
# * Other formatting methods
#-------------------------------------------------------------------------------

# Wrap and return a feature's id inside quotes ("") if it's a factor level or a
# string. Otherwise, return it as is.
wrapFeatureId <- function(id) {
  if (is.factor(id)) {
    # cat('host id is factor\n')
    id <- sprintf('"%s"', levels(id)[id])
  } else if (is.character(id)) {
    # cat('host id is character\n')
    id <- sprintf('"%s"', id)
  } else {
    # cat('convert host id to character\n')
  }

  id
}

# Convert vector values into a string, and wrap and separate them with the
# given strings. If there are no values, return empty brackets '[]'.
wrapProperty <- function(property, wrap = '[%s]', sep = ',') {
  if (is.null(property) || length(property) == 0) {
    return('[]')
  }
  sprintf(wrap, paste(property[!is.na(property)], sep = '', collapse = sep))
}

# Convert a vector of times represented as dates, timestamps, julian days or
# strings, and return them as characters reprensting timestamps in ms.
# Julian days are refered to the current day.
# If times are received as strings but contain only digits, convert them to
# numeric, for them to be processed as such.
formatTimes <- function(times) {
  if (is.null(times) || length(times) == 0) {
    return(c())
  }

  # Convert times to numbers if they are string and contain only digits
  if (!grepl("\\D", times[1])) {
    times <- as.numeric(times)
  }

  if (!is.timestamp(times[1])) {
    if (is.juliendate(times[1])) {
      # Calculate julien days refering to the current day
      times <- as.Date(times, origin = Sys.Date())
    }
    times <- as.numeric(as.POSIXct(times, origin = "1970-01-01", tz = "GMT"))
  }

  # R have timestamps in seconds, so we convert them in ms for js
  times <- sapply(times, FUN = "*", 1000)
  # Prevent R from using scientific notation
  times <- format(times, scientific = FALSE)

  times
}

# Convert coordinates into a string
formatCoordinates <- function(featureAtTime) {
  coordinates <- paste(featureAtTime$X[1], featureAtTime$Y[1], sep = ",", collapse = '')
  coordinates <- sprintf('[%s]', coordinates)
}

#-------------------------------------------------------------------------------
# * Utils
#-------------------------------------------------------------------------------

is.timestamp <- function(time) {
  if (!is.numeric(time)) return(FALSE)
  formattedTime <- as.numeric(format(as.Date(as.POSIXct(as.numeric(time), origin="1970-01-01", tz = "GMT")), format="%Y"))
  return(formattedTime >= 1971)
}

is.juliendate <- function(time) {
  return(is.numeric(time) && !is.timestamp(time))
}

# Replace NA values inside the given vector with the previous value
replaceNA <- function(vector) {
  if (length(vector) == 0) { return(vector) }

  old = 0

  for (i in 1:length(vector)) {
    if (is.na(vector[i])) {
      vector[i] = old
    } else {
      old = vector[i]
    }
  }

  vector
}
