# IO.R

######################################################################
# Classes and functions to load data                                 #
######################################################################

IO$slurp <- function(file, variable = NULL) {
  if (!file.exists(file)) {
    stop('File not found: ', file)
  }
  if (is.null(variable)) {
    # text format
    data <- dget(file)
  } else {
    # RData
    load(file)
    data <- get(variable)
    rm(variable)
  }
  data
}

# DOC:
# - To read from file.attribute.200x, call with variable = NULL
# some values:
# - $times.pickup is represented as the offset from weather$startDate
# - $times.of.day is the time in the day in min [0:1439]
#
# The data are quite wild..
# - $hour ranges from 1 to 24 (not 0 to 23)
# - $min ranges from 1 to 60 (not 0 to 59)
# - there are times like 24:01 (!!) (eg rindex 401, file 'satoAM00')
# - lsd of mins are 1 iff ssd is 0 (01, 10, 20, .., 50, 60)
#   - most of the data have min = 01
# So, probably the fact is:
# - $time is the offset of weather.data
# - hh:01 means: picked up within hh:00-59
# - hh:10 means: picked up within hh:00-09min
# - hh:20 means: picked up within hh:10-19min etc.
#
# - $time is the offset from 2008-05-01 00:00:00 weather$startDate
# - 2008/06/19 24:01 seems to come after 2008/06/19/ 22:01 (by 120min)
# - 2008/08/07 19:60 seems to come after 2008/08/07 19:50 (by 10min)
#
# file.data:
# $time       : num [1:461] 148741 148741 148741 148741 148741 ...
# $year       : int [1:461] 2008 2008 2008 2008 2008 2008 2008 2008 2008 2008 ...
# $month      : int [1:461] 8 8 8 8 8 8 8 8 8 8 ...
# $day        : int [1:461] 12 12 12 12 12 12 12 12 12 12 ...
# $hour       : int [1:461] 7 7 7 7 7 7 7 7 9 9 ...
# $min        : int [1:461] 1 1 1 1 1 1 1 1 1 1 ...
# $age        : int [1:461] 68 68 61 61 54 54 47 47 68 68 ...
# $file       : chr [1:461] "izawa24hr" "izawa24hr" "izawa24hr" "izawa24hr" ...
# $type       : int [1:461] 1 1 1 1 1 1 1 1 1 1 ...
#
# secondaries:
# $date       : POSIXct[1:461], format: "2008-08-12 07:00:00" "2008-08-12 07:00:00" ...
# $times.of.day: num [1:461] 421 421 421 421 421 421 421 421 541 541 ...
# $times.pickup: num [1:461] 148741 148741 148741 148741 148741 ...
# $age.norm   : num [1:461] 0.00428 0.00428 -0.05767 -0.05767 -0.11961 ...

# load attributes of samples
IO$Attribute <- setRefClass(
  Class = 'attribute',

  fields = list(
    data      = 'data.frame', # data
    startDate = 'POSIXct',    
    entries   = 'character'   # entries
  ),

  methods = list(
    initialize = function(file.data, sample = NULL) {
      cat('# Preparing attribute data..')

      entries <<- c('age', 'type')

      if (!is.null(sample)) file.data <- file.data[, sample, drop = FALSE]

      d <- file.data[, c('year', 'month', 'day', 'hour', 'min')]
      times.of.day <- (d$hour * 60 + d$min) %% 1440
      times.pickup <- file.data[, 'time']

      # normalize age
      age.norm <- file.data[['age']]
      if (max(age.norm) > min(age.norm)) {
        age.norm <- age.norm / (max(age.norm) - min(age.norm))
      }
      age.norm <- age.norm - mean(age.norm)

      data <<- data.frame(times.of.day=times.of.day, times.pickup=times.pickup,
                          file.data, age.norm=age.norm)
      
      start <- d[1,]
      startDate <<- as.POSIXct(
        strptime(sprintf('%d/%d/%d', start[[1]], start[[2]], start[[3]]),
                 '%Y/%m/%d') +
          start[[4]] * 3600 + (start[[5]] - file.data$time[[1]]) * 60,
        origin = '1970-1-1')
      
      cat('done.\n')
    }
  )
)

# load meteorological data
# DOC:
# - min ranges [0:59] but hour [0:24]! (but hour=24 only if min=0)
# - a day starts at 00:01 and ends at 24:00 (there is no 00:00 in the data)
# - to read from file.weather.200x, call with variable = 'weather'
#
# file.data slurped from file.weather.2008
# 'data.frame':	264960 obs. of  13 variables:
# $time         : num  1 2 3 4 5 6 7 8 9 10 ...
# $wind         : num  2.1 2.3 2.1 2.6 2.8 2.8 2.6 2.3 1.9 2.1 ...
# $temperature  : num  14.5 14.4 14.4 14.7 14.7 14.4 14.7 14.8 14.6 14.5 ...
# $humidity     : num  90 90 90 91 91 90 91 91 90 90 ...
# $atmosphere   : num  1015 1015 1015 1015 1015 ...
# $precipitation: num  0 0 0 0 0 0 0 0 0 0 ...
# $radiation    : num  0 0 0 0 0 0 0 0 0 0 ...
# secondaries:
IO$weather.entries <- c('wind', 'temperature', 'humidity',
                        'atmosphere', 'precipitation', 'radiation')

IO$Weather <- setRefClass(
  Class = 'Weather',

  fields = list(
    data      = 'data.frame',
    entries   = 'character',  
    data.step = 'integer'     # sampling interval
  ),

  methods = list(
    initialize = function(file.data, entries = IO$weather.entries) {
      cat('# Preparing weather data..')
      
      entries <<- entries
      # DOC: file.data should have at least two rows
      data.step <<- as.integer(file.data$time[[2]] - file.data$time[[1]])

      file.data <- file.data[, entries, drop = FALSE]
      data <<- data.frame(file.data)

      cat('done.\n')
    }
  )
)

# load gene expression data (nsamples * ngenes)
IO$Expression <- setRefClass(
  Class = 'Expression',

  fields = list(
    rawdata = 'matrix',
    entries = 'character' # gene names
  ),

  methods = list(
    initialize = function(data, entries = NULL) {
      cat('# Preparing expression data..')

      rawdata <<- data
      if (!is.null(entries)) rawdata <<- rawdata[, entries, drop = FALSE]
      entries <<- colnames(rawdata)
      cat('done.\n')
    }
  )
)

# load weights (nsamples * ngenes)
IO$Weights <- setRefClass(
  Class = 'Weights',

  fields = list(
    rawdata   = 'matrix',
    entries   = 'character' # gene names
  ),

  methods = list(
    initialize = function(data, entries = NULL) {
      cat('# Preparing weight data..')
      
      rawdata <<- data
      if (!is.null(entries)) rawdata <<- rawdata[, entries, drop = FALSE]
      entries <<- colnames(rawdata)
      cat('done.\n')
    }
  )
)

IO$trivialWeights <- function(samples.n, genes) {
  genes.n <- length(genes)
  weights <- matrix(1, nrow=samples.n, ncol=genes.n)
  colnames(weights) <- genes
  list(rawdata=weights, entries=genes)
}
