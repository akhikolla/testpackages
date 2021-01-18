# Norm.R - normalization of input data

## > system.time( csvs <- Norm$make.weather.dataframe('jma1min', 'Tateno', dir, start.date_, end.date_) )
##    user  system elapsed
##   6.503   0.089   6.708

# Note:
# - convention: when a POSIXct is used to represent a day in a timezone,
#   its hh:mm:ss must be 00:00:00 in the timezone
# - so if there are two POSIXct's that represent days in possibly different timezones,
#   we can check the equality of the days by as.integer(d1 - d2) %/% 24 == 0
# - but taking complete measures againt tz issues appears to be a nightmare..
Norm$timezone <- 'Japan'
Norm$date <- function(year, month, day, hour = 0, min = 0) {
  s <- sprintf('%04d-%02d-%02d %02d:%20d', year, month, day, hour, min)
  as.POSIXct(strptime(s, '%Y-%m-%d %H:%M', tz = Norm$timezone))
}

# DOC:
# - start.date and end.date: POSIXct at yyyy/mm/dd 00:00 JST
# - csv: <dir>/yyyy/*.csv
# XXX: add cache.dir?
Norm$make.weather.data <- function(src, station, start.date, end.date, dir = NULL) {
  if (src != 'jma1min' && src != 'web10min' && src != 'web60min')
    stop('Unsupported input src: ', src)
  if (src == 'jma1min' && is.null(dir)) stop('Data dir missing')

  diff <- as.integer(end.date - start.date) # ndays - 1
  # 1 csv read into a df is ~40KB at 1min sampling, so keeping the list of them
  # is not an issue.
  dfs <- if (src == 'jma1min')
           lapply(0:diff,
                  function(i) {
                    date <- start.date + i*86400 # POSIX ignores leap secs
                    # bind time (min) from 1970-01-01 00:00:00 JST (-9h)
                    # (not the offset from start.date)
                    beg <- as.integer(date)/60 # min
                    cbind(time = beg:(beg+1439),
                          Jma$read.csv(station, date, dir))
                  })
         else if (src == 'web10min')
           lapply(0:diff,
                  function(i) {
                    date <- start.date + i*86400 # POSIX ignores leap secs
                    beg <- as.integer(date)/60 # min
                    cbind(time = seq(beg, beg+1439, 10),
                          Jma$read.from.web(station, date, data.step = 10))
                  })
         else if (src == 'web60min')
           lapply(0:diff,
                  function(i) {
                    date <- start.date + i*86400 # POSIX ignores leap secs
                    beg  <- as.integer(date)/60 # min
                    cbind(time = seq(beg, beg+1439, 60),
                          Jma$read.from.web(station, date, data.step = 60))
                  })
         else stop('cannot happen')
  data <- do.call(`rbind`, dfs)
  data
}
