# JmaCsv.R - 1min data from CD-R (csv)

# JMA uses new format since this date (inclusive):
# Jma$boundary.date <- Norm$date(2008, 6, 25)
Jma$boundary.date <- as.POSIXct(strptime('2008-06-25', '%Y-%m-%d', 'Japan'))

# Note: station.id is a string (of length 5) rather than a number.
# (Could have 0 padded at print time..)
Jma$filename <- function(station.id, date) { # POSIXct
  d <- as.POSIXlt(date)
  if (date < Jma$boundary.date)
    sprintf('s1m_%s_%04d%02d%02d.csv', station.id, 1900+d$year, 1+d$mon, d$mday)
  else
    sprintf('sfc_1min_%04d%02d%02d_%s.csv', 1900+d$year, 1+d$mon, d$mday, station.id)
}

Jma$read.old.csv <- function(date, file) {
  # a day: 00:01 -- 24:00 (1min windows labeled by the time at the end)
  cols <- c(3,4,5,6,7, 22, 29, 38, 45, 63, 76)
  cols.header <- c('year', 'month', 'day', 'hour', 'min',
                   'wind', 'temperature', 'humidity', 'atmosphere',
                   'precipitation', 'radiation')
  sel <- rep('NULL', 90)
  sel[cols] <- 'numeric'

  # Hmm, this is slow..
  # XXX: shall read.table from pipe('cat ..')?
  csv <- read.table(file, header = FALSE, sep = ',', colClasses = sel,
                    skip = 2, fileEncoding = 'Shift-JIS')
  names(csv) <- cols.header
  file.date <- Norm$date(csv$year[[1]], csv$month[[1]], csv$day[[1]])
  if (nrow(csv) != 1440) stop(file, ' corrupt?')
  if (as.integer(date - file.date) %/% 24 != 0) stop(file, ': ', date, ' != ', file.date)

  # flag(ST) is encoded in the least significant decimal digit of each entry
  # 0: unavailable
  # 1-3: very bad
  # 4-6: bad
  # 7-9: okay
  w.flag <- csv$wind          %% 10
  t.flag <- csv$temperature   %% 10
  h.flag <- csv$humidity      %% 10
  a.flag <- csv$atmosphere    %% 10
  p.flag <- csv$precipitation %% 10
  r.flag <- csv$radiation     %% 10
  csv$wind         [w.flag < 7] <- NA
  csv$temperature  [t.flag < 7] <- NA
  csv$humidity     [h.flag < 7] <- NA
  csv$atmosphere   [a.flag < 7] <- NA
  csv$precipitation[p.flag < 7] <- NA
  csv$radiation    [r.flag < 7] <- NA
  data <- data.frame(
    time.of.day    = as.integer(csv$hour*60 + csv$min - 1), # min
    wind           = (csv$wind-w.flag)/100.0,               # m/s
    temperature    = (csv$temperature-t.flag)/100.0 - 50.0, # celsius
    humidity       = (csv$humidity-h.flag)/10.0,            # %
    atmosphere     = (csv$atmosphere-a.flag)/100.0,         # hPa
    precipitation  = (csv$precipitation-p.flag)/10.0,       # mm
    radiation      = (csv$radiation-r.flag)/10.0            # sec
    )
  if (!all(data$time.of.day == 0:1439)) stop(file, ' contains invalid time info')
  data
}

Jma$read.new.csv <- function(date, file) {
  cols <- c(15,16,17,18,19, 22,23, 50,51, 53,54, 61,62, 73,74, 80,81)
  cols.header <- c('year', 'month', 'day', 'hour', 'min',
                   'precipitation', 'p.flag', # 22,23
                   'wind', 'w.flag',          # 50,51
                   'temperature', 't.flag',   # 53,54
                   'radiation', 'r.flag',     # 61,62
                   'atmosphere', 'a.flag',    # 73,74
                   'humidity', 'h.flag'       # 80,81
                   )
  sel <- rep('NULL', 92)
  sel[cols] <- 'numeric'

  csv <- read.table(file, header = FALSE, sep = ',', colClasses = sel,
                    skip = 1, fileEncoding = 'Shift-JIS')
  names(csv) <- cols.header
  file.date <- Norm$date(csv$year[[1]], csv$month[[1]], csv$day[[1]])
  if (nrow(csv) != 1440) stop(file, ' corrupt?')
  if (as.integer(date - file.date) %/% 24 != 0) stop(file, ': ', date, ' != ', file.date)
  # +-------+-------------+--------+--------------------------------+--------+
  # |0,1    |normal       |numeric |normal quality                  |normal  |
  # +-------+             +--------+                                |        |
  # |2,3    |             |NA      |                                |        |
  # +-------+-------------+--------+--------------------------------+--------+
  # |8,9    |doubtful     |numeric |quasi-normal                    |        |
  # +-------+             +--------+                                |        |
  # |10,11  |             |NA      |                                |        |
  # +-------+-------------+--------+--------------------------------+--------+
  # csv$wind         [csv$w.flag!=0 & csv$w.flag!=1] <- NA
  csv$wind         [csv$w.flag >= 4] <- NA
  csv$temperature  [csv$t.flag >= 4] <- NA
  csv$humidity     [csv$h.flag >= 4] <- NA
  csv$atmosphere   [csv$a.flag >= 4] <- NA
  csv$precipitation[csv$p.flag >= 4] <- NA
  csv$radiation    [csv$r.flag >= 4] <- NA

  data <- data.frame(
    time.of.day    = as.integer(csv$hour*60 + csv$min - 1), # min
    wind           = csv$wind*10.0,                         # m/s? XXX
    temperature    = csv$temperature/10.0,                  # celsius
    humidity       = csv$humidity,                          # %
    atmosphere     = csv$atmosphere/10.0,                   # hPa
    precipitation  = csv$precipitation,                     # mm
    radiation      = csv$radiation                          # sec
    )
  if (!all(data$time.of.day == 0:1439)) stop(file, ' contains invalid time info')
  data
}

# date here should be POSIXct
Jma$read.csv <- function(station, date, dir) {
  file <- paste(dir, 1900+as.POSIXlt(date)$year, Jma$filename(station$id, date), sep='/')
  if (date < Jma$boundary.date)
    Jma$read.old.csv(date, file)
  else
    Jma$read.new.csv(date, file)
}
