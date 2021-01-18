# JmaWeb.R - 10min/60min data from Jma webpage
# - http://www.data.jma.go.jp/obd/stats/etrn/

library('XML')

# station: group (prec), id, type('a1', 's1')
# page.fmt: format: ntabs, nth, heads, cols (col.names)
# cols, col.names
# data.step: 10, 60

# columns of interest (from 140705_GetJMAdata.R)
#
# JMA stations:
#  hourly:                   | 10min:
# ---------------------------+-------------------------
#  1:time of day,            | 1:time of day
#  2:atmosphere local        | 2:atmosphere local
#  3:atmosphere sea level    | 3:atmosphere sea level
#  4:precipitation,          | 4:precipitation
#  5:air temperature,        | 5:air temperature
#  8:relative humidity,      | 6:relative humidity
#  9:wind intensity          | 7:wind intensity
#  11:duration of sunshine   | 11:duration of sunshine
#  12:global solar radiation |
#  13:snow fall              |
#  14:snow depth             |
# amedas:
#  hourly:                | 10min:
# ------------------------+------------------------
#  1:time of day          | 1:time of day
#  2:precipitation        | 2:precipitation
#  3:air temperature      | 3:air temperature
#  4:wind intensit        | 4:wind intensity
#  6:duration of sunshine | 8:duration of sunshine

# Schema:
# - ntabs: num tables in a page
# - nth: the table to look at (1..ntabs)
# - heads: num header lines to skip
# - ncols: num columns (= length(cols), to double check)
# - cols: col names
Jma$page.formats <- list(
  's1' = list(
    '10' = list(ntabs = 2, nth = 2, heads = 1:2, ncols = 11, data.step = 10,
      cols = c('clock.hm', 'pressure.hPa', 'pressure.sea.hPa',
        'precipitation.mm', 'temperature.C', 'humidity.perc',
        'wind.m.s', 'wind.direction', 'wind.max.m.s', 'wind.max.direction',
        'radiation.m')),
    '60' = list(ntabs = 3, nth = 2, heads = 1:1, ncols = 17, data.step = 60,
      cols = c('clock.h', 'pressure.hPa', 'pressure.sea.hPa',
        'precipitation.mm', 'temperature.C', 'dew.point.C',
        'steam.pressure.hPa', 'humidity.perc', 'wind.m.s', 'wind.direction',
        'radiation.h', 'global.irradiance.MJ.m2', 'snow.fall.cm', 'snow.acc.cm',
        'weather', 'coludiness', 'visibility.km'))
    ),
  'a1'= list(
    '10' = list(ntabs = 2, nth = 2, heads = 1:1, ncols = 8, data.step = 10,
      cols = c('clock.hm', 'precipitation.mm', 'temperature.C',
        'wind.m.s', 'wind.direction', 'wind.max.m.s', 'wind.max.direction',
        'radiation.m')),
    '60' = list(ntabs = 2, nth = 2, heads = 1:1, ncols = 8, data.step = 60,
      cols = c('clock.h', 'precipitation.mm', 'temperature.C',
        'wind.m.s', 'wind.direction', 'radiation.h',
        'snow.fall.cm', 'snow.acc.cm'))
    )
  )

# http://www.data.jma.go.jp/obd/stats/data/mdrr/man/remark.html
Jma$f2n <- function(x) {
  fs <- levels(x)
  fs[fs == '--'] <- 0
  fs[fs == ''] <- 0 # XXX
  suppressWarnings(as.numeric(fs)[x])
}
Jma$hm2min <- function(x) {
  hms <- strsplit(as.character(x), ':')
  vapply(hms, function(hm) as.numeric(hm[[1]])*60 + as.numeric(hm[[2]]), 0)
}

Jma$page.fmt <- function(station, data.step) {
  fmt <- Jma$page.formats[[station$type]][[as.character(data.step)]]
  if (fmt$ncols != length(fmt$cols)) stop('inconsistent format: ', fmt)
  fmt
}

Jma$url <- function(station.type, data.type, prec, block, year, month, day) {
  paste('http://www.data.jma.go.jp/obd/stats/etrn/view',
        sprintf('%s_%s.php?prec_no=%s&block_no=%s&year=%s&month=%s&day=%s',
                data.type, station.type, prec, block, year, month, day),
        sep='/')
}

# Jma$slurp <- function(url_) { lines <- readLines(url(url_)) }

Jma$get1 <- function(station, year, month, day, fmt, sel) {
  prec <- station$group
  block <- station$id
  data.type <- if (fmt$data.step == 10) '10min'
               else if (fmt$data.step == 60) 'hourly'
               else stop('Unsupported data.step: ', fmt$data.step)
  url <- Jma$url(station$type, data.type, prec, block, year, month, day)
  # cache pages?
  tabs <- XML::readHTMLTable(url)
  if (length(tabs) != fmt$ntabs) stop('JMA has changed the file format for: ', url)

  tab <- tabs[[fmt$nth]]
  names(tab) <- fmt$cols
  # head <- tab[fmt$heads, sel]
  body <- tab[-fmt$heads, sel]
  body
}

# Modify 'sel' to save more columns
Jma$get.s1.10min <- function(station, year, month, day) {
  fmt <- Jma$page.fmt(station, 10)
  # c(1,7,5,6,2,4,11,3)
  sel <- c('clock.hm', 'wind.m.s', 'temperature.C', 'humidity.perc',
           'pressure.hPa', 'precipitation.mm', 'radiation.m',
           'pressure.sea.hPa')
  tab <- Jma$get1(station, year, month, day, fmt, sel)
  # XXX: invalid vs NA
  df <- data.frame(
    'time.of.day'   = Jma$hm2min(tab$clock.hm) - 10,
    'wind'          = Jma$f2n(tab$wind.m.s),
    'temperature'   = Jma$f2n(tab$temperature.C),
    'humidity'      = Jma$f2n(tab$humidity.perc),
    'atmosphere'    = Jma$f2n(tab$pressure.hPa),
    'precipitation' = Jma$f2n(tab$precipitation.mm),
    'radiation'     = Jma$f2n(tab$radiation.m),
    'atmosphere.sea.level' = Jma$f2n(tab$pressure.sea.hPa)
    )
}

Jma$get.s1.60min <- function(station, year, month, day) {
  fmt <- Jma$page.fmt(station, 60)
  # c(1,9,5,8,2,4,11,3,12,13,14)
  sel <- c('clock.h', 'wind.m.s', 'temperature.C', 'humidity.perc',
           'pressure.hPa', 'precipitation.mm', 'radiation.h',
           'pressure.sea.hPa', 'global.irradiance.MJ.m2', 'snow.fall.cm', 'snow.acc.cm')
  tab <- Jma$get1(station, year, month, day, fmt, sel)
  # XXX: invalid vs NA
  df <- data.frame(
    'time.of.day'   = (Jma$f2n(tab$clock.h) - 1) * 60,
    'wind'          = Jma$f2n(tab$wind.m.s),
    'temperature'   = Jma$f2n(tab$temperature.C),
    'humidity'      = Jma$f2n(tab$humidity.perc),
    'atmosphere'    = Jma$f2n(tab$pressure.hPa),
    'precipitation' = Jma$f2n(tab$precipitation.mm),
    'radiation'     = Jma$f2n(tab$radiation.h) * 60,
    'atmosphere.sea.level'   = Jma$f2n(tab$pressure.sea.hPa),
    'global.solar.radiation' = Jma$f2n(tab$global.irradiance.MJ.m2),
    'snow.fall'              = Jma$f2n(tab$snow.fall.cm),
    'snow.depth'             = Jma$f2n(tab$snow.acc.cm)
    )
}

Jma$get.a1.10min <- function(station, year, month, day) {
  fmt <- Jma$page.formats[[station$type]][['10']]
  # c(1,4,3,2,8)
  sel <- c('clock.hm', 'wind.m.s', 'temperature.C', 'precipitation.mm', 'radiation.m')
  tab <- Jma$get1(station, year, month, day, fmt, sel)
  df <- data.frame(
    'time.of.day' =  Jma$hm2min(tab$clock.hm) - 10,
    'wind' = Jma$f2n(tab$wind.m.s),
    'temperature' = Jma$f2n(tab$temperature.C),
    'precipitation' = Jma$f2n(tab$precipitation.mm),
    'radiation' = Jma$f2n(tab$radiation.m) # XXX: radiation NA -> 0
    )
}

Jma$get.a1.60min <- function(station, year, month, day) {
  fmt <- Jma$page.formats[[station$type]][['60']]
  # c(1,4,3,2,6)
  sel <- c('clock.h', 'wind.m.s', 'temperature.C', 'precipitation.mm', 'radiation.h')
  tab <- Jma$get1(station, year, month, day, fmt, sel)
  # XXX: normalize
  df <- data.frame(
    'time.of.day'   = (Jma$f2n(tab$clock.h) - 1) * 60,
    'wind'          = Jma$f2n(tab$wind.m.s),
    'temperature'   = Jma$f2n(tab$temperature.C),
    'precipitation' = Jma$f2n(tab$precipitation.mm),
    'radiation'     = Jma$f2n(tab$radiation.h) * 60 # XXX: radiation NA -> 0
    )
}

# read from web
Jma$read.from.web <- function(station, date, data.step) {
  d <- as.POSIXlt(date)
  year <- 1900+d$year; month <- 1+d$mon; day <- d$mday

  df <- if (station$type == 's1' && data.step == 10) {
    # 10min_s1 seems to be available only from 2008-06-25
    if  (date < Jma$boundary.date)
      stop('10min_s1 data before 2008-06-25 are unavailable (at the moment).')
    Jma$get.s1.10min(station, year, month, day)
  } else if (station$type == 's1' && data.step == 60)
      Jma$get.s1.60min(station, year, month, day)
    else if (station$type == 'a1' && data.step == 10)
      Jma$get.a1.10min(station, year, month, day)
    else if (station$type == 'a1' && data.step == 60)
      Jma$get.a1.60min(station, year, month, day)
    else stop('Unknown web resource: ', station, ', ', data.step)
}
