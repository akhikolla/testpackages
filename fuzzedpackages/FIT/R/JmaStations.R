# JmaStations.R

# Get prec_no from:
# - http://www.data.jma.go.jp/obd/stats/etrn/select/prefecture00.php
# - prec_no are two digit numbers (44: tokyo, 99: antarctic etc.)
# Then get station/amedas ids, names, etc. from
# - http://www.data.jma.go.jp/obd/stats/etrn/select/prefecture.php?prec_no=NN
#
# Note: Amedas station names are not quite unique..

# JMA stations
# - 1661 stations (including Amedas stations)
# Some anomalies:
# - Fujisan belongs to groups 49 (Yamanashi) and 50(Shizuoka)?
# - names are not quite unique
# - Renamed: 'Tsukuba (Tateno)' to 'Tateno'
# - Renamed: 'Okunikkou (Nikkou)' to 'Okunikkou'

# > names(Jma$stations)
#  [1] "id"         "group"      "name"       "type"       "kanji.name"
#  [6] "kana.name"  "lat.d"      "lat.m"      "lon.d"      "lon.m"
# [11] "height"     "rem1"       "rem2"       "rem3"       "rem4"
# [16] "rem5"
#  
# - a simplified view of Jma$stations would be:
# Jma$stations <- data.frame(
#   id     = c('47646', '0324'),         # observation station id   (block_no)
#   group  = c('40', '40'),              # prefecture number (prec_no)
#   name   = c('Tateno', 'Tsuchiura'),   # for convenience
#   type   = c('s1', 'a1'),              # type: meteorological station or AMEDAS
#   stringsAsFactors = FALSE
#   )
# Jma$stations <- read.csv('Jma.stations.csv', stringsAsFactors = FALSE)
Jma$stations.colspec <-
  c('character', 'character', 'character', 'character', 'character',
    'character', 'numeric', 'numeric', 'numeric', 'numeric',
    'numeric', 'character', 'character', 'character', 'character',
    'character')

Jma$stations <- read.csv('inst/extdata/Jma.stations.csv',
                         colClasses=Jma$stations.colspec,
                         stringsAsFactors = FALSE,
                         fileEncoding="UTF-8")

Jma$station <- function(x) {
  if (is.character(x))
    if (!is.na(suppressWarnings(as.numeric(x))))
      Jma$stations[Jma$stations$id == x, ]
    else
      Jma$stations[Jma$stations$name == x, ]
  else stop('Jma$station(x): x must be either a string (some station ids start from 0).')
}
