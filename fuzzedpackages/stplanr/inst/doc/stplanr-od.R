## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  eval = TRUE
)
has_webshot <- "webshot" %in% installed.packages()

## ----setup, message=FALSE-----------------------------------------------------
library(stplanr)
library(dplyr)
od <- stplanr::od_data_sample %>%
  select(-matches("rail|name|moto|car|tax|home|la_")) %>%
  top_n(n = 14, wt = all)
class(od)
od

## -----------------------------------------------------------------------------
od[1:3]

## -----------------------------------------------------------------------------
od_matrix <- od_to_odmatrix(od[1:3])
class(od_matrix)
od_matrix

## -----------------------------------------------------------------------------
lapply(c("all", "bicycle"), function(x) od_to_odmatrix(od[c("geo_code1", "geo_code2", x)]))

## -----------------------------------------------------------------------------
odmatrix_to_od(od_matrix)

## -----------------------------------------------------------------------------
(od_inter <- od %>% filter(geo_code1 != geo_code2))
(od_intra <- od %>% filter(geo_code1 == geo_code2))

## -----------------------------------------------------------------------------
(od_min <- od_data_sample[c(1, 2, 9), 1:6])
(od_oneway <- od_oneway(od_min))

## -----------------------------------------------------------------------------
z <- zones_sf
class(z)
l <- od2line(flow = od_inter, zones = z)

## -----------------------------------------------------------------------------
class(l)
nrow(od) - nrow(l)
ncol(l) - ncol(od)

## -----------------------------------------------------------------------------
plot(l$geometry)

## -----------------------------------------------------------------------------
plot(l)

## ---- eval=has_webshot--------------------------------------------------------
library(leaflet)
leaflet() %>%
  addTiles() %>%
  addPolygons(data = l)

## ---- error=TRUE--------------------------------------------------------------
od$geo_code2[3] <- "nomatch"
od2line(od, z)

## ---- eval=FALSE--------------------------------------------------------------
#  library(dplyr)
#  
#  # get nationwide OD data
#  od_all <- pct::get_od()
#  nrow(od_all)
#  # > 2402201
#  od_all$Active <- (od_all$bicycle + od_all$foot) /
#    od_all$all * 100
#  centroids_all <- pct::get_centroids_ew() %>% sf::st_transform(4326)
#  nrow(centroids_all)
#  # > 7201
#  london <- pct::pct_regions %>% filter(region_name == "london")
#  centroids_london <- centroids_all[london, ]
#  od_london <- od_all %>%
#    filter(geo_code1 %in% centroids_london$msoa11cd) %>%
#    filter(geo_code2 %in% centroids_london$msoa11cd)
#  od_london <- od_all[
#    od_all$geo_code1 %in% centroids_london$msoa11cd &
#      od_all$geo_code2 %in% centroids_london$msoa11cd,
#  ]

## ---- eval=FALSE, echo=FALSE--------------------------------------------------
#  # aim: create a reproducible OD dataset
#  od_lnd <- od_london %>%
#    select(-matches("rail|name|moto|car|tax|home")) %>%
#    filter(geo_code2 == "E02000001") %>%
#    top_n(4, wt = all)
#  z_lnd <- centroids_london %>%
#    filter(msoa11cd %in% c(od$geo_code1, od$geo_code2))

## ---- eval=FALSE--------------------------------------------------------------
#  desire_lines_london <- od2line(od_london, centroids_london)
#  nrow(desire_lines_london)
#  # > 352654

## ---- eval=FALSE--------------------------------------------------------------
#  min_trips_threshold <- 20
#  desire_lines_inter <- desire_lines_london %>% filter(geo_code1 != geo_code2)
#  desire_lines_intra <- desire_lines_london %>% filter(geo_code1 == geo_code2)
#  desire_lines_top <- desire_lines_inter %>% filter(all >= min_trips_threshold)
#  nrow(desire_lines_top)
#  # > 28879

## ---- eval=FALSE--------------------------------------------------------------
#  nrow(desire_lines_top) / nrow(desire_lines_london)
#  # > 0.08189046
#  sum(desire_lines_top$all) / sum(desire_lines_london$all)
#  # > 0.557343

## ---- eval=FALSE--------------------------------------------------------------
#  plot(desire_lines_top["all"])

## ---- echo=FALSE--------------------------------------------------------------
knitr::include_graphics("https://user-images.githubusercontent.com/1825120/61058906-030a5c80-a3f0-11e9-90b5-d216964e9681.png")

## ---- eval=FALSE--------------------------------------------------------------
#  lwd <- desire_lines_top$all / mean(desire_lines_top$all) / 10
#  desire_lines_top$percent_dont_drive <- 100 - desire_lines_top$car_driver / desire_lines_top$all * 100
#  plot(desire_lines_top["percent_dont_drive"], lwd = lwd, breaks = c(0, 50, 70, 80, 90, 95, 100))

## ---- echo=FALSE--------------------------------------------------------------
knitr::include_graphics("https://user-images.githubusercontent.com/1825120/62073083-e5ceee00-b237-11e9-9cc7-8bf62d0e9b3f.png")

## ---- eval=FALSE--------------------------------------------------------------
#  library(tmap)
#  desire_lines_top <- desire_lines_top %>%
#    arrange(Active)
#  tm_shape(london) + tm_borders() +
#    tm_shape(desire_lines_top) +
#    tm_lines(
#      palette = "plasma", breaks = c(0, 5, 10, 20, 40, 100),
#      lwd = "all",
#      scale = 9,
#      title.lwd = "Number of trips",
#      alpha = 0.5,
#      col = "Active",
#      title = "Active travel (%)",
#      legend.lwd.show = FALSE
#    ) +
#    tm_scale_bar() +
#    tm_layout(
#      legend.bg.alpha = 0.5,
#      legend.bg.color = "white"
#    )

## ---- echo=FALSE--------------------------------------------------------------
# tmap_save(.Last.value, "tmap-london.png")
knitr::include_graphics("https://user-images.githubusercontent.com/1825120/61066243-12dc6d80-a3fd-11e9-8805-826a47c553f6.png")

## ---- eval=FALSE, echo=FALSE--------------------------------------------------
#  saveRDS(od_all, "od_all.Rds")
#  piggyback::pb_upload("od_all.Rds")

## ---- eval=FALSE--------------------------------------------------------------
#  zones_london <- pct::get_pct_zones("london") %>%
#    select("geo_code")
#  origin_attributes <- desire_lines_top %>%
#    sf::st_drop_geometry() %>%
#    group_by(geo_code1) %>%
#    summarize_if(is.numeric, sum) %>%
#    dplyr::rename(geo_code = geo_code1)
#  # origin_attributes <-
#  zones_origins <- left_join(zones_london, origin_attributes, by = "geo_code")
#  plot(zones_origins, border = NA)

## ---- echo=FALSE--------------------------------------------------------------
knitr::include_graphics("https://user-images.githubusercontent.com/1825120/61067619-e7a74d80-a3ff-11e9-8c15-7467717b36ec.png")

## ---- eval=FALSE--------------------------------------------------------------
#  destination_attributes <- desire_lines_top %>%
#    sf::st_drop_geometry() %>%
#    group_by(geo_code2) %>%
#    summarize_if(is.numeric, sum) %>%
#    dplyr::rename(geo_code = geo_code2) %>%
#    mutate_at(vars(-matches("geo_|all")), funs(. / all)) %>%
#    left_join(zones_london, ., by = "geo_code")
#  
#  plot(destination_attributes, border = NA)

## ---- echo=FALSE--------------------------------------------------------------
knitr::include_graphics("https://user-images.githubusercontent.com/1825120/61069409-27703400-a404-11e9-9c83-1cd5f2397260.png")

## ---- out.width="100%", warning=FALSE, eval=FALSE, echo=FALSE-----------------
#  u <- "https://github.com/ropensci/stplanr/releases/download/0.2.9/lines_cars.Rds"
#  f <- file.path(tempdir(), "lines_cars.Rds")
#  download.file(u, f)
#  lines_cars <- readRDS(f)
#  plot(lines_cars["car_km"], lwd = lines_cars$car_km / 1000)

## ---- eval=FALSE, echo=FALSE--------------------------------------------------
#  sum(lines_cars$car_km * 2.5 * 200) / 1e9

## ---- echo=FALSE, eval=FALSE--------------------------------------------------
#  # out-takes and test code
#  u <- "https://www.rita.dot.gov/bts/sites/rita.dot.gov.bts/files/publications/commodity_flow_survey/2007/zip/origin_destination_files.zip"
#  f <- file.path(tempdir(), "origin_destination_files.zip")
#  download.file(u, f)
#  unzip(f)
#  
#  # demonstrate bug/feature in sf
#  library(sf)
#  m <- matrix(c(
#    0, 0,
#    1, 0,
#    0, 1,
#    0, 0
#  ), ncol = 2)
#  p <- st_polygon(list(m))
#  
#  m <- matrix(c(
#    0, 0,
#    1, 0,
#    0, NA,
#    0, 0
#  ), ncol = 2)
#  p <- st_polygon(list(m))
#  plot(p)
#  
#  l <- st_linestring(m)
#  plot(l)
#  plot(p)
#  m <- matrix(c(0, 0, 0, NA), ncol = 2)
#  l <- st_linestring(m)
#  plot(l)

## ---- echo=FALSE, eval=FALSE--------------------------------------------------
#  usethis::use_data(od_data_sample)
#  # aim: get top flows by car use multiplied by distance
#  # subset flows with more than n people driving:
#  od_cars <- od_data_all[od_data_all$car_driver >= 50, ]
#  cents_ew <- pct::get_centroids_ew()
#  od_cars <- od_cars[
#    od_cars$geo_code1 %in% cents_ew$msoa11cd &
#      od_cars$geo_code2 %in% cents_ew$msoa11cd,
#  ]
#  desire_lines_cars <- od2line(od_cars, cents_ew)
#  plot(desire_lines_cars[1:999, ])
#  desire_lines_cars$euclidean_distance_m <- as.numeric(sf::st_length(desire_lines_cars)) / 1000
#  desire_lines_cars$car_km <- desire_lines_cars$car_driver * desire_lines_cars$euclidean_distance_m
#  lines_cars <- dplyr::top_n(desire_lines_cars, 20000, car_km)
#  summary(lines_cars$car_driver)
#  plot(lines_cars["car_km"])
#  saveRDS(lines_cars, "lines_cars.Rds")
#  piggyback::pb_upload("lines_cars.Rds")

