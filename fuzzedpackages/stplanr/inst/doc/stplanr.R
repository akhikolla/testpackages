## ---- include=FALSE-----------------------------------------------------------
library(stplanr)
build <- curl::has_internet()
knitr::opts_chunk$set(eval = build)

## ---- eval=FALSE--------------------------------------------------------------
#  install.packages("stplanr")

## ---- eval=FALSE--------------------------------------------------------------
#  remotes::install_github("ropensci/stplanr")

## -----------------------------------------------------------------------------
library(stplanr)

## -----------------------------------------------------------------------------
data(package = "stplanr")$result[, "Item"]

## -----------------------------------------------------------------------------
od_eg <- read.csv(
  text =
  "origin, destination, V1, V2
  1, 2, 100, 3
  1, 3, 50, 5"
)
knitr::kable(od_eg)

## -----------------------------------------------------------------------------
head(flow[c(1:3, 12)])

## -----------------------------------------------------------------------------
head(cents_sf)

## -----------------------------------------------------------------------------
library(sf)
class(cents_sf)
plot(cents_sf)

## -----------------------------------------------------------------------------
flow_single_line <- od_data_sample[2:3, ] # select only the first line
desire_line_single <- od2line(flow = flow_single_line, zones = cents_sf)

## -----------------------------------------------------------------------------
plot(desire_line_single$geometry, lwd = 5)
plot(cents_sf, add = TRUE, cex = 5)

## -----------------------------------------------------------------------------
l <- od2line(flow = flow, zones = cents_sf)
# identify 'intrazone flows'
sel_intra <- l$Area.of.residence == l$Area.of.workplace
# find distances
l_distances <- geo_length(l)
summary(l_distances)
sel_dist <- l_distances > 2000
sel <- !sel_intra & sel_dist
l <- l[sel, ]

## ---- eval=FALSE--------------------------------------------------------------
#  plot(l)

## ---- echo=FALSE--------------------------------------------------------------
l_bb <- sf::st_bbox(l)
# l_bb[1] <- NA
no_na_in_bb <- !any(is.na(as.numeric(l_bb)))

## ---- eval=no_na_in_bb--------------------------------------------------------
lwd <- l$All / mean(l$All)
plot(st_geometry(l), lwd = lwd)

## ---- eval=no_na_in_bb--------------------------------------------------------
plot(l["Bicycle"], lwd = lwd)

## ---- message=FALSE, warning=FALSE--------------------------------------------
# if the next line returns FALSE the code will not run
(has_internet <- curl::has_internet())
(cs_key <- nchar(Sys.getenv("CYCLESTREETS")))
if (has_internet & cs_key == 16) {
  r <- route(l = l, route_fun = cyclestreets::journey)
  r <- aggregate(r[c(3, 12)], by = list(r[[1]], r[[2]]), FUN = mean)
} else {
  r <- routes_fast_sf[sel, ]
}

## ---- out.width="500", out.height="500", eval=FALSE---------------------------
#  plot(r$geometry, lwd = lwd * 3, reset = FALSE)

## ---- out.width="500", out.height="500", echo=FALSE, eval=FALSE---------------
#  # alternative showing buildings:
#  r_sf <- st_sf(l, geometry = st_as_sfc(r))
#  if (require(osmdata)) {
#    buildings <- opq(st_bbox(l)) %>%
#      add_osm_feature(key = "building", value = "industrial") %>%
#      osmdata_sf()
#  }
#  plot(r_sf["Bicycle"], lwd = lwd * 3, reset = FALSE)
#  plot(st_geometry(buildings$osm_polygons), col = "grey", add = TRUE)

