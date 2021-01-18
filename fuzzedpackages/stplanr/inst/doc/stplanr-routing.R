## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(stplanr)

## ---- eval=FALSE--------------------------------------------------------------
#  remotes::install_github("itsleeds/geofabrik")
#  library(geofabrik)
#  roads_iow = get_geofabrik(name = "Isle of Wight")
#  f = gf_filename("Isle of Wight")
#  file.copy(f, "iow.pbf")
#  options(osrm.server = "https://0.0.0.0:5000/", osrm.profile = "driving")

## -----------------------------------------------------------------------------
library(osrm)

## ---- eval=FALSE--------------------------------------------------------------
#  l = pct::wight_lines_30
#  p = line2points(l)
#  r = osrm::osrmRoute(src = p[1, ], dst = p[2, ], returnclass = "sf", overview = "full")
#  plot(r)

## ---- echo=FALSE--------------------------------------------------------------
knitr::include_graphics("https://user-images.githubusercontent.com/1825120/86902789-577d1080-c106-11ea-91df-8d0180931562.png")

## ---- eval=FALSE--------------------------------------------------------------
#  routes_osrm = route(l = l, route_fun = osrmRoute, returnclass = "sf", overview = "full")
#  rnet_osrm = overline(routes_osrm, attrib = "bicycle")
#  mapview::mapview(rnet_osrm, lwd = rnet_osrm$bicycle / 10)

## ---- eval=FALSE, echo=FALSE--------------------------------------------------
#  system.time({
#    routes_osrm = route(l = l, route_fun = osrmRoute, returnclass = "sf", overview = "full")
#  })
#  30 / 0.9 # around 30 routes per second
#  saveRDS(routes_osrm, "routes_osrm.Rds")
#  piggyback::pb_upload("routes_osrm.Rds")

## ---- echo=FALSE--------------------------------------------------------------
knitr::include_graphics("https://user-images.githubusercontent.com/1825120/86858225-2970df80-c0b8-11ea-8394-07f98f1c8e8a.png")

## ---- eval=FALSE--------------------------------------------------------------
#  # tidy up
#  f = list.files(pattern = "iow")
#  unlink(x = f, recursive = TRUE)

