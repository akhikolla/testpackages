## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(intkrige)
data(ohtemp)
head(ohtemp)

## ----ohioMap, eval = FALSE, echo = FALSE--------------------------------------
#  
#  library(ggmap)
#  # # Load the Ohio River Basin Shapefile
#  # ohMap <- rgdal::readOGR(dsn = "ohMap", layer = "ohMap")
#  #
#  # # Simplify the shapefile
#  # ohMap <- rgeos::gSimplify(ohMap, tol = .05, topologyPreserve = TRUE)
#  data(ohMap)
#  
#  # key has been removed for security reasons
#  # register_google(key = "insert key here")
#  register_google(key = "<my key>")
#  center = apply(ohMap@bbox, 1, mean)
#  
#  oh.google <- ggmap::get_googlemap(center = center, zoom = 6,
#                                      maptype = "terrain",
#                                      color = "bw",
#                                    style = "feature:all|element:labels|visibility:off|
#                             &style=feature:road|element:all|visibility:off|
#                             &style=feature:landscape.man_made|element:all|visibility:off|")
#  
#  map <- ggmap(oh.google) +
#    xlim(ohMap@bbox[1, 1], ohMap@bbox[1, 2]) +
#    xlab("Longitude (easting)") +
#    ylim(ohMap@bbox[2, 1], ohMap@bbox[2, 2]) +
#    ylab("Latitude (northing)") +
#    geom_polygon(data = ohMap, aes(x = long, y = lat),
#                 colour = "gray50", fill = NA, lwd = 1.1, lty = 1) +
#    geom_point(data = as.data.frame(ohtemp), aes(x = LONGITUDE, y = LATITUDE),
#               pch = 16, alpha = 0.9) +
#    theme(legend.title = element_blank(),
#          legend.text = element_text(size = 16),
#          axis.text = element_text(size = 16),
#          axis.title = element_text(size = 16))
#  
#  pdf("vignettes/ohmap.png", width = 4, height = 4)
#  map
#  dev.off()

## -----------------------------------------------------------------------------
# First, create a SpatialPointsDataFrame in the usual way
sp::coordinates(ohtemp) <- c("LONGITUDE", "LATITUDE")
sp::proj4string(ohtemp) <- CRS("+proj=longlat +ellps=WGS84")
interval(ohtemp) <- c("minm", "maxm")

head(ohtemp)

## -----------------------------------------------------------------------------
interval(ohtemp) <- log(interval(ohtemp))
head(ohtemp)

## ---- fig.align='center'------------------------------------------------------
# Revert back to the standard interval
interval(ohtemp) <- exp(interval(ohtemp))
varios <- intvariogram(ohtemp, cutoff = 500)

plot(varios)

## ---- fig.align='center'------------------------------------------------------
varioFit <- fit.intvariogram(varios, models = gstat::vgm(c("Lin", "Sph", "Sph")))
varioFit
intvCheck(varios, varioFit)

## ---- fig.align='center'------------------------------------------------------
# Replace non-convergent variogram fit with a surrogate that 
# contains a reasonable range and sill. 
varioFit[[1]] <- gstat::vgm(psill = 350, nugget = 4.608542, 
                            range = 600, model = "Sph")

intvCheck(varios, varioFit)

## ---- fig.align = 'center'----------------------------------------------------
# Include the Ohio river basin shapefile
data(ohMap)

# New location data preparation
lon <- seq(-89.26637, -77.83937, length.out = 10)
lat <- seq(35.31332, 42.44983, length.out = 10)
newlocations <- expand.grid(lon, lat)
colnames(newlocations) <- c("lon", "lat")
sp::coordinates(newlocations) <- c("lon", "lat")
sp::proj4string(newlocations) <- sp::proj4string(ohtemp)
sp::gridded(newlocations) <- TRUE

# Adjust r and theta to ensure answers remain in feasible region.
preds <- intkrige(ohtemp, newlocations, varioFit, 
                  A = c(1, 1, 0.5), r = 200, eta = 0.9, maxp = 225)

plot(preds, beside = FALSE, circleCol = "gray") + 
  latticeExtra::layer(sp.lines(ohMap, col = "white"))
plot(preds, beside = TRUE)

