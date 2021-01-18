## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ---- message=FALSE-----------------------------------------------------------
library(SpatialKDE)
library(sp)
library(sf)
library(dplyr)
library(tmap)

## -----------------------------------------------------------------------------
data(meuse)

meuse <- meuse %>% 
  st_as_sf(coords = c("x", "y"), dim = "XY") %>% 
  st_set_crs(28992) %>% 
  select()

## -----------------------------------------------------------------------------
cell_size <- 100
band_width <- 150

## -----------------------------------------------------------------------------
grid_meuse <- meuse %>% 
  create_grid_rectangular(cell_size = cell_size, side_offset = band_width)

## -----------------------------------------------------------------------------
kde <- meuse %>% 
  kde(band_width = band_width, kernel = "quartic", grid = grid_meuse)

## -----------------------------------------------------------------------------
tm_shape(kde) + 
  tm_polygons(col = "kde_value", palette = "viridis", title = "KDE Estimate") +
tm_shape(meuse) +
  tm_bubbles(size = 0.1, col = "red")

## -----------------------------------------------------------------------------
raster_meuse <- meuse %>% 
  create_raster(cell_size = cell_size, side_offset = band_width)

## -----------------------------------------------------------------------------
kde <- meuse %>% 
  kde(band_width = band_width, kernel = "triweight", grid = raster_meuse)

## -----------------------------------------------------------------------------
tm_shape(kde) + 
  tm_raster(palette = "viridis", title = "KDE Estimate") +
tm_shape(meuse) +
  tm_bubbles(size = 0.1, col = "red") + 
tm_layout(legend.outside = TRUE)

