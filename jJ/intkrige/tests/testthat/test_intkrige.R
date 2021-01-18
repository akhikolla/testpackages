context("Intkrige predictions")
library(intkrige)

test_that("prediction for sample case is as expected (c++ version)", {
  # First, define the location and elevation of interest.
  # (In this case we pick coordinates of Utah State University)
  templocs <- data.frame(lat = 41.745, long = -111.810, ELEVATION = 1456)
  sp::coordinates(templocs) <- c("long", "lat")
  sp::proj4string(templocs) <- "+proj=longlat +ellps=WGS84
  +datum=WGS84 +no_defs +towgs84=0,0,0"

  # Load the Utah Snow Load Data
  data(utsnow)
  utsnow.sp <- utsnow

  # Convert to an 'intsp' object that inherits a SpatialPointsDataFrame
  sp::coordinates(utsnow.sp) <- c("LONGITUDE", "LATITUDE")
  sp::proj4string(utsnow.sp) <- sp::proj4string(templocs)
  interval(utsnow.sp) <- c("minDL", "maxDL")

  # analyze interval on log scale
  interval(utsnow.sp) <- log(interval(utsnow.sp))

  # Define the formulas we will use to define the intervals.
  temp_formulas <- list(center ~ ELEVATION,
                        radius*(ELEVATION/median(ELEVATION)) ~ 1)

  # Define, fit and check the variogram fits.
  varios <- intvariogram(utsnow.sp,
                         formulas = temp_formulas)
  varioFit <- fit.intvariogram(varios, models = gstat::vgm(c("Sph", "Sph", "Gau")))

  preds <- intkrige::intkrige(locations = utsnow.sp,
                              newdata = templocs,
                              models = varioFit,
                              formulas = temp_formulas)
  preds2 <- intkrige::intkrige(locations = utsnow.sp,
                              newdata = templocs,
                              models = varioFit,
                              formulas = temp_formulas, useR = FALSE)


  target <- data.frame(lower = -0.0361, upper = 0.9703)
  # The final results are predicted intervals after removing the effect of elevation.
  expect_equal(round(interval(preds), 4), as.matrix(target))
  expect_equal(round(interval(preds2), 4), as.matrix(target))
})
