## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

library(landscapetools)

## ----fig.retina=2, message=FALSE, warning=FALSE--------------------------
# Plot continous landscapes
show_landscape(gradient_landscape)

# Plot continous landscapes 
show_landscape(classified_landscape, discrete = TRUE)

# RasterStack/RasterBrick
show_landscape(raster::stack(gradient_landscape, random_landscape), discrete = TRUE)

# Plot a list of raster (list names become facet text)
show_landscape(list("Gradient landscape" = gradient_landscape,
                    "Random landscape" = random_landscape))

# Plot multiple raster with unique scales
show_landscape(raster::stack(gradient_landscape, random_landscape, classified_landscape), unique_scales = TRUE)

## ----fig.retina=2--------------------------------------------------------
# Binarize the landscape into habitat and matrix
binarized_raster <- util_binarize(fractal_landscape, breaks = 0.31415)
show_landscape(binarized_raster)

# You can also provide a vector with thresholds and get a RasterStack with multiple binarized maps
binarized_raster <- util_binarize(fractal_landscape, breaks = c(0.25, 0.5, 0.7))
show_landscape(binarized_raster)

## ----fig.retina=2--------------------------------------------------------
# Mode 1: Classify landscape into 3 classes based on the Fisher-Jenks algorithm:
mode_1 <- util_classify(fractal_landscape, n = 3)

# Mode 2: Classify landscapes into landscape with exact proportions:
mode_2 <- util_classify(fractal_landscape, weighting = c(0.5, 0.25, 0.25))

# Mode 3: Classify landscapes based on a real dataset (which we first create here)
#         and the distribution of values in this real dataset
mode_3 <- util_classify(gradient_landscape, n = 3)

## Mode 3a: ... now we just have to provide the "real landscape" (mode_3)
mode_3a <- util_classify(fractal_landscape, real_land = mode_3)

## Mode 3b: ... and we can also say that certain values are not important for our classification:
mode_3b <- util_classify(fractal_landscape, real_land = mode_3, mask_val = 1)

landscapes <- list(
'Mode 1'  = mode_1,
'Mode 2'  = mode_2,
'Mode 3'  = mode_3,
'Mode 3a' = mode_3a,
'Mode 3b' = mode_3b
)

show_landscape(landscapes, unique_scales = TRUE, nrow = 1)

# ... you can also name the classes:
classified_raster <- util_classify(fractal_landscape,
                                   n = 3,
                                   level_names = c("Land Use 1",
                                                   "Land Use 2",
                                                   "Land Use 3"))
show_landscape(classified_raster, discrete = TRUE)

## ----fig.retina=2--------------------------------------------------------
library(raster) 
landscape <- raster(matrix(1:100, 10, 10))
summary(landscape)

scaled_landscape <- util_rescale(landscape)
summary(scaled_landscape)

## ----fig.retina=2--------------------------------------------------------
# Merge all maps into one
merg <- util_merge(fractal_landscape, c(gradient_landscape, random_landscape), scalingfactor = 1)

# Plot an overview
merge_vis <- list(
    "1) Primary" = fractal_landscape,
    "2) Secondary 1" = gradient_landscape,
    "3) Secondary 2" = random_landscape,
    "4) Result" = merg
)
show_landscape(merge_vis)

## ----eval=FALSE----------------------------------------------------------
#  util_rescale(fractal_landscape, "fractal.asc")

