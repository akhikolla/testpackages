## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

required <- c("viridis")
if (!all(sapply(required, requireNamespace, quietly = TRUE))) {
  knitr::opts_chunk$set(eval = FALSE)
}
options(max.print = 20)

## ---- message = FALSE---------------------------------------------------------
# First step is to load the libraries. Not all of these libraries are stricly
# needed; some are used for convenience and visualization for this tutorial.
library("samc")
library("raster")
library("viridis")


# "Load" the data. In this case we are using data built into the package.
# In practice, users will likely load raster data using the raster() function
# from the raster package.
res_data <- samc::ex_res_data
abs_data <- samc::ex_abs_data
occ_data <- samc::ex_occ_data


# Create a samc object using the resistance and absorption data. We use the 
# recipricol of the arithmetic mean for calculating the transition matrix. Note,
# the input data here are matrices, not RasterLayers. If using RasterLayers, the
# `latlon` parameter must be set.
samc_obj <- samc(res_data, abs_data, tr_fun = function(x) 1/mean(x))


## -----------------------------------------------------------------------------
# First, let's specify the time steps we are interested in
time_steps <- c(10, 100, 1000, 10000)

## ---- fig.width = 6.5, out.width = '100%', fig.align = "center"---------------
# First, calculate the metric for the entire vector of time steps
results <- distribution(samc_obj, origin = 1, time = time_steps)

# The result is a list of vectors. Note that the list is named with the time steps
str(results)

# We can take this list, and use map() to convert it to a list of RasterLayers.
results_map <- map(samc_obj, results)
str(results_map, max.level = 1) # max.level is to hide a lot of gory details

# Conveniently, this is easy to turn into a RasterStack
results_stack <- raster::stack(results_map)

# A caveat: RasterStacks cannot have numeric names, so it prepends an 'X'
names(results_stack)

# RasterStacks are convenient for a lot of different things, like processing the
# rasters all at once using spatial methods. But we're just going to plot them
plot(results_stack, xlab = "x", ylab = "y", col = viridis(256))

## ---- max.print=5, fig.width = 6.5, out.width = '100%', dpi = 300, fig.align = "center"----
# The results of indivudal time steps can be retrieved by either index or by name 
results[[3]]
results[["1000"]]

# The latter is particularly useful in making reliable for loops because indexed
# for loops can be harder to troubleshoot if something goes wrong. 
for (ts in time_steps) {
  name <- as.character(ts)
  r <- results[[name]]
  r_map <- map(samc_obj, r)
  plot(r_map, main = paste("Individual Location at Time", ts), xlab = "x", ylab = "y", col = viridis(256))
}

# For comparison, here is an indexed for loop that does the same thing
for (i in 1:length(time_steps)) {
  name <- as.character(time_steps[i])
  r <- results[[i]] # or we could replace `i` with `name` on this line
  # ...
}

## ---- fig.width = 6.5, out.width = '100%', dpi = 300, fig.align = "center"----
for (ts in time_steps) {
  dist <- distribution(samc_obj, origin = 1, time = ts)
  dist_map <- map(samc_obj, dist)
  plot(dist_map, main = paste("Individual Location at Time", ts), xlab = "x", ylab = "y", col = viridis(256))
}


