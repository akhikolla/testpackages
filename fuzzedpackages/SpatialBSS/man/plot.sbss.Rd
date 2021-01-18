\name{plot.sbss}
\alias{plot.sbss}

\title{
Plot Method for an Object of Class 'sbss'
}

\description{
\code{plot.sbss} is an interface to the standard plot method for the class of the estimated source random field. 
}

\usage{
\method{plot}{sbss}(x, which = 1:ncol(x$s), \dots)
}

\arguments{
  \item{x}{
object of class \code{'sbss'}. Usually result of \code{\link{sbss}}.
}
  \item{which}{
a numeric vector indicating which components of the latent field should be plotted.
}
  \item{\dots}{
further arguments to the plot method of \code{class(x$s)}, which is either \code{\link[sp]{spplot}} or \code{\link[sf]{plot}}.
}
}

\details{
This method calls the corresponding plot method of \code{class(x$s)}. Either \code{\link[sp]{spplot}} for \code{class(x$s)} is \code{\link[sp]{SpatialPointsDataFrame}} or \code{\link[sf]{plot.sf}} for \code{class(x$s)} is \code{\link[sf]{sf}}. If \code{x$s} is a matrix then it is internally cast to \code{\link[sp]{SpatialPointsDataFrame}} and \code{\link[sp]{spplot}} is used for plotting. Arguments to the corresponding plot functions can be given through \code{\dots}.
}

\seealso{
\code{\link{sbss}}, \code{\link[sp]{spplot}}, \code{\link[sf]{plot.sf}}
}

\examples{
# simulate coordinates
coords <- runif(1000 * 2) * 20
dim(coords) <- c(1000, 2)

# simulate random field
if (!requireNamespace('RandomFields', quietly = TRUE)) {
  stop('Please install the package RandomFields to run the example code.')
} else {
  RandomFields::RFoptions(spConform = FALSE)
  field_1 <- RandomFields::RFsimulate(model = RandomFields::RMexp(), 
                                      x = coords)
  field_2 <- RandomFields::RFsimulate(model = RandomFields::RMspheric(), 
                                      x = coords)
  field_3 <- RandomFields::RFsimulate(model = RandomFields::RMwhittle(nu = 2),  
                                      x = coords)
  field <- cbind(field_1, field_2, field_3)
}

# compute ring kernel matrices
kernel_parameters <- c(0, 1, 1, 2, 2, 3)
ring_kernel_list <- spatial_kernel_matrix(coords, 'ring', kernel_parameters)

# apply sbss SpatialPointsDataFrame object
field_sp <- sp::SpatialPointsDataFrame(coords = coords, data = data.frame(field))
res_sp <- sbss(field_sp, kernel_list = ring_kernel_list)

# plot with SpatialPointsDataFrame object
plot(res_sp)

# plot with SpatialPointsDataFrame object 
# and additional arguments for spplot function
plot(res_sp, colorkey = TRUE, as.table = TRUE, cex = 1)

# apply sbss with sf object
if (!requireNamespace('sf', quietly = TRUE)) {
  stop('Please install the package sf to run the example code.')
} else {
  field_sf <- sf::st_as_sf(data.frame(coords = coords, field), 
                           coords = c(1,2))
  res_sf <- sbss(x = field_sf, kernel_list = ring_kernel_list)
}

# plot with sf object
plot(res_sf)

# plot with sf object 
# and additional arguments for plot.sf function
plot(res_sf, axes = TRUE, key.pos = 4)

}
