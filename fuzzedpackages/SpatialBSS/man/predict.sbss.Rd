\name{predict.sbss}
\alias{predict.sbss}

\title{
Predict Method for an Object of Class 'sbss'
}

\description{
\code{predict.sbss} predicts the estimated source random field on a grid with Inverse Distance Weighting (IDW) and plots these predictions.
}

\usage{
\method{predict}{sbss}(object, p = 2, n_grid = 50, which = 1:ncol(object$s), \dots)
}

\arguments{
  \item{object}{
object of class \code{'sbss'}. Usually result of \code{\link{sbss}}.
}
  \item{p}{
numeric. The positive power parameter for IDW. Default is 2.
}
  \item{n_grid}{
numeric. Each dimension of the spatial domain is divided by this integer to derive a grid for IDW predictions. Default is 50.
}
  \item{which}{
a numeric vector indicating which components of the latent field should be predicted.
}
  \item{\dots}{
further arguments to the plot method of \code{class(x$s)}, which is either \code{\link[sp]{spplot}} or \code{\link[sf]{plot}}.
}
}

\details{
\href{https://en.wikipedia.org/wiki/Inverse_distance_weighting}{IDW} predictions are made on a grid. The side lengths of the rectangular shaped grid cells are derived by the differences of the rounded maximum and minimum values divided by the \code{n_grid} argument for each column of \code{object$coords}. Hence, the grid contains a total of \code{n_grid ^ 2} points. The power parameter of the IDW predictions is given by \code{p} (default: 2).

The predictions are plotted with the corresponding plot method of \code{class(x$s)}. Either \code{\link[sp]{spplot}} for \code{class(x$s)} is \code{\link[sp]{SpatialPointsDataFrame}} or \code{\link[sf]{plot.sf}} for \code{class(x$s)} is \code{\link[sf]{sf}}. If \code{x$s} is a matrix then it is internally cast to \code{\link[sp]{SpatialPointsDataFrame}} and \code{\link[sp]{spplot}} is used for plotting. Arguments to the corresponding plot functions can be given through \code{\dots} as it is done by the method \code{\link{plot.sbss}}.
}

\value{
The return is dependent on the class of the latent field in the \code{'sbss'} object.
If \code{class(object$s)} is a matrix then a list with the following entries is returned:
\item{vals_pred_idw}{a matrix of dimension \code{c(n,p)} (when \code{which} is default or less than \code{p} columns according to the selected components with the \code{which} argument) with the IDW predictions of the estimated source random field.}
\item{coords_pred_idw }{a matrix of dimension \code{c(n,2)} with the grid coordinates for the IDW predictions.}

If \code{class(object$s)} is \code{\link[sp]{SpatialPointsDataFrame}} or \code{\link[sf]{sf}} then the predicted values and their coordinates are returned as an object of the corresponding class.

The return is invisible.

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

# apply sbss with three ring kernels
kernel_borders <- c(0, 1, 1, 2, 2, 4)
res_sbss <-  sbss(field, coords, 'ring', kernel_borders)

# predict latent fields on grid with default settings
predict(res_sbss)

# predict latent fields on grid with custom plotting settings
predict(res_sbss, colorkey = TRUE, as.table = TRUE, cex = 1)

# predict latent fields on a 60x60 grid
predict(res_sbss, n_grid = 60, colorkey = TRUE, as.table = TRUE, cex = 1)

# predict latent fields with a higher IDW power parameter
predict(res_sbss, p = 10, colorkey = TRUE, as.table = TRUE, cex = 1)

# predict latent fields and save the predictions
predict_list <- predict(res_sbss, p = 5, colorkey = TRUE, as.table = TRUE, cex = 1)
}

\seealso{
\code{\link{sbss}}, \code{\link{plot.sbss}}, \code{\link[sp]{spplot}}, \code{\link[sf]{plot.sf}}
}

