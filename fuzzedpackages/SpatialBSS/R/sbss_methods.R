#-------------------------------------------#
# coef
#-------------------------------------------#
coef.sbss <- function(object, ...) {
  object$w
}

#-------------------------------------------#
# print
#-------------------------------------------#
print.sbss <- function(x, ...) {
  print.listof(x[!(names(x) %in% c('s', 'coords', 'w_inv', 'x_mu', 'cov_inv_sqrt'))], ...)
}

#-------------------------------------------#
# predict
#-------------------------------------------#
predict.sbss <- function(object, p = 2, n_grid = 50, which = 1:ncol(object$s), ...) {
  if ('matrix' %in% class(object$s)) {
    if (!is.null(object$coords)) {
      vals_idw <- predict_idw(vals = object$s[, which, drop = FALSE], coords = object$coords, p = p, n_grid = n_grid)
      plot(sp::spplot(sp::SpatialPointsDataFrame(coords = vals_idw$coords_pred_idw, 
                                                 data = data.frame(vals_idw$vals_pred_idw)), ...))
      invisible(vals_idw)
    } else {
      stop('Please call the function sbss with the argument coords given.')
    }    
  } else if ('SpatialPointsDataFrame' %in% class(object$s)) {
    vals_idw <- predict_idw(vals = as.matrix(object$s[, which]@data), coords = object$s@coords, p = p, n_grid = n_grid)
    vals_idw <- sp::SpatialPointsDataFrame(coords = vals_idw$coords_pred_idw, 
                                           data = data.frame(vals_idw$vals_pred_idw))
    plot(sp::spplot(vals_idw, ...))
    invisible(vals_idw)    
  } else if ('sf' %in% class(object$s)) {
    if (!requireNamespace('sf', quietly = TRUE)) {
      stop('Please install the package sf to use this function.')
    } else {
      vals_idw <- predict_idw(vals = as.matrix(sf::st_drop_geometry(object$s[, which])), 
                              coords = sf::st_coordinates(object$s), p = p, n_grid = n_grid)
      vals_idw <- sf::st_as_sf(data.frame(coords = vals_idw$coords_pred_idw, vals_idw$vals_pred_idw), 
                               coords = c(1,2))
      plot(vals_idw, ...)
      invisible(vals_idw)
    }
  } else {
    stop('Unknown class of the latent field.')
  }
}

#-------------------------------------------#
# plot
#-------------------------------------------#
plot.sbss <- function(x, which = 1:ncol(x$s), ...) {
  if ('matrix' %in% class(x$s)) {
    if (!is.null(x$coords)) {
      plot(sp::spplot(sp::SpatialPointsDataFrame(coords = x$coords, data = data.frame(x$s[, which])), ...))
    } else {
      stop('Please call the function sbss with the argument coords given.')
    }    
  } else if ('SpatialPointsDataFrame' %in% class(x$s)) {
    plot(sp::spplot(x$s[, which], ...))
  } else if ('sf' %in% class(x$s)) {
    if (!requireNamespace('sf', quietly = TRUE)) {
      stop('Please install the package sf to use this function.')
    } else {
      plot(x$s[, which], ...)
    }
  } else {
    stop('Unknown class of the latent field.')
  }
}

