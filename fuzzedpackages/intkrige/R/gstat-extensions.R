#=============================================================================
# This file contains wrappers to packages on which intkrige depends:
# gstat, sp, and raster. These functions simplify the variogram calculation
# and spatial visualizations of spatial data.
#=============================================================================
#' Function to create a variogram object for interval-valued data.
#'
#' @param x An object of class \code{intvariogram}.
#' @param models an object of class \code{variogramModelList}. The user must
#'   specify at least two variogram models to fit (for center and radius).
#'   If less than three models are specified then the method does not fit a
#'   variogram for the center/radius interaction.
#' @param ... Additional arguments for \code{\link[gstat]{fit.variogram}}.
#'
#' @return A list of variograms objects from the gstat package.
#'
#' @export
fit.intvariogram <-
  function(x, models = gstat::vgm(rep("Sph", 3)), ...){
    if(class(x)[1] != "intvariogram"){
      stop("Function only defined for objects of class intvariogram")
    }
    if(class(models)[1] != "variogramModelList" && class(models)[1] != "list"){
      stop("models must be specified as a list object")
    }
    if(length(models) < 2 || length(models) > 3){
      stop("between 2 and three variogram models must be specified")
    }
    if(class(models[[1]])[1] != "variogramModel" ||
       class(models[[2]])[1] != "variogramModel"){
      stop("each component of the model argument must be a variogram model")
    }

    # Fit the variogram models in the following order:
    # 1. center, 2. radius, 3. center/radius.
    vorder <- c("center", "radius", "center.radius")
    varioFit <- vector("list", length(models))
    for(i in 1:length(models)){
      x_sub <- x[x$id == vorder[i], ]

      varioFit[[i]] <- gstat::fit.variogram(x_sub,
                                            model = models[[i]], ...)
    }

    return(varioFit)
  }

#' Function to visualize the three variograms from an interval valued
#' spatial data frame.
#'
#' @param x An object of class \code{intvariogram}.
#' @param models A list of fitted variogram models, typically an output of
#' \code{\link{fit.intvariogram}}.
#' @param ... Additional arguments to
#'   \code{\link[gstat]{plot.gstatVariogram}}.
#' @export
intvCheck <- function(x, models, ...){
  p1 <- plot(x[x$id == "center", ], model = models[[1]],
             main = "center", ...)
  p2 <- plot(x[x$id == "radius", ], model = models[[2]],
             main = "radius", ...)
  if(length(models) > 2){
    p3 <- plot(x[x$id == "center.radius", ], model = models[[3]],
               main = "center.radius", ...)
  }else{
    p3 <- plot(x[x$id == "center.radius", ])
  }

  return(gridExtra::grid.arrange(p1, p2, p3, nrow = 2, ncol = 2))
}




