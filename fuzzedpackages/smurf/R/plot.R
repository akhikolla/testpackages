###############################################
#
# Plot model coefficients
#
###############################################



#' @export
#' @title Plot Coefficients of Estimated Model
#' 
#' @description Function to plot the coefficients of the estimated model. 
#' 
#' @param x An object of class '\code{\link[=glmsmurf-class]{glmsmurf}}', typically the result of a call to \code{\link{glmsmurf}} or \code{\link{glmsmurf.fit}}.
#' @param xlab Label for the x-axis, default is \code{"Index"}.
#' @param ylab Label for the y-axis, default is \code{"Estimated coefficients"}.
#' @param basic Logical indicating if the basic lay-out is used for the plot, default is \code{FALSE}.
#' @param ... Additional arguments for the \code{\link[graphics:plot.default]{plot}} function.
#' 
#' @details When \code{basic=FALSE}, the improved lay-out for the plot is used. Per predictor, groups of equal coefficients are indicated
#' in the same color (up to 8 colors), and zero coefficients are indicated by grey squares. 
#'  
#' @seealso \code{\link{plot_reest}}, \code{\link{coef.glmsmurf}}, \code{\link{summary.glmsmurf}}, \code{\link{glmsmurf}}, \code{\link{glmsmurf-class}}
#'
#' @examples ## See example(glmsmurf) for examples
plot.glmsmurf <- function(x, xlab = "Index", ylab = "Estimated coefficients", basic = FALSE, ...) {
  
  .plot.glmsmurf.aux(x = x , xlab = xlab, ylab = ylab, basic = basic, reest = FALSE, ...)
  
}


#' @export
#' @title Plot Coefficients of Re-estimated Model
#' 
#' @description Function to plot the coefficients of the re-estimated model. 
#' 
#' @param x An object for which the extraction of model coefficients is meaningful. 
#'          E.g. an object of class '\code{\link[=glmsmurf-class]{glmsmurf}}', typically the result of a call to \code{\link{glmsmurf}} or \code{\link{glmsmurf.fit}}.
#' @param ylab Label for the y-axis, default is \code{"Re-estimated coefficients"}.
#' @inheritParams plot.glmsmurf
#' 
#' @details When the re-estimated model is not included in \code{x}, 
#'          the coefficients of the estimated model in \code{x} are plotted with a warning.
#'          
#'          See \code{\link{plot.glmsmurf}} for more details.
#'         
#' @seealso \code{\link{plot.glmsmurf}}, \code{\link{coef_reest}}, \code{\link{summary.glmsmurf}}, \code{\link{glmsmurf}}, \code{\link{glmsmurf-class}}
#' 
#' @examples ## See example(glmsmurf) for examples
plot_reest <- function(x, ...) UseMethod("plot_reest", x)

#' @export
#' @rdname plot_reest
plot_reest.glmsmurf <- function(x, xlab = "Index", ylab = "Re-estimated coefficients", basic = FALSE, ...) {
  
    .plot.glmsmurf.aux(x = x , xlab = xlab, ylab = ylab, basic = basic, reest = TRUE, ...)
}

# Auxiliary function to plot coefficients of estimated model (reest=FALSE) or re-estimated model (reest=TRUE)
.plot.glmsmurf.aux <- function(x, xlab, ylab, reest, basic, ...) {
  
  # Handle reest
  if (!exists("coefficients.reest", x) & reest) {
    warning("Coefficients of the re-estimated model are not present in 'x', the coefficients of the estimated model are plotted.")
    reest <- FALSE
    
    if (ylab == "Re-estimated coefficients") {
      ylab <- "Estimated coefficients"
    }
  }
  
  # Select coefficients to use
  if (reest) {
    # Use re-estimated coefficients
    coefs <- coef_reest(x)
    
  } else {
    # Use estimated coefficients
    coefs <- coef(x)
  }

  if (basic) {
    
    # New plot  
    plot(coefs, xlab = xlab, ylab = ylab, ...)
    
  } else {
    
    # Number of estimated parameters per covariate
    n.par.cov <- x$n.par.cov
    # Penalty type per covariate
    pen.cov <- x$pen.cov
    # Number of covariates
    n.cov <- length(n.par.cov)
    
    # Get color palette
    col.palette <- brewer.pal(8, "Dark2")
    
    # New empty plot  
    plot(coefs, xlab = xlab, ylab = ylab, type = "n", ...)
    
    # Add horizontal line at 0
    abline(h = 0, lty = 2)
    
    
    # Point types, default is filled circle
    pchs <- rep(16, length(coefs))
    # Colors, default is black
    cols <- rep("black", length(coefs))
    
    # Set point type for zero coefficients (filled square with edges in different color)
    pchs[coefs == 0] <- 22
    
    for (j in 1:n.cov) {
      
      # Index of first coefficient corresponding to this predictor
      ind.start <- ifelse(j == 1, 1L, sum(unlist(n.par.cov[1:(j-1L)])) + 1L)
      # Index of last coefficient corresponding to this predictor
      ind.end <- sum(unlist(n.par.cov[1:j]))
      
      # Get unique non-zero coefficients for this covariate
      unique.coefs <- unique(coefs[ind.start:ind.end])
      unique.coefs <- unique.coefs[unique.coefs != 0]
      
      # Color index
      l2 <- 1L
      
      if (pen.cov[[j]] %in% c("flasso", "gflasso", "ggflasso", "2dflasso")) {
        
        for (l in 1:length(unique.coefs)) {
          
          # Get indices of coefficients of this predictor that are equal to selected coefficient
          ind <- which(coefs == unique.coefs[l])
          ind <- ind[ind >= ind.start & ind <= ind.end]
          
          # Color per bin
          if (pen.cov[[j]] == "flasso") {
            
            cols[ind] <- col.palette[l2]
            
            if (l2 >= length(col.palette)) {
              # Reset color index
              l2 <- 0L
            }
            
          } else {
            # Do not reset color index for these penalty types, just return black dot instead
            cols[ind] <- ifelse(l2 <= length(col.palette), col.palette[l2], "black")
          }
          
          l2 <- l2 + 1L
        }
      }
      
      # Issue warning if too many bins. Can only happen for "gflasso", "ggflasso", "2dflasso" penalty types
      if (l2 > length(col.palette) + 1) {
        warning(paste0("There are more than ", length(col.palette), " coefficient bins for the predictor '", 
                       names(n.par.cov)[j],"'. ",
                       "Only the first ", length(col.palette), " coefficient bins are indicated in color."))
      }
      
    }
    
    # Add points, bg is needed for zero coefficients
    points(coefs, col = cols, pch = pchs, bg = 'grey', ...)
    
    # Add vertical lines between predictors (if more than one predictor)
    if (n.cov > 1) {
      cs <- cumsum(unlist(n.par.cov))
      abline(v = cs[-length(cs)] + 0.5, lty = 3)
    }
  }
  
} 
