
#' Persistence Image
#' 
#' Compute the Persistence Image for a given diagram, using piecewise linear weight functions and Gaussian baseline distribution.
#'
#' @param d1 A persistence diagram, in the form of a matrix with 3 columns (first one is the dimension, second is the birth-time, last one is the death-time).
#' @param nbins Number of bins for the discretization of the Persistence Surface into the Persistence Image.
#' @param dimension Dimension of the topological features of interest (0 for connected components, 1 for cycles etc).
#' @param h Standard deviation of the Gaussian baseline used to compute the Persistence Surface.
#' @return a \code{nbins} x \code{nbins} matrix containing the Persistence Image.
#' @author Tullia Padellini
#' @references 
#' \insertRef{adams2017persistence}{kernelTDA}
#' @examples
#' diag1 <- matrix(c(1,1,1,0,2,3,2,2.5,4), ncol = 3, byrow = FALSE)
#' colnames(diag1) <- c("dimension", "birth", "death")
#' pi1 <- pers.image(d1 = diag1, nbins = 20, dimension = 1, h = 1)
#' image(pi1) 
#' @export
pers.image <- function(d1, nbins, dimension, h){
  
  d1 = d1[d1[,1]==dimension, 2:3, drop = F]
  d1[,2] = d1[,2] - d1[,1]

  maxD = max(d1)
  maxP = max(d1[,2])
  
  minD = min(0, min(d1))
  dx = maxD / nbins
  
  
  x_lower = seq(minD, maxD, length.out = nbins)
  x_upper = x_lower + dx
  
  y_lower = seq(0, maxD, length.out = nbins)
  y_upper = y_lower + dx
  
  PSurface = function(point, maxP){
    
    x = point[1]
    y = point[2]
    
    out1 = pnorm(x_upper, mean = x, sd = h) - pnorm(x_lower, mean = x, sd = h)
    out2 = pnorm(y_upper, mean = y, sd = h) - pnorm(y_lower, mean = y, sd = h)
    
    wgt = y / maxP * (y < maxP) + 1 * (y>= maxP)
    return(out1 %o% out2 * wgt)
    
  }

  
  Psurf_mat = apply(d1,1, PSurface, maxP= maxP)
  out = apply(Psurf_mat, 1, sum)

  return(matrix(out, nrow = nbins))

  
}
