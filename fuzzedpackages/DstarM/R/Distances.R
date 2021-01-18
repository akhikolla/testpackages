#' Calculates the distance between two probability densities.
#'
#'
#' @param tt the time grid on which the densities are evaluated.
#' @param a a vector with values of the first density.
#' @param b a vector with values of the second density.
#' @return The distance between densities \code{a} and \code{b}.
#'
#' @examples
#' # Lets simulate a bunch of parameters and compare the three distance measures.
#'
#' tt = seq(0, 5, .001)
#' parsMatV = cbind(.8, seq(0, 5, .5), .5, .5, .5) # differ only in drift speed
#' parsMatA = cbind(seq(.5, 2, .15), 2, .5, .5, .5)# differ only in boundary
#' # calculate densities for all these parameters
#' dV = apply(parsMatV, 1, function(x, tt) Voss.density(tt, x, boundary = 'upper'), tt = tt)
#' dA = apply(parsMatA, 1, function(x, tt) Voss.density(tt, x, boundary = 'upper'), tt = tt)
#' # make plots of the densities
#' matplot(tt, dA, xlim = c(0, .6), main = 'Densities with different Boundary',
#'         col = rainbow(ncol(dA)),type = 'l', lty = 1, las = 1, bty = 'n',
#'         xlab = 'Time', ylab = 'Density')
#' legend('topright', lty = 1, bty = 'n', col = rainbow(ncol(dA)),
#'        legend = paste('a = ', parsMatA[, 1]))
#' matplot(tt, dV, xlim = c(0, .6), main = 'Densities with different Drift Speed',
#'         col = rainbow(ncol(dV)), type = 'l', lty = 1, las = 1, bty = 'n',
#'         xlab = 'Time', ylab = 'Density')
#' legend('topright', lty = 1, bty = 'n', col = rainbow(ncol(dV)),
#'        legend = paste('v = ',parsMatV[, 2]))
#' # empty matrices for data storage
#' distMatV = matrix(NA, nrow = ncol(dV) - 1, ncol = 3,
#'                   dimnames = list(NULL, c('Chisq', 'Bhattacharyya', 'Hellinger')))
#' distMatA = matrix(NA, nrow = ncol(dA) - 1, ncol = 3,
#'                   dimnames = list(NULL, c('Chisq', 'Bhattacharyya', 'Hellinger')))
#' # calculate distances between densities in column i and i + 1.
#' # this is done using three different distance measures
#' for (i in 1:(ncol(dA) - 1)) {
#'   distMatV[i, ] = c(chisq(tt, dV[, i], dV[, i + 1]),
#'                     battacharyya(tt, dV[, i], dV[, i + 1]),
#'                     hellinger(tt, dV[, i], dV[, i + 1]))
#'   distMatA[i, ] = c(chisq(tt, dA[, i], dA[, i + 1]),
#'                     battacharyya(tt, dA[, i], dA[, i + 1]),
#'                     hellinger(tt, dA[, i], dA[, i + 1]))
#' }
#' # The three distance measures correlate highly for differences in Boundary
#' cor(distMatA)
#' # The battacharyya distance measures does not correlate with the others
#' # when calculating differences in drift speed
#' cor(distMatV)

#' @export
# calc chisq dist via simpsons rule; scalar | simpson
chisq <- function(tt, a, b) {
  # Calc chi-square dist
  vals <- (a - b)^2/(a + b + 1e-10)  # + small float to avoid divide by 0
  return(simpson(tt, vals))
}

#' @rdname chisq
#' @export
# Bhattacharyya distance; scalar | simpson
battacharyya <- function(tt, a, b) {
  a[a < 0L] <- 0L
  b[b < 0L] <- 0L
  return(abs(log(simpson(tt, sqrt(a * b)))))
}

#' @rdname chisq
#' @export
# Hellinger distance; scalar | simpson
hellinger <- function(tt, a, b) {
  # browser()
  a[a < 0L] <- 0L
  b[b < 0L] <- 0L
  # return(abs(1 - simpson(tt, sqrt(a * b))))
  return(0.5 * simpson(tt, (sqrt(a) - sqrt(b))^2))
}

