#' @title Incomplete data generator
#'
#' @description Generate a matrix with missing values,
#' where the indices of missing values are uniformly randomly distributed in the matrix.
#'
#' @param m the rows of the matrix.
#' @param n the columns of the matrix.
#' @param r the rank of the matrix.
#' @param snr the signal-to-noise ratio in generating the matrix. Default \code{snr = 3}.
#' @param prop the proportion of missing observations. Default \code{prop = 0.5}.
#' @param seed the random seed. Default \code{seed = 1}.
#'
#' @details We generate the matrix by \eqn{UV + \epsilon}, where \eqn{U}, \eqn{V} are \eqn{m} by \eqn{r}, \eqn{r} by \eqn{n} matrix satisfy standard normal
#' distribution. \eqn{\epsilon} has a normal distribution with mean 0 and variance \eqn{\frac{r}{snr}}.
#'
#' @return A matrix with missing values.
#'
#' @export
#' @importFrom stats rnorm
#'
#' @examples
#' m <- 100
#' n <- 100
#' r <- 10
#' x_na <- incomplete.generator(m, n, r)
#' head(x_na[, 1:6])
incomplete.generator <- function(m, n, r, snr = 3, prop = 0.5, seed = 1) {
  prop <- 1 - prop
  set.seed(seed)
  U <- matrix(rnorm(m * r), nrow = m)
  V <- matrix(rnorm(n * r), ncol = n)
  epsilon <- matrix(rnorm(m * n, sd = sqrt(r) / snr), nrow = m)
  X_true <- U %*% V
  X <- X_true + epsilon
  s <- sample(0 : (m * n - 1), round(m * n * prop)) # index of observed value
  s_c <- (0 : (m * n - 1))[-(s + 1)] # index of filling in the missing value
  Xna <- X
  Xna[s_c + 1] <- NA
  Xna
}
