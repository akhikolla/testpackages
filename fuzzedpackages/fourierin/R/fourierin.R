## -------------------------------------------------------------------
#' Univariate Fourier integrals
#'
#' It computes Fourier integrals of functions of one and two
#' variables on a regular grid.
#'
#' See vignette for more detailed examples.
#'
#' @param f function or a vector of size m. If a function is
#'     provided, it must be able to be evaluated at vectors. If a
#'     vector of values is provided, such evaluations must have been
#'     obtained on a regular grid and the Fourier integral is faster
#'     is m is a power of 2.
#' @param lower_int Lower integration limit(s).
#' @param upper_int Upper integration limit(s).
#' @param lower_eval Lower evaluation limit(s). It can be NULL if an
#'     evaluation grid is provided.
#' @param upper_eval Upper evaluation limit(s). It can be NULL if an
#'     evaluation grid is provided.
#' @param const_adj Factor related to adjust definition of Fourier
#'     transform. It is usually equal to 0, -1 or 1.
#' @param freq_adj Constant to adjust the exponent on the definition
#'     of the Fourier transform. It is usually equal to 1, -1, 2pi or
#'     -2pi.
#' @param resolution A vector of integers (faster if powers of two)
#'     determining the resolution of the evaluation grid. Not
#'     required if f is a vector.
#' @param eval_grid Optional matrix with d columns with the points
#'     where the Fourier integral will be evaluated. If it is
#'     provided, the FFT will not be used.
#' @param use_fft Logical value specifying whether the FFT will be
#'     used.
#' @return If w is given, only the values of the Fourier integral are
#'     returned, otherwise, a list with the elements \item{w}{A vector
#'     of size m where the integral was computed.}  \item{values}{A
#'     complex vector of size m with the values of the integral}
#'
#' @examples
#' ##--- Example 1 ---------------------------------------------------
#' ##--- Recovering std. normal from its characteristic function -----
#' library(fourierin)
#'
#' #' Function to to be used in integrand
#' myfun <- function(t) exp(-t^2/2)
#'
#'                                         # Compute Foueien integral
#' out <- fourierin_1d(f = myfun,
#'                     lower_int = -5, upper_int = 5,
#'                     lower_eval = -3, upper_eval = 3,
#'                     const_adj = -1, freq_adj = -1,
#'                     resolution = 64)
#'
#' ## Extract grid and values
#' grid <- out$w
#' values <- Re(out$values)
#'
#' plot(grid, values, type = "l", col = 3)
#' lines(grid, dnorm(grid), col = 4)
#'
#' ##--- Example 2 -----------------------------------------------
#' ##--- Computing characteristic function of a gamma r. v. ------
#'
#' library(fourierin)
#'
#' ## Function to to be used in integrand
#' myfun <- function(t) dgamma(t, shape, rate)
#'
#' ## Compute integral
#' shape <- 5
#' rate <- 3
#' out <- fourierin_1d(f = myfun, lower_int = 0, upper_int = 6,
#'                     lower_eval = -4, upper_eval = 4,
#'                     const_adj = 1, freq_adj = 1, resolution = 64)
#'
#' grid <- out$w                           # Extract grid
#' re_values <- Re(out$values)             # Real values
#' im_values <- Im(out$values)             # Imag values
#'
#'                                      # Now compute the real and
#'                                      # imaginary true values of the
#'                                         # characteric function.
#' true_cf <- function(t, shape, rate) (1 - 1i*t/rate)^-shape
#' true_re <- Re(true_cf(grid, shape, rate))
#' true_im <- Im(true_cf(grid, shape, rate))
#'
#'                                      # Compare them. We can see a
#'                                      # slight discrepancy on the
#'                                      # tails, but that is fixed
#'                                      # when resulution is
#'                                      # increased.
#' plot(grid, re_values, type = "l", col = 3)
#' lines(grid, true_re, col = 4)
#'
#'                                         # Same here
#' plot(grid, im_values, type = "l", col = 3)
#' lines(grid, true_im, col = 4)
#' @export
fourierin_1d <- function(f, lower_int, upper_int,
                         lower_eval = NULL, upper_eval = NULL,
                         const_adj, freq_adj, resolution = NULL,
                         eval_grid = NULL, use_fft = TRUE) {
    ## Condensed notation
    a <- lower_int
    b <- upper_int
    c <- lower_eval
    d <- upper_eval
    r <- const_adj
    s <- freq_adj
    resol <- resolution
    w <- eval_grid

## fourierin_1d <- function(f, a, b, c = NULL, d = NULL,
##                          r, s, resol = NULL, w = NULL,
##                          use_fft = TRUE) {
    ## Flag to determine wheter a list or a vector will be returned
    w_given <- !is.null(w)

    ## If function values are provided, then the resolution
    ## is the length of the vector of values.
    if (!is.function(f)) resol <- length(f)

    ## Increment in the frequency domain.
    gam <- (d - c)/resol

    ## Freq. dom. vector. If w is provided, FFT will NOT be used.
    if (is.null(w)) {
        if (is.null(c) | is.null(d)) {
            stop("c and d must be provided.")
        }
        w <- seq(c, d - gam, length.out = resol)
    } else {
        use_fft <- FALSE
    }

    ## If f is the function, it needs to be evaluated in
    ## the time domain values.
    if (is.function(f)) {

        del <- (b - a)/resol # Increment in the time domain.
        t <- seq(a + del/2, b - del/2,
                 length.out = resol)    # Freq. dom. vector.
        f_t <- f(t)                     # Function values
        ## Rutinary check
        if(is.null(f_t)) stop("Function f is null.")
    } else {
        f_t <- f
    }

    if (!use_fft) {
        out <- switch(is.complex(f_t) + 1,
                      fourierin_1d_nonregular_cpp(f_t, a, b, w,
                                                  resol, r, s),
                      fourierin_cx_1d_nonregular_cpp(f_t,
					 a, b, w, resol, r, s))
    } else {
        out <- switch(is.complex(f_t) + 1,
                      fourierin_1d_cpp(f_t, a, b, c, d, r, s),
                      fourierin_cx_1d_cpp(f_t, a, b, c, d, r, s))
    }

    ## If w is given, return only the values of the integral,
    ## otherwise alse return w.
    if(w_given) return(out)

    return(list(w = w,                  # Return list.
                values = drop(out)))
}

#' Bivariate Fourier integrals
#'
#' It computes Fourier integrals for functions of one and two
#' variables.
#'
#' @inheritParams fourierin_1d
#' @return If w is given, only the values of the Fourier integral are
#'     returned, otherwise, a list with three elements
#'     \item{w1}{Evaluation grid for first entry} \item{w2}{Evaluation
#'     grid for second entry} \item{values}{m1 x m2 matrix of complex
#'     numbers, corresponding to the evaluations of the integral}
#'
#' @examples
#' ##--- Recovering std. normal from its characteristic function -----
#' library(fourierin)
#'
#' ##-Parameters of bivariate normal distribution
#' mu <- c(-1, 1)
#' sig <- matrix(c(3, -1, -1, 2), 2, 2)
#'
#' ##-Multivariate normal density
#' ##-x is n x d
#' f <- function(x) {
#'     ##-Auxiliar values
#'     d <- ncol(x)
#'     z <- sweep(x, 2, mu, "-")
#'
#'     ##-Get numerator and denominator of normal density
#'     num <- exp(-0.5*rowSums(z * (z %*% solve(sig))))
#'     denom <- sqrt((2*pi)^d*det(sig))
#'
#'     return(num/denom)
#' }
#'
#' ##-Characteristic function
#' ##-s is n x d
#' phi <- function(s) {
#'     complex(modulus = exp(- 0.5*rowSums(s*(s %*% sig))),
#'             argument = s %*% mu)
#' }
#'
#' ##-Approximate cf using Fourier integrals
#' eval <- fourierin_2d(f, lower_int = c(-8, -6), upper_int = c(6, 8),
#'                      lower_eval = c(-4, -4), upper_eval = c(4, 4),
#'                      const_adj = 1, freq_adj =  1,
#'                      resolution = c(128, 128))
#'
#' ## Extract values
#' t1 <- eval$w1
#' t2 <- eval$w2
#' t <- as.matrix(expand.grid(t1 = t1, t2 = t2))
#' approx <- eval$values
#' true <- matrix(phi(t), 128, 128)        # Compute true values
#'
#' ##-This is a section of the characteristic functions
#' i <- 65
#' plot(t2, Re(approx[i, ]), type = "l", col = 2,
#'      ylab = "",
#'      xlab = expression(t[2]),
#'      main = expression(paste("Real part section at ",
#'                              t[1], "= 0")))
#' lines(t2, Re(true[i, ]), col = 3)
#' legend("topleft", legend = c("true", "approximation"),
#'        col = 3:2, lwd = 1)
#'
#' ##-Another section, now of the imaginary part
#' plot(t1, Im(approx[, i]), type = "l", col = 2,
#'      ylab = "",
#'      xlab = expression(t[1]),
#'      main = expression(paste("Imaginary part section at ",
#'                              t[2], "= 0")))
#' lines(t1, Im(true[, i]), col = 3)
#' legend("topleft", legend = c("true", "approximation"),
#'        col = 3:2, lwd = 1)
#' @export
fourierin_2d <- function(f, lower_int, upper_int,
                         lower_eval = NULL, upper_eval = NULL,
                         const_adj, freq_adj,
                         resolution = NULL, eval_grid = NULL,
                         use_fft = TRUE){
    ## Condensed notation
    a <- lower_int
    b <- upper_int
    c <- lower_eval
    d <- upper_eval
    r <- const_adj
    s <- freq_adj
    resol <- resolution
    w <- eval_grid

    ## If function values are provided, then the resolution is the
    ## length of the vector of values.
    if(!is.function(f)) resol <- dim(f)

    ## Increment in the frequency domain.
    gam <- (d - c)/resol

    ## Freq. dom. vectors.
    if (is.null(w)){
        if (is.null(c) | is.null(d)) {
            stop("c and d must be provided.")
            }
        w1 <- seq(c[1], d[1] - gam[1], length.out = resol[1])
        w2 <- seq(c[2], d[2] - gam[2], length.out = resol[2])
    }

    ## If f is the function, it needs to be evaluated in the time
    ## domain values.
    if(is.function(f)){
        del <- (b - a)/resol # Increment in the time domain.
        t1 <- seq(a[1] + del[1]/2, b[1] - del[1]/2,
                  length.out = resol[1]) # Freq. dom. vector.
        t2 <- seq(a[2] + del[2]/2, b[2] - del[2]/2,
                  length.out = resol[1]) # Freq. dom. vector.
        t <- as.matrix(expand.grid(t1, t2))
        f_vals <- matrix(f(t), resol[1], resol[2])

                                        # Rutinary check
        if(is.null(f_vals)) stop("Function f is null.")

    } else{
        f_vals <- f
    }

    if (!use_fft) {
        w_temp <- switch(is.null(w) + 1,
                         w,
                         as.matrix(expand.grid(w1, w2))
                         )

        out <- switch(is.complex(f_vals) + 1,
                      fourierin_2d_nonregular_cpp(f_vals, a, b,
                                                  w_temp,
                                                  resol, r, s),
                      fourierin_cx_2d_nonregular_cpp(f_vals, a, b,
                                                     w_temp, resol,
                                                     r, s))
        ## If no grid was provided, the values of the Fourier integral
        ## are put into a matrix.
        if(is.null(w)) out <- matrix(out, resol[1], resol[2])
    } else {
        out <- switch(is.complex(f_vals) + 1,
                      fourierin_2d_cpp(f_vals, a, b, c, d, r, s),
                      fourierin_cx_2d_cpp(f_vals, a, b, c, d, r, s))
    }


    ## Return grid if it was not provided.
    if (is.null(w)) {
      return(list(w1 = w1,
                  w2 = w2,
                  values = out))

    } else {
      return(out)
    }
}



#' Compute Fourier integrals
#'
#' It computes Fourier integrals for functions of one and two
#' variables.
#'
#' See plenty of detailed examples in the vignette.
#'
#' @inheritParams fourierin_1d
#' @return A list with the elements n-dimensional array and n vectors
#'     with their corresponding resolution. Specifically,
#'     \item{values}{A n-dimensional (resol_1 x resol_2 x ... x
#'     resol_n) complex array with the values.}  \item{w1}{A vector of
#'     size resol_1} \item{...}{ } \item{wn}{A vector of size resol_n}
#'
#' @examples
#' ##--- Example 1 ---------------------------------------------------
#' ##--- Recovering std. normal from its characteristic function -----
#' library(fourierin)
#'
#' ## Function to be used in the integrand
#' myfnc <- function(t) exp(-t^2/2)
#'
#' ## Compute integral
#' out <- fourierin(f = myfnc, lower_int = -5, upper_int = 5,
#'                  lower_eval= -3, upper_eval = 3, const_adj = -1,
#'                  freq_adj = -1, resolution = 64)
#'
#' ## Extract grid and values
#' grid <- out$w
#' values <- Re(out$values)
#'
#' ## Compare with true values of Fourier transform
#' plot(grid, values, type = "l", col = 3)
#' lines(grid, dnorm(grid), col = 4)
#'
#'
#' ##--- Example 2 ---------------------------------------------------
#' ##--- Computing characteristic function of a gamma r. v. ----------
#'
#' library(fourierin)
#'
#' ## Function to be used in integrand
#' myfnc <- function(t) dgamma(t, shape, rate)
#'
#' ## Compute integral
#' shape <- 5
#' rate <- 3
#' out <- fourierin(f = myfnc, lower_int = 0, upper_int = 6,
#'                  lower_eval = -4, upper_eval = 4,
#'                  const_adj = 1, freq_adj = 1, resolution = 64)
#'
#' ## Extract values
#' grid <- out$w                           # Extract grid
#' re_values <- Re(out$values)             # Real values
#' im_values <- Im(out$values)             # Imag values
#'
#' ## Now compute the real and imaginary true values of the
#' ## characteric function.
#' true_cf <- function(t, shape, rate) (1 - 1i*t/rate)^-shape
#' true_re <- Re(true_cf(grid, shape, rate))
#' true_im <- Im(true_cf(grid, shape, rate))
#'
#' ## Compare them. We can see a slight discrepancy on the tails,
#' ## but that is fixed when resulution is increased.
#' plot(grid, re_values, type = "l", col = 3)
#' lines(grid, true_re, col = 4)
#'
#'                                         # Same here
#' plot(grid, im_values, type = "l", col = 3)
#' lines(grid, true_im, col = 4)
#'
#' ##--- Example 3 -------------------------------------------------
#' ##--- Recovering std. normal from its characteristic function ---
#' library(fourierin)
#'
#' ##-Parameters of bivariate normal distribution
#' mu <- c(-1, 1)
#' sig <- matrix(c(3, -1, -1, 2), 2, 2)
#'
#' ##-Multivariate normal density
#' ##-x is n x d
#' f <- function(x) {
#'     ##-Auxiliar values
#'     d <- ncol(x)
#'     z <- sweep(x, 2, mu, "-")
#'     ##-Get numerator and denominator of normal density
#'     num <- exp(-0.5*rowSums(z * (z %*% solve(sig))))
#'     denom <- sqrt((2*pi)^d*det(sig))
#'     return(num/denom)
#' }
#'
#' ## Characteristic function
#' ## s is n x d
#' phi <- function(s) {
#'     complex(modulus = exp(- 0.5*rowSums(s*(s %*% sig))),
#'             argument = s %*% mu)
#' }
#'
#' ##-Approximate cf using Fourier integrals
#' eval <- fourierin(f, lower_int = c(-8, -6), upper_int = c(6, 8),
#'                   lower_eval = c(-4, -4), upper_eval = c(4, 4),
#'                   const_adj = 1, freq_adj =  1,
#'                   resolution = c(128, 128))
#'
#' ## Extract values
#' t1 <- eval$w1
#' t2 <- eval$w2
#' t <- as.matrix(expand.grid(t1 = t1, t2 = t2))
#' approx <- eval$values
#' true <- matrix(phi(t), 128, 128)        # Compute true values
#'
#'
#' ## This is a section of the characteristic function
#' i <- 65
#' plot(t2, Re(approx[i, ]), type = "l", col = 2,
#'      ylab = "",
#'      xlab = expression(t[2]),
#'      main = expression(paste("Real part section at ",
#'                              t[1], "= 0")))
#' lines(t2, Re(true[i, ]), col = 3)
#' legend("topleft", legend = c("true", "approximation"),
#'        col = 3:2, lwd = 1)
#'
#' ##-Another section, now of the imaginary part
#' plot(t1, Im(approx[, i]), type = "l", col = 2,
#'      ylab = "",
#'      xlab = expression(t[1]),
#'      main = expression(paste("Imaginary part section at ",
#'                              t[2], "= 0")))
#' lines(t1, Im(true[, i]), col = 3)
#' legend("topleft", legend = c("true", "approximation"),
#'        col = 3:2, lwd = 1)
#' @export
fourierin <- function(f, lower_int, upper_int,
                      lower_eval = NULL, upper_eval = NULL,
                      const_adj, freq_adj,
                      resolution = NULL, eval_grid = NULL,
                      use_fft = TRUE){

    ## Get dimension of function from lower integration limit.
    n <- length(lower_int)

    switch(n,
           return(fourierin_1d(f, lower_int, upper_int,
                               lower_eval, upper_eval,
                               const_adj, freq_adj, resolution,
                               eval_grid, use_fft)),
           return(fourierin_2d(f, lower_int, upper_int,
                               lower_eval, upper_eval,
                               const_adj, freq_adj,
                               resolution, eval_grid, use_fft))
           ) # End switch
}

#' @useDynLib fourierin, .registration=TRUE
#' @importFrom Rcpp sourceCpp
NULL


.onUnload <- function (libpath) {
  library.dynam.unload("fourierin", libpath)
}
