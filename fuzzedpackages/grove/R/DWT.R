#' Discrete wavelet transform  
#'
#' This function performs the discrete wavelet transform (DWT) according to
#' Mallat's pyramidal algorithm (Mallat, 1989).
#'
#' @param data A matrix of data, where each row is an observation. 
#' The number of columns must be a power of two.
#' @param filter.number The smoothness of the wavelet to use in the 
#' decomposition. 
#' @param family The family of wavelets. The two most common options 
#' are \code{DaubExPhase} and \code{DaubLeAsymm}.
#' @details
#' See function \code{wd} from package \code{wavethresh} for more details.
#' @return A \code{DWT} object. This object is a list with the following 
#' components: 
#' @export
#' @examples
#' data <- GenerateSyntheticAnova(st.dev = 5, n.replicates = 10)
#' W <- DWT(data$noisy.Y)

DWT <- function(data, 
                filter.number = 10, 
                family = "DaubLeAsymm") {
  
  if (sum(methods::is(data) == c("numeric", "vector")) == 2) {
    data <- matrix(data, nrow = 1)
  }
  
  J <- ncol(data)
  n <- nrow(data)
  D <- matrix(NA, nrow = n, ncol = J - 1 )
  C <- rep(NA, n)
  for (i in 1:n) {
    temp <- wd(data[i, ], filter.number = filter.number, family = family)
    D[i, ] <- rev(temp$D)
    C[i] <- accessC(temp, level = 0)
  }
  output <- list(C = C, 
                 D = D, 
                 J = log2(J), 
                 filter.number = filter.number, 
                 family = family)
  class(output) <- "DWT"
  return(output)
}