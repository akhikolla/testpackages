#' Concatenate strings
#' 
#' Concatenate Strings
#' 
#' The same result may be achieved with \code{paste}, but in some circumstance this function is more user friendly.
#'
#' @rdname ConcatenateStrings
#' @return A character string of the concatenated values.
#' @seealso \code{\link{paste}}
#' @param x asd
#' @param y asd 
#' @export
#' @examples
#' paste('I ', 'love ', 'R.', sep='')
#' 'I ' %.% 'love ' %.% 'R.'
#'
#' x = c(2,1,6,7,9)
#' paste('The length of vector (', paste(x , sep='', collapse =','), ') is ', length(x) , sep='')
#' 'The length of vector (' %.% paste(x , sep='', collapse =',') %.% ') is ' %.% length(x)
`%.%` <- function(x, y) paste(x, y, sep = "")



######### Generating functions #########

#' Random process generators
#' 
#' Generate a trajectory of random processes.  
#' 
#' \code{rwiener} generate Wiener process via partial sums process and 
#' \code{rbridge} generate Brownian bridge via \code{rwiener}.
#' The original code of \code{rwiener} and \code{rbridge} was written in the package \code{e1071}. 
#' In this package these functions was modified to
#' include leading zero in the beginning of the sample. 
#' 
#' \code{rcumbin} generate partial sums process from random variables with values \code{-1,  0, 1}.
#' 
#' @rdname ProcessGenerators
#' @return A time series containing a simulated realization of random processes. 
#' The length of time series is \code{frequency+1}, since zero is always included in the beginning of the sample.
#' @param frequency a number specifying the size of trajectory vector. The trajectory will start at point 0 
#' and will have \code{frequency} more observations. The length of the results will be \code{frequency+1} .  
#' @param end a number. The end point of the process in the 'time' scale.
#' @export
rwiener <- function(frequency = 1000, end = 1) {
  z <- c(0, cumsum(rnorm(end * frequency)/sqrt(frequency)))
  ts(z, start = 0, end = 1, frequency = frequency)
}


#' @rdname ProcessGenerators
#' @export
rbridge <- function(frequency = 1000, end = 1) {
  z <- rwiener(frequency = frequency, end = end)
  ts(z - time(z) * as.vector(z)[frequency], start = 0, frequency = frequency)
}

#' @rdname ProcessGenerators
#' @export
rcumbin <- function(frequency = 1000, end = 1) {
  z <- c(0, cumsum(sample(c(-1, 0, 1), frequency, replace = TRUE)))
  ts(z, start = 0, end = 1, frequency = frequency)
} 
