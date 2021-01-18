# Part of the landsepi R package.
# Copyright (C) 2017 Loup Rimbaud <loup.rimbaud@inrae.fr>
#                    Julien Papaix <julien.papaix@inrae.fr>
#                    Jean-François Rey <jean-francois.rey@inrae.fr>
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software Foundation, Inc.,i
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
#


#' @title setSeedValue
#' @name setSeedValue
#' @description Set RNG seed to seed value if not NULL, otherwise set it to timestamps value
#' @details Sets seed for "Mersenne-Twister" algorithm using Inversion generation
#' @param seed an interger as seed value or NULL
#' @return the new seed value for RNG
#' @examples
#' setSeedValue(seed = 10)
#' @export
setSeedValue <- function(seed = NULL) {
  if (is.null(seed)) {
    newseed <- as.numeric(Sys.time())
  }
  else {
    newseed <- seed
  }
  
  ## insure newseed is an integer
  newseed <- trunc(newseed)
  
  # set seed for RNG
  # if( getRversion >= "3.6.0" ) RNGkind(sample.kind = "Rounding")
  set.seed(newseed, kind = "Mersenne-Twister", normal.kind = "Inversion")
  message("Seed for RGN set to : ", newseed)

  return(newseed)
}


#' @title is.wholenumber
#' @name is.wholenumber
#' @description Tests if a number or vector is a whole number
#' @param x a number or vector or matrix
#' @param tol double tolerance
#' @return a logical of the same format as x
#' @examples
#' is.wholenumber(-5)
#' is.wholenumber(10)
#' is.wholenumber(2.5)
#' is.wholenumber(matrix(1:9, nrow=3))
#' @export
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) {
    return( abs(x - round(x)) < tol )
}


#' @title is.positive
#' @name is.positive
#' @description Tests if a number or vector is positive (including 0)
#' @param x a number or vector or matrix
#' @return a logical of the same size as x
#' @examples
#' is.positive(-5)
#' is.positive(10)
#' is.positive(2.5)
#' is.positive(matrix(1:9, nrow=3))
#' @export
is.positive <- function(x) {
  start <- 1
  end <- length(x) + 1
  result <- c()
  while (start < end) {
    y <- x[start]
    if (y >= 0) {
      result <- c(result, TRUE)
    } else {
      result <- c(result, FALSE)
    }
    start <- start + 1
  }
  return(result)
}

#' @title is.strict.positive
#' @name is.strict.positive
#' @description Tests if a number or vector is strictly positive (i.e. excluding 0)
#' @param x a number or vector or matrix
#' @return a logical of the same size as x
#' @examples
#' is.strict.positive(-5)
#' is.strict.positive(10)
#' is.strict.positive(2.5)
#' is.strict.positive(matrix(1:9, nrow=3))
#' @export
is.strict.positive <- function(x) {
  start <- 1
  end <- length(x) + 1
  result <- c()
  while (start < end) {
    y <- x[start]
    if (y > 0) {
      result <- c(result, TRUE)
    } else {
      result <- c(result, FALSE)
    }
    start <- start + 1
  }
  return(result)
}


#' @title is.in.01
#' @name is.in.01
#' @description Tests if a number or vector is in the interval \[0,1\]
#' @param x a number or vector or matrix
#' @param exclude0 TRUE is 0 is excluded, FALSE otherwise (default)
#' @return a logical of the same size as x
#' @examples
#' is.in.01(-5)
#' is.in.01(0)
#' is.in.01(1)
#' is.in.01(0, exclude0 = TRUE)
#' is.in.01(2.5)
#' is.in.01(matrix(5:13/10, nrow=3))
#' @export
is.in.01 <- function(x, exclude0 = FALSE){
  if (exclude0){
    return (x > 0 & x <= 1)
  } else {
    return (x >= 0 & x <= 1)
  }
}



