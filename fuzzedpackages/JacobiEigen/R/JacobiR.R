##' The Jacobi Algorithm
##'
##' Eigenvalues and optionally, eigenvectore of a real symmetric matrix using the
##' classical Jacobi algorithm, (Jacobi, 1854)
##' @title The Jacobi Algorithm in Pure R
##' @param x a real symmetric matrix
##' @param symmetric a logical value.  Is the matrix symmetric?  (Only symmetric matrices are allowed.)
##' @param only.values A logical value: Do you want only the eigenvalues?
##' @param eps a small positive error tolerance
##' @export JacobiR
##' @examples
##' V <- crossprod(matrix(rnorm(25), 5))
##' JacobiR(V)
##' identical(Jacobi(V), JacobiR(V))
##' all.equal(Jacobi(V)$values, base::eigen(V)$values)
##' @return a list of two components as for \code{base::eigen}
JacobiR <- function(x, symmetric = TRUE, only.values = FALSE,
                    eps = if(!only.values) .Machine$double.eps else
                      sqrt(.Machine$double.eps)) {
  if(!symmetric) 
    stop("only real symmetric matrices are allowed")
  n <- nrow(x)
  H <- if(only.values) NULL else diag(n)
  eps <- max(eps, .Machine$double.eps)

  if(n > 1) {
    lt <- which(lower.tri(x))

    repeat {
      k <- lt[which.max(abs(x[lt]))]  ## the matrix element
      j <- floor(1 + (k - 2)/(n + 1)) ## the column
      i <- k - n * (j - 1)            ## the row

      if(abs(x[i, j]) < eps) break

      Si <- x[, i]
      Sj <- x[, j]

      theta <- 0.5*atan2(2*Si[j], Sj[j] - Si[i])
      c <- cos(theta)
      s <- sin(theta)

      x[i, ] <- x[, i] <- c*Si - s*Sj
      x[j, ] <- x[, j] <- s*Si + c*Sj
      x[i,j] <- x[j,i] <- 0
      x[i,i] <- c^2*Si[i] - 2*s*c*Si[j] + s^2*Sj[j]
      x[j,j] <- s^2*Si[i] + 2*s*c*Si[j] + c^2*Sj[j]
      if(!only.values) {
        Hi <- H[, i]
        H[, i] <- c*Hi - s*H[, j]
        H[, j] <- s*Hi + c*H[, j]
      }
    }
  }
  out <- structure(list(values = as.vector(diag(x)), vectors = H), 
                   class = c("Jacobi", "eigen"))
  o <- order(out$values, decreasing = TRUE)
  out$values <- out$values[o]
  if(!only.values) out$vectors <- out$vectors[, o]
  out
}

##' The Classical Jacobi Algorithm
##'
##' Eigenvalues and optionally, eigenvectore, of a real symmetric matrix using the
##' classical Jacobi algorithm, (Jacobi, 1854)
##' @import Rcpp
##' @useDynLib JacobiEigen
##' @title The Jacobi Algorithm using Rcpp
##' @param x A real symmetric matrix
##' @param symmetric a logical value.  Is the matrix symmetric?  (Only symmetric matrices are allowed.)
##' @param only.values A logical value: do you want only the eigenvalues?
##' @param eps an error tolerance. 0.0 implies \code{.Machine$double.eps} and
##'   \code{sqrt(.Machine$double.eps)} if \code{only.values = TRUE}
##' @export Jacobi
##' @examples
##' V <- crossprod(matrix(runif(40, -1, 1), 8))
##' Jacobi(V)
##' identical(Jacobi(V), JacobiR(V))
##' all.equal(Jacobi(V)$values, base::eigen(V)$values)
##' @return a list of two components as for \code{base::eigen}
Jacobi <- function(x, symmetric = TRUE, only.values = FALSE, eps = 0.0) {
  if(!symmetric) 
    stop("only real symmetric matrices are allowed")
  out <- structure(JacobiCpp(x, only.values, eps), 
                   class = c("Jacobi", "eigen"))
  o <- order(out$values, decreasing = TRUE)
  out$values <- out$values[o]
  if(!only.values) out$vectors <- out$vectors[, o]
  out
}

##' The Classical Jacobi Algorithm with a stagewise protocol
##'
##' Eigenvalues and optionally, eigenvectore, of a real symmetric matrix using the
##' classical Jacobi algorithm, (Jacobi, 1846) using a stagewise rotation protocol
##' @title The Jacobi Algorithm using Rcpp with a stagewise rotation protocol 
##' @param x A real symmetric matrix
##' @param symmetric a logical value.  Is the matrix symmetric?  (Only symmetric matrices are allowed.)
##' @param only.values A logical value: do you want only the eigenvalues?
##' @param eps an error tolerance. 0.0 implies \code{.Machine$double.eps} and
##'   \code{sqrt(.Machine$double.eps)} if \code{only.values = TRUE}
##' @export JacobiS
##' @examples
##' V <- crossprod(matrix(runif(40, -1, 1), 8))
##' JacobiS(V)
##' all.equal(JacobiS(V)$values, Jacobi(V)$values)
##' zapsmall(crossprod(JacobiS(V)$vectors, Jacobi(V)$vectors))
##' @return a list of two components as for \code{base::eigen}
JacobiS <- function(x, symmetric = TRUE, only.values = FALSE, eps = 0.0) {
  if(!symmetric) 
    stop("only real symmetric matrices are allowed")
  out <- structure(JacobiSCpp(x, only.values, eps), 
                   class = c("Jacobi", "eigen"))
  o <- order(out$values, decreasing = TRUE)
  out$values <- out$values[o]
  if(!only.values) out$vectors <- out$vectors[, o]
  out
}

##' @export
print.Jacobi <- function (x, ...) {
  cat("eigen decomposition (Jacobi algorithm)\n")
  print(unclass(x), ...)
  invisible(x)
}

.onUnload <- function(libpath) {
  library.dynam.unload("JacobiEigen", libpath)
}
