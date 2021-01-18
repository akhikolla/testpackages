
###################################################################################
##' Generator
##'
##' @description
##' This function generates sequence from a MHMM model
##' 
##' @param n numeric, number of subjects.
##' @param nbT numeric, length of the sequence.
##' @param A list, matrices of the transition probabilities per class
##' @param delta numeric, proportions of the classes.
##' @param a numeric, shapes of the gamma distributions.
##' @param b numeric, rates of the gamma distributions.
##' @param eps numeric, proportions of the zero-inflated in the emission laws.
##' @return Returns a list of sequences.
##'
##'
##' @examples
##' ech <- rdata.mhmm(25, 10)
##' res <- mhmm(ech$y, K = 2, M = 4, nbcores = 1, nbinit = 5, iterSmall = 2)
##' @export rdata.mhmm
rdata.mhmm <- function(n = 50,
                       nbT = 40,
                       A = list(matrix(c(.7, .3, .1, .2, .6, .1, .1, .1, .8), 3, 3), matrix(1/3,3,3)),
                       delta = rep(1/length(A), length(A)),
                       a = c(1,2,3)*10,
                       b = c(1,1,1),
                       eps = c(.5,.1,.2)*0){
  M <- ncol(A[[1]])
  z <- sample(1:length(delta), n, replace = TRUE, prob = delta)
  x <- y <- list()
  for (i in 1:n){
    x[[i]] <- y[[i]] <- rep(0, nbT)
    pistart <- abs(eigen(A[[2]])$vectors[,1])
    pistart <- pistart / sum(pistart)
    x[[i]][1] <- sample(1:M, 1, prob = pistart)
    y[[i]][1] <- rnorm(1, a[x[[i]][1]], b[x[[i]][1]])
    for (t in 2:nbT){
      x[[i]][t] <- sample(1:M, 1, prob = A[[z[i]]][x[[i]][t-1],])
      y[[i]][t] <- rgamma(1, a[x[[i]][t]], b[x[[i]][t]]) * (runif(1)>eps[x[[i]][t]])
    }
  }
  param <- new("mhmmparam", K = length(A), M = M, A = A,
               delta = delta, pi = matrix(1/2, length(A), length(a)),
               lambda=list(eps= as.numeric(eps), a = as.numeric(a), b = as.numeric(b)))
  list(x=x, y=y, z=z, param=param)
}
