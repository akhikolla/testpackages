lowpassFilter <- function(type = c("bessel"), param, sr = 1, len = NULL, shift = 0.5) {
  type <- match.arg(type)
  
  if (!is.list(param)) {
    stop("param must be a list")
  }
  
  if (!is.numeric(sr) || length(sr) != 1 || !is.finite(sr) || sr <= 0) {
    stop("sr must be a single positive finite numeric")
  }
  
  if (!is.null(len)) {
    if (!is.numeric(len) || length(len) != 1 || !is.finite(len)) {
      stop("len must be a single positive integer")
    }
    
    if (!is.integer(len)) {
      len <- as.integer(len + 1e-6)
    }
    
    if (len <= 0L) {
      stop("len must be a single positive integer")
    }
  }
  
  if (!is.numeric(shift) || length(shift) != 1 || !is.finite(shift) || shift < 0 || shift > 1) {
    stop("shift must be a single finite numeric between 0 and 1")
  }
  
  ret <- list(type = type, param = param, sr = sr, len = len)
  
  switch(type,
         bessel = {
           if (!all(c("pole", "cutoff") %in% names(param))) {
             stop("for type 'bessel' param must be a list with entries 'pole' and 'cutoff'")
           }
           
           if (!all(names(param) %in% c("pole", "cutoff"))) {
             stop(paste("param contains unused entries, ",
                        "for type 'bessel' only the entries 'pole' and 'cutoff' should be given",
                        sep = ""))
           }
           
           if (!is.numeric(param$pole) || length(param$pole) != 1 || !is.finite(param$pole)) {
             stop("param$pole must be a single positive integer")
           }
           
           if (!is.integer(param$pole)) {
             param$pole <- as.integer(param$pole + 1e-6)
           }
           
           if (param$pole <= 0L) {
             stop("param$pole must be a single positive integer")
           }
           
           if (!is.numeric(param$cutoff) || length(param$cutoff) != 1 || !is.finite(param$cutoff) ||
               param$cutoff <= 0 || param$cutoff > 1) {
             stop("param$cutoff must be a single positive finite numeric smaller than or equal to 1")
           }
           
           # ensure a fixed order
           ret$param <- list(pole = param$pole, cutoff = param$cutoff)

           # filter coefficients
           a <- .BesselPolynomial(param$pole, reverse = TRUE)
           # zeros
           r <- polyroot(a)
           # coefficients of Heaviside functions
           p <- sapply(seq(along = r), function(i) 1 / prod(r[i] - r[-i]))
           # power coefficients
           A2 <- a * 1i^(seq(along = a) - 1)
           A <- sapply(1:(2 * length(A2) - 1), function(i) {
             j <- max(1, i - length(A2) + 1):min(i, length(A2))
             sum(A2[j] * Conj(A2[i + 1 - j]))
           })
           # compute cut-off frequency of "default" filter, i.e. where power is halved
           omega0 <- polyroot(A / a[1]^2 - c(2, rep(0, 2 * length(A2) - 2)))
           omega0 <- Re(omega0[which.min(abs(Arg(omega0)))])
           tau <- param$cutoff / omega0 * 2 * pi
           # kernel function
           ret$kernfun <- function(t)
             a[1] * tau * sr * Re(sapply(t * tau * sr, function(s) {
               if (s <= 0) {return(0)}
               sum(p * exp(r * s))
             }
             ))
           # roots and coefficients for step response
           rs <- c(0, r)
           ps <- sapply(seq(along = rs), function(i) 1 / prod(rs[i] - rs[-i]))
           # step response
           ret$stepfun <- function(t) a[1] *
             Re(sapply(t * tau * sr, function(s) {
               if (s <= 0) {return(0)}
               sum(ps * exp(rs * s))
             }))
           # auto-correlation, note that the auto-correlation at lag s > 0 of Heavy(t) p exp(r t)
           # are given by the sum over all i,j of
           # integral Heavy(t) pi exp(ri t) Heavy(t+s) pj exp(rj (t+s)) dt = 
           # integral Heavy(t) pi pj exp(rj s) exp( (ri + rj) t) = pi pj exp(rj s) / (ri + rj)
           # and for s < 0 similarly pi pj exp(ri s) / (ri + rj)
           K2 <- -Re(sum(outer(seq(along = p), seq(along = p), function(i, j) p[i] * p[j] / (r[i] + r[j]))))
           ret$acfun <- function(t) {
             acf <- sapply(t * tau * sr,
                           function(s) -Re(sum(outer(seq(along = p), seq(along = p), function(i, j) 
                             p[i] * p[j] * exp(r[j] * abs(s)) / (r[i] + r[j])))))
             acf / K2
           }
           ret$acAntiderivative <- function(t, lag) {
             if (length(t) == 1L) {
               t <- rep(t, length(lag))
             }
             if (length(lag) == 1L) {
               lag <- rep(lag, length(t))
             }
             
             acf <- numeric(length(t))
             compute <- t > 0
             
             for (k in seq(along = lag)[compute]) {
               acf[k] <- Re(sum(outer(seq(along = p), seq(along = p), function(i, j) 
                 p[i] * p[j] * exp(r[j] * tau * sr * abs(lag[k])) / (r[i] + r[j]) * 
                   (exp((r[i] + r[j]) * t[k] * tau * sr) - 1))))
             }
             
             acf / K2
           }
           
           # truncation
           if (is.null(len)) {
             len <- max(which(abs(ret$acfun(1:as.integer(10 / param$cutoff) / sr)) > 1e-3)) + 1L
             ret$len <- len
           }
           
           ret$truncatedKernfun <- function(t) {
             result <- numeric(length(t))
             compute <- t > 0 & t <= len / sr
             if (any(compute)) {
               result[compute] <- ret$kernfun(t[compute]) / ret$stepfun(len / sr)
             }
             result
           }
           ret$truncatedStepfun <- function(t) {
             result <- numeric(length(t))
             compute <- t > 0 & t < len / sr
             if (any(compute)) {
               result[compute] <- ret$stepfun(t[compute]) / ret$stepfun(len / sr)
             }
             result[t >= len / sr] <- 1
             result
           }
           K2truncated <- -a[1]^2 / ret$stepfun(len / sr)^2 * tau * sr * 
             Re(sum(outer(seq(along = p), seq(along = p), function(i, j) p[i] * p[j] / (r[i] + r[j]) * 
                            (exp((r[i] + r[j]) * len * tau) - 1))))
           ret$truncatedAcfun <- function(t) {
             acf <- numeric(length(t))
             compute <- t >= -len / sr & t <= len / sr
             
             if (any(compute)) {
               acf[compute] <- sapply(t[compute] * tau * sr,
                                      function(s) -a[1]^2 / ret$stepfun(len / sr)^2 * tau * sr *
                                        Re(sum(outer(seq(along = p), seq(along = p), function(i, j) 
                                          p[i] * p[j] * exp(r[j] * abs(s)) / (r[i] + r[j]) * (exp((r[i] + r[j]) *
                                            (len * tau - abs(s))) - 1))
                                        ))
               )
             }

             acf / K2truncated
           }
           ret$truncatedAcAntiderivative <- function(t, lag) {
             if (length(t) == 1L) {
               t <- rep(t, length(lag))
             }
             if (length(lag) == 1L) {
               lag <- rep(lag, length(t))
             }
             
             acf <- numeric(length(t))
             compute <- lag >= -len / sr & lag <= len / sr & t > 0
             trunc <- t[compute] > len / sr - abs(lag[compute])
             t[compute][trunc] <- len / sr - abs(lag[compute][trunc])
             
             for (k in seq(along = lag)[compute]) {
               acf[k] <- -a[1]^2 / ret$stepfun(len / sr)^2 * tau * sr *
                 Re(sum(outer(seq(along = p), seq(along = p), function(i, j) 
                   p[i] * p[j] * exp(r[j] * tau * sr * abs(lag[k])) / (r[i] + r[j]) * 
                     (exp((r[i] + r[j]) * abs(t[k] * tau * sr)) - 1))))
             }
             acf / K2truncated
           }
           
           
           shift <- (0:len + shift) / sr
           ret$kern <- ret$kernfun(shift) / sum(ret$kernfun(shift))
           ret$step <- ret$stepfun(shift)
           ret$acf <- ret$acfun(0:len / sr)
           
           indicesRightHalf <- which(ret$step >= 0.5)
           if (length(indicesRightHalf) > 0L) {
             ret$jump <- min(indicesRightHalf) - 1L # last index of left half
           } else {
             warning("jump could not be computed, since all step values are below 0.5, jump == len is returned")
             ret$jump <- ret$len
           }
           
           ret$number <- 0L
           ret$list <- list(truncation = as.numeric(ret$len / ret$sr),
                            C = a[1] / ret$stepfun(ret$len / ret$sr),
                            timescaling = ret$param$cutoff / omega0 * 2 * pi * ret$sr,
                            A = a[1] / ret$stepfun(ret$len / ret$sr) * 
                              (-1)^ret$param$pole * Re(1 / prod(r)),
                            a = Re(r),
                            b = Im(r),
                            c = Re(p),
                            d = Im(p))
         },
         stop("Unknown type")
  )
  
  class(ret) <- c("lowpassFilter", class(ret))
  ret
}

.BesselPolynomial <- function(n, reverse = FALSE) {
  k <- 0:n
  y.2 <- 1L
  y.1 <- c(1L, 1L)
  if (n == 0L) {
    y <- y.2
  } else if (n == 1L) {
    y <- y.1
  } else {
    for (i in 2:n) {
      y <- (2L * i - 1L) * c(0L, y.1) + c(y.2, 0L, 0L)
      y.2 <- y.1
      y.1 <- y
    }
  }
  if (reverse) rev(y) else y # if reverse return coefficients from highest to lowest
}

print.lowpassFilter <- function(x, ...) {
  cat("\n")
  switch(x$type,
         bessel = {
           cat(x$param$pole, "-pole Bessel filter\n\n", sep = "")
           cat("cut-off frequency:", x$param$cutoff, "\n")
         }
  )
  cat("length:", x$len, "\n")
  cat("sampling rate:", x$sr, "\n")
  cat("\n")
}
