# standard single-interval Simpson rule for numerical integration
# implemented as a reference
simpson.standard <- function(f, a, b, arg) {
  c <- (a + b) / 2.0
  h3 <- (b - a) / 6.0
  return (h3*(f(a, arg) + 4.0*f(c, arg) + f(b, arg)))
}

# numerical integration "\int_a^b f(t, ...) dg(t, ...)" over a single interval, using generalized Simpson's rule
simpson.generalized <- function(f, g = function(x, ...) x, a, b, ...) { # can be abbreviated if f,g are assumed to be vectorized
  m <- (a + b) / 2
  ga <- g(a, ...)
  gm <- g(m, ...)
  gb <- g(b, ...)
  fa <- f(a, ...)
  fm <- f(m, ...)
  fb <- f(b, ...)
  
  h <- gb - ga
  d <- 2 * gm - (ga + gb)
  r <- d / h
  
  ret <- (h / 6) * ((fa + 4 * fm + fb) + 2 * r * (fa - fb) - 3 * r^2 * (fa + fb)) / (1 - r^2)
  return (ret)
}

plot.summary.cfc <- function(
  x
  , which = c(1, 2)
  , ...
  ) {
  tvec <- x$tout
  ci.q <- x$ci
  s.q <- x$s
  quantiles <- x$quantiles

  dimCI <- dim(ci.q)
  K <- dimCI[length(dimCI)]
  
  if (1 %in% which) {
    # cumulative incidence function
    for (k in 1:K) {
      if (quantiles) {
        ci.ylim <- range(ci.q[, , k])
        plot(tvec, ci.q[2, , k], type = "l"
             , xlab = "time from index", ylab = "CI"
             , ylim = ci.ylim, main = paste0("cause ", k, " - cumulative incidence"))
        lines(tvec, ci.q[1, , k], lty = 2)
        lines(tvec, ci.q[3, , k], lty = 2)
      } else {
        plot(tvec, ci.q[, k], type = "l"
             , xlab = "time from index", ylab = "CI"
             , main = paste0("cause ", k, " - cumulative incidence"))
      }
    }
  }
  
  if (2 %in% which) {
    # survival probability
    for (k in 1:K) {
      if (quantiles) {
        s.ylim <- range(s.q[, , k])
        plot(tvec, s.q[2, , k], type = "l"
             , xlab = "time from index", ylab = "survival probability"
             , ylim = s.ylim, main = paste0("cause ", k, " - survival"))
        lines(tvec, s.q[1, , k], lty = 2)
        lines(tvec, s.q[3, , k], lty = 2)
      } else {
        plot(tvec, s.q[, k], type = "l"
             , xlab = "time from index", ylab = "survival probability"
             , main = paste0("cause ", k, " - survival"))
      }
    }
  }
}

summary.cfc <- function(
  object
  , f.reduce = function(x) x # transformation function applied to values at each time-point and cause
  , pval = 0.05
  , ... # other arguments passed to 'f'
) {
  # check if f.reduce output is a vector or a scalar
  # if vector: we treat it as samples and compute quantiles for it
  # if scalar: just report the single number; no quantile calculation needed
  isVector <- length(f.reduce(object$ci[1, 1, ], ...)) > 1
  
  s <- NA
  ci <- NA
  if (isVector) {
    ci <- apply(object$ci, MARGIN = c(1, 2), FUN = function(x) {
      quantile(f.reduce(x, ...), probs = c(pval / 2, 0.5, 1 - pval / 2))
    })
    s <- apply(object$s, MARGIN = c(1, 2), FUN = function(x) {
      quantile(f.reduce(x, ...), probs = c(pval / 2, 0.5, 1 - pval / 2))
    })
  } else {
    ci <- apply(object$ci, MARGIN = c(1, 2), FUN = f.reduce, ...)
    s <- apply(object$s, MARGIN = c(1, 2), FUN = f.reduce, ...)
  }
  
  ret <- list(tout = object$tout, ci = ci, s = s, quantiles = isVector)
  class(ret) <- c("summary.cfc", class(ret))

  return (ret)
}

print.summary.cfc <- function(x, ...) {}

cfc <- function(f.list, args.list, n, tout, Nmax = 100L, rel.tol = 1e-5, ncores = 1) {
  if (!is.list(f.list)) stop("f.list must be a list")
  K <- length(f.list) # number of causes
  if (K != length(args.list)) stop("length of args.list must match that of f.list")
  if (is.function(f.list[[1]])) { # R path
    ret <- cscr.samples.R(f.list, args.list, tout, Nmax, rel.tol, n, ncores)
  } else { # Cpp path
    func_list <- list()
    init_list <- list()
    free_list <- list()
    for (k in 1:K) {
      func_list[[k]] <- f.list[[k]][[1]]
      init_list[[k]] <- f.list[[k]][[2]]
      free_list[[k]] <- f.list[[k]][[3]]
    }
    ret <- cscr_samples_Cpp(func_list, init_list, free_list, args.list, tout, Nmax, rel.tol, n, ncores)
    ret$is.maxiter <- as.vector(ret$is.maxiter) # TODO: find better way to covert to vector (from matrix) inside Cpp code
  }
  ret$tout <- tout # we use this in plotting inside summary.cfc
  class(ret) <- c("cfc", class(ret))
  return (ret)
}

cscr.samples.R <- function(f.list, args.list, tout, Nmax = 100L, rel.tol = 1e-6, nsmp, ncores = 1) {
  trapezoidal.step <- function(fvec, gvec) {
    h <- gvec[3] - gvec[1]
    ret <- (h / 2) * (fvec[1] + fvec[3])
    return (ret)
  }
  simpson.step <- function(fvec, gvec) {
    gdiff <- diff(gvec)
    #h <- gvec[3] - gvec[1]
    #if (abs(h) < .Machine$double.eps) return (0.0)
    if (abs(gdiff[1]) < .Machine$double.eps) {
      return (trapezoidal.step(c(fvec[2], NA, fvec[3]), c(gvec[2], NA, gvec[3])))
    } else if (abs(gdiff[2]) < .Machine$double.eps) {
      return (trapezoidal.step(c(fvec[1], NA, fvec[2]), c(gvec[1], NA, gvec[2])))
    }
    #r <- (2* gvec[2] - gvec[1] - gvec[3]) / h
    h <- gdiff[1] + gdiff[2]
    r <- (gdiff[1] - gdiff[2]) / h
    ret <- (h / 6) * ((fvec[1] + 4 * fvec[2] + fvec[3]) + 2 * r * (fvec[1] - fvec[3]) - 3 * r^2 * (fvec[1] + fvec[3])) / (1 - r^2)
    return (ret)
  }
  
  tmax <- max(tout)
  nout <- length(tout)
  K <- length(f.list)
  I_out_value <- array(NA, dim = c(nout, K, nsmp))
  scube <- array(NA, dim = c(nout, K, nsmp)) # unadjusted survival probabilities
  smat <- array(NA, dim = c(nout, K))
  
  is.maxiter <- rep(F, nsmp)
  registerDoParallel(ncores)
  retsink <- foreach(j=1:nsmp) %dopar% {
  #for (j in 1:nsmp) {
    N <- 1
    
    xvec.new <- c(0.0, tmax/2, tmax)
    fmat.new <- array(NA, dim = c(3, K))
    for (k in 1:K) fmat.new[, k] <- f.list[[k]](xvec.new, args.list[[k]], j)
    fprodmat.new <- t(apply(fmat.new, 1, function(x) prod(x)/x)) # TODO: handle division by zero
    
    I.trap.int.new <- -sapply(1:K, function(k) {
      trapezoidal.step(fprodmat.new[, k], fmat.new[, k])
    })
    
    I.simp.int.new <- -sapply(1:K, function(k) {
      simpson.step(fprodmat.new[, k], fmat.new[, k])
    })
    
    xvec <- xvec.new
    f.mat <- fmat.new
    fprod.mat <- fprodmat.new
    I.trap.int <- matrix(I.trap.int.new, ncol = K)
    I.simp.int <- matrix(I.simp.int.new, ncol = K)
    
    I.trap.cum <- rbind(0, I.trap.int)
    I.simp.cum <- rbind(0, I.simp.int)
    
    err.abs.int <- abs(I.simp.int - I.trap.int)
    err.abs.cum <- abs(I.simp.cum - I.trap.cum)
    err.rel.cum <- err.abs.cum / abs(I.simp.cum)
    err.abs <- err.abs.cum[N + 1, ]
    err.rel <- err.rel.cum[N + 1, ]
    err.rel.max <- max(err.rel)
    
    while(max(err.rel) > rel.tol && N < Nmax) { # add Nmin
      idx <- which.max(apply(err.abs.int, 1, max))
      
      xvec.new <- c(mean(xvec[(2*idx - 1):(2*idx)]), xvec[2*idx], mean(xvec[(2*idx):(2*idx + 1)]))
      xvec <- c(xvec[1:(2*idx - 1)], xvec.new, xvec[(2*idx + 1):(2*N + 1)])
      
      for (k in 1:K) {
        fmat.new[, k] <- f.list[[k]](xvec.new, args.list[[k]], j)
      }
      fprodmat.new <- t(apply(fmat.new, 1, function(x) prod(x)/x)) # TODO: handle division by zero
      
      f.mat <- rbind(f.mat[1:(2*idx - 1), ], fmat.new, f.mat[(2*idx + 1):(2*N + 1), ])
      fprod.mat <- rbind(fprod.mat[1:(2*idx - 1), ], fprodmat.new, fprod.mat[(2*idx + 1):(2*N + 1), ])
      
      I.trap.int.1 <- -sapply(1:K, function(k) {
        trapezoidal.step(fprod.mat[2*idx + (-1:1), k], f.mat[2*idx + (-1:1), k])
      })
      I.trap.int.2 <- -sapply(1:K, function(k) {
        trapezoidal.step(fprod.mat[2*idx + (1:3), k], f.mat[2*idx + (1:3), k])
      })
      I.simp.int.1 <- -sapply(1:K, function(k) {
        simpson.step(fprod.mat[2*idx + (-1:1), k], f.mat[2*idx + (-1:1), k])
      })
      I.simp.int.2 <- -sapply(1:K, function(k) {
        simpson.step(fprod.mat[2*idx + (1:3), k], f.mat[2*idx + (1:3), k])
      })
      
      if (idx == 1) {
        I.trap.int <- rbind(I.trap.int.1, I.trap.int.2, I.trap.int[-1, ])
        I.simp.int <- rbind(I.simp.int.1, I.simp.int.2, I.simp.int[-1, ])
      } else if (idx == N) {
        I.trap.int <- rbind(I.trap.int[-N, ], I.trap.int.1, I.trap.int.2)
        I.simp.int <- rbind(I.simp.int[-N, ], I.simp.int.1, I.simp.int.2)
      } else {
        I.trap.int <- rbind(I.trap.int[1:(idx - 1), ], I.trap.int.1, I.trap.int.2, I.trap.int[(idx + 1):N, ])
        I.simp.int <- rbind(I.simp.int[1:(idx - 1), ], I.simp.int.1, I.simp.int.2, I.simp.int[(idx + 1):N, ])
      }
      I.trap.cum <- rbind(0, apply(I.trap.int, 2, cumsum)) # room for efficiency
      I.simp.cum <- rbind(0, apply(I.simp.int, 2, cumsum)) # room for efficiency
      
      err.abs.cum <- abs(I.simp.cum - I.trap.cum)
      err.rel.cum <- err.abs.cum / abs(I.simp.cum)
      err.abs.int <- abs(I.simp.int - I.trap.int)
      err.rel.int <- err.abs.int / abs(I.simp.int)
      
      err.abs <- err.abs.cum[N + 2, ]
      err.rel <- err.rel.cum[N + 2, ]
      err.rel.max <- max(err.rel)
      
      N <- N + 1
    }
    
    idx.nodes <- seq(from = 1, to = 2*N + 1, by = 2)
    I.trap.out <- apply(I.trap.cum, 2, function(x) approx(xvec[idx.nodes], x, tout)$y) # TODO: consider nonlinear interpolations to make the curves look smoother
    I.simp.out <- apply(I.simp.cum, 2, function(x) approx(xvec[idx.nodes], x, tout)$y) # TODO: see above
    
    for (k in 1:K) {
      smat[, k] <- f.list[[k]](tout, args.list[[k]], j)
    }
    
    return (list(N = N, I.simp.out = I.simp.out, smat = smat))
  }
  stopImplicitCluster()

  for (j in 1:nsmp) {
    is.maxiter[j] <- 1*(retsink[[j]]$N == Nmax)
    I_out_value[, , j] <- retsink[[j]]$I.simp.out
    scube[, , j] <- retsink[[j]]$smat
  }

  n.maxiter <- sum(is.maxiter)

  if (n.maxiter > 0) warning(paste0(n.maxiter, " of ", nsmp, " integrals did not converge after reaching maximum iterations"))
  
  return (list(ci = I_out_value, s = scube, is.maxiter = is.maxiter, n.maxiter = n.maxiter))
}
