# functions only used to find start values

init_d <- function(d) {
  # d[, , k] is the square upper-triangular matrix of d_ij for component k with
  #   NA for missing values
  #
  # Possible speed improvements:
  #   make m a sparse matrix and move it inside the parentheses
  #   make d0 a sparse matrix
  #   replace abs and sign with smooth approximations
  #   the loop is trivially parallelizable
  #
  # Invariants:
  #   should work for any positive semidefinite censored rank 1 positive
  #   matrices the loss sum(m*(d0 - x %*% t(x))^2) should be approximately zero
  #   at optimum
  #
  K <- dim(d)[3]
  N <- dim(d)[1]
  return_d <- matrix(NA, K, N)
  lambda <- 1e-6
  for (k in 1:K) {
    d0 <- d[, , k]
    if (all(is.na(d0))) {
      return_d[k, ] <- rep(0, N)
      next
    }
    signs <- init_signs(d0)
    m <- !is.na(d0)
    d0[is.na(d0)] <- 0
    f <- function(x) {
      sum(m*(d0 - (signs * x) %*% t(signs * x))^2) + lambda*sum(abs(x))
    }
    df <- function(x) {
      t0 <- m * (d0 - (signs * x) %*% t(signs * x))
      -2*(t0 %*% (signs * x) + t(t0) %*% (signs * x))*signs + lambda*sign(x)
    }
    startx <- sqrt(diag(d[, , k]))
    startx[is.na(startx)] <- 0.1
    sol <- stats::optim(startx, f, df, method='L-BFGS-B', lower=0,
      control=list(maxit=10000, factr=1e2))
    if (sol$convergence != 0) {
      if (!(sol$value < 1e-5 && sol$convergence == 52)) {
        message(d0)
        message(list(k=k))
        message(sol)
        stop('No convergence!')
      }
    }
    return_d[k, ] <- signs * sol$par
  }
  return(return_d)
}

init_signs <- function(d) {
  n <- ncol(d)
  signs <- rep(NA, n)
  d[lower.tri(d)] <- t(d)[lower.tri(d)]
  while (any(is.na(signs))) {
    signs[which(is.na(signs))[1]] <- 1
    while (TRUE) {
      nacount <- sum(is.na(signs))
      for (i in 1:n) {
        if (!is.na(signs[i])) {
          signs[!is.na(d[i, ])] <- sign(signs[i] / d[i, !is.na(d[i, ])])
        }
      }
      if (nacount == sum(is.na(signs))) {
        break
      }
    }
  }
  return(signs)
}

init_singular_values <- function(V, X, inds) {
  K <- ncol(V[[1]])
  N <- length(V)
  M <- length(X)
  d <- array(NA, c(N, N, K))
  for (m in 1:M) {
    i <- inds[m, 1]
    j <- inds[m, 2]
    d[i, j, ] <- diag(t(V[[i]]) %*% X[[m]] %*% V[[j]])
  }
  return(d)
}

init_v <- function(X, inds, K) {
  N <- max(inds)
  V <- list()
  transpose_all <- function(X) lapply(X, function(x) t(x))
  for (i in 1:N) {
    x <- do.call(cbind, c(X[inds[, 1] == i], transpose_all(X[inds[, 2] == i])))
    if (nrow(x) < K) {
      V[[i]] <- cbind(svd(x)$u, matrix(0, nrow(x), K-nrow(x)))
    } else {
      V[[i]] <- svd(x)$u[, 1:K]
    }
  }
  return(V)
}

# find dimension for each view
init_view_dimensions <- function(X, inds) {
  n <- max(inds)
  p <- rep(NA, n)
  for (i in 1:n) {
    if (i %in% inds[, 1]) {
      p[i] <- nrow(X[[which(i == inds[, 1])[1]]])
    } else {
      p[i] <- ncol(X[[which(i == inds[, 2])[1]]])
    }
  }
  return(p)
}

# ref: David K. Hoffman, Richard C. Raffenetti, and Klaus Ruedenberg,
# Generalization of Euler Angles to N-Dimensional Orthogonal Matrices (1972)
init_inv_v <- function(V) {
  orth.compl <- function(x) {
    qr.Q(qr(x), complete=TRUE)[, (ncol(x)+1):nrow(x), drop=FALSE]
  }
  ref_invVinner <- function(V, Vorth) {
    n <- nrow(V)
    xi <- matrix(0, n, n)
    t <- t(cbind(V, Vorth))
    for (v in n:2) {
      xi[v, v-1] <- atan2(t[v-1, v], t[v, v])
      if (v > 2) {
        for (k in (v-2):1) {
          s <- sin(xi[v, k+1])
          if (abs(t[k+1, v] - s) < 1e-17) {
            if (v == n && abs(s) < 1e-17) {
              xi[v, k] <- asin(t[k, v])
            } else {
              xi[v, k] <- 0
            }
          } else {
            xi[v, k] <- atan2(t[k, v], t[k+1, v]/s)
          }
          if (xi[v, k] > pi/2) xi[v, k] <- xi[v, k] - pi
          if (xi[v, k] < -pi/2) xi[v, k] <- xi[v, k] + pi
        }
      }
      tnew <- diag(n)
      f <- t
      for (k in (v-1):1) {
        tnew[k, ] <- cos(xi[v, k]) * t[k, ] - sin(xi[v, k]) * f[k+1, ]
        f[k, ] <- sin(xi[v, k]) * t[k, ] + cos(xi[v, k]) * f[k+1, ]
      }
      t <- tnew
    }
    return(-xi)
  }
  if (ncol(V) < nrow(V)) {
    Vorth <- orth.compl(V)
    xi <- c_invVinner(V, Vorth)[, 1:ncol(V), drop=FALSE]
    s <- sign(V[1, 1])
    if (s == 0) {
      return(xi)
    }
    shat <- sign(ref_Vxi(xi)[1, 1])
    if (shat == 0 || s == shat) {
      return(xi)
    }
    Vorth[, ncol(Vorth)] <- -Vorth[, ncol(Vorth)]
    xi <- c_invVinner(V, Vorth)[, 1:ncol(V), drop=FALSE]
    return(xi)
  } else if (ncol(V) > nrow(V)) {
    xi <- c_invVinner(V[, 1:nrow(V)], matrix(NA, nrow(V), 0))
    return(cbind(xi, matrix(0, nrow(xi), ncol(V)-nrow(V))))
  } else {
    return(c_invVinner(V, matrix(NA, nrow(V), 0)))
  }
}
