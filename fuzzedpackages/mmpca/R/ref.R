# reference implementations for functions implemented in C++

# small constant used in soft abs and soft sign
eps <- 1e-6

# g differentiated with respect to (the vector) D1
ref_dlossdD1 <- function(xi1, D1, xi2, D2, X, mask) {
  V1 <- ref_Vxi(xi1)
  V2 <- ref_Vxi(xi2)
  diag(-2*t(V1) %*% (mask * (X - V1 %*% diag(D1 * D2) %*% t(V2))) %*% V2 %*%
    diag(D2))
}

# g differentiated with respect to (the matrix) xi1
ref_dlossdxi1 <- function(xi1, D1, xi2, D2, X, mask) {
  n <- nrow(xi1)
  k <- ncol(xi1)
  T2 <- ref_Vxi(xi2) %*% diag(D2 * D1)
  T1 <- t((mask * (X-ref_Vxi(xi1)%*%t(T2)))%*%T2)
  dgdxi1 <- matrix(0, n, k)
  V <- diag(n)
  for (j in 1:min(k, n-1)) {
    for (i in (j+1):n) {
      dgdR <- -2*ref_VpartR(i, j, xi1) %*% T1 %*% V
      dgdxi1[i, j] <- sum(t(dgdR) * ref_dR(i, j, n, xi1[i, j]))
      V <- V %*% ref_R(i, j, n, xi1[i, j])
    }
  }
  return(dgdxi1)
}

# R differentiated with respect to xi
ref_dR <- function(i, j, n, xi) {
  dR <- matrix(0, n, n)
  dR[i, i] <- dR[j, j] <- -sin(xi)
  dR[j, i] <- cos(xi)
  dR[i, j] <- -dR[j, i]
  return(dR)
}

ref_gradient <- function(x, X, masks, inds, k, p, lambda) {
  n <- max(inds)
  tmp <- ref_unvectorize(x, k, n, p)
  xi <- tmp$xi
  D <- tmp$D
  dxi <- list()
  dD <- matrix(0, k, n)
  for (i in 1:n) {
    dxi[[i]] <- matrix(0, nrow(xi[[i]]), ncol(xi[[i]]))
  }
  for (i in 1:length(X)) {
    row <- inds[i, 1]
    col <- inds[i, 2]
    dxi[[row]] <- dxi[[row]] +
      ref_dlossdxi1(xi[[row]], D[, row], xi[[col]], D[, col], X[[i]],
      masks[[i]])
    dxi[[col]] <- dxi[[col]] +
      ref_dlossdxi1(xi[[col]], D[, col], xi[[row]], D[, row], t(X[[i]]),
      t(masks[[i]]))
    dD[, row] <- dD[, row] +
      ref_dlossdD1(xi[[row]], D[, row], xi[[col]], D[, col], X[[i]], masks[[i]])
    dD[, col] <- dD[, col] +
      ref_dlossdD1(xi[[col]], D[, col], xi[[row]], D[, row], t(X[[i]]),
      t(masks[[i]]))
  }
  for (view in 1:n) {
    V <- ref_Vxi(xi[[view]])
    vd <- V %*% diag(D[, view])
    for (j in 1:min(k, p[view]-1)) {
      for (i in (j+1):p[view]) {
        A <- ref_VpartL(i, j, xi[[view]])
        B <- ref_VpartR(i, j, xi[[view]])
        dxi[[view]][i, j] <- dxi[[view]][i, j] + lambda[3] * sum(
          t(A)%*%ref_soft_sign(vd)%*%diag(D[, view])%*%t(B)*
          ref_dR(i, j, p[view], xi[[view]][i, j]))
        for (row in 1:p[view]) {
          denom <- sqrt(sum(vd[row, ]^2))
          if (denom > 1e-8) {
            L <- matrix(0, 1, p[view])
            L[1, row] <- 1
            dxi[[view]][i, j] <- dxi[[view]][i, j] + lambda[4] * sum(
              t(A)%*%t(L)%*%L%*%V%*%diag(D[, view]^2)%*%t(B)*
              ref_dR(i, j, p[view], xi[[view]][i, j])
              ) / denom
          }
        }
      }
    }
  }
  dD <- dD + matrix(ref_penalty_gradient(D, xi) %*% lambda, k, n)
  return(ref_vectorize(dxi, dD))
}

ref_penalty_gradient <- function(D, xi) {
  a <- sqrt(rowSums(D^2))
  find_k_penalty <- D / (a %*% matrix(1, 1, ncol(D)))
  find_k_penalty[a < 1e-8, ] <- 0
  sparsity_penalty <- matrix(NA, nrow(D), ncol(D))
  varsel_penalty <- matrix(0, nrow(D), ncol(D))
  for (i in 1:length(xi)) {
    V <- ref_Vxi(xi[[i]])
    vd <- V %*% diag(D[, i])
    sparsity_penalty[, i] <- diag(t(V) %*% ref_soft_sign(vd))
    denom <- sqrt(rowSums(vd^2))
    for (j in 1:nrow(V)) {
      if (denom[j] > 1e-8) {
        varsel_penalty[, i] <- varsel_penalty[, i] + (V[j, ] * vd[j, ])/denom[j]
      }
    }
  }
  cbind(c(ref_soft_sign(D)), c(find_k_penalty), c(sparsity_penalty),
    c(varsel_penalty))
}

ref_objective <- function(x, X, masks, inds, k, p, lambda) {
  n <- max(inds)
  tmp <- ref_unvectorize(x, k, n, p)
  xi <- tmp$xi
  V <- lapply(xi, ref_Vxi)
  D <- tmp$D
  loss <- 0
  for (i in 1:length(X)) {
    row <- inds[i, 1]
    col <- inds[i, 2]
    loss <- loss +
      sum(masks[[i]] *
        (X[[i]] - V[[row]] %*% diag(D[, row] * D[, col]) %*% t(V[[col]]))^2)
  }
  penalties <- rep(NA, 3)
  penalties[1] <- sum(ref_soft_abs(D)) # integration penalty
  penalties[2] <- sum(apply(D, 1, function(x) sqrt(sum(x^2)))) # rank penalty
  penalties[3] <- sum(sapply(1:n,
    function(i) sum(abs(V[[i]] %*% diag(D[, i])))))# sparsity penalty
  penalties[4] <- 0 # variable selection penalty
  for (i in 1:n) {
    vd <- V[[i]] %*% diag(D[, i])
    for (j in 1:p[i]) {
      penalties[4] <- penalties[4] + sqrt(sum(vd[j, ]^2))
    }
  }
  return(loss + sum(lambda * penalties))
}

ref_optim_mmpca <- function(start, x, inds, k, p, lambda, trace) {
  f <- function(theta) ref_objective(theta, x, inds, k, p, lambda)
  df <- function(theta) ref_gradient(theta, x, inds, k, p, lambda)
  sol <- stats::optim(start, f, df, method='BFGS',
    control=list(maxit=1e6, reltol=.Machine$double.eps))
  if (trace) {
    message(sol[-1])
  }
  list(sol$par, NA, NA, NA, NA, NA)
}

# a rotation in a 2D subspace
ref_R <- function(i, j, n, xi) {
  R <- diag(n)
  R[i, i] <- R[j, j] <- cos(xi)
  R[j, i] <- sin(xi)
  R[i, j] <- -R[j, i]
  return(R)
}

ref_soft_abs <- function(x) {
  eps * gsl::lncosh(1/eps * x)
}

ref_soft_sign <- function(x) {
  tanh(1/eps * x)
}

ref_vectorize <- function(xi, D) {
  do.call(base::c, lapply(c(xi, D), function(x) x[]))
}

# the part of V (R_1...R_a) to the left of R_{ii,jj}
ref_VpartL <- function(ii, jj, xi) {
  n <- nrow(xi)
  k <- ncol(xi)
  V <- diag(n)
  for (j in k:1) {
    for (i in n:(j+1)) {
      if (j < jj || (j == jj && i < ii)) {
        V <- ref_R(i, j, n, xi[i, j]) %*% V
      }
    }
  }
  return(V)
}

# the part of V (R_a...R_mI_{nk}) to the right of R_{ii,jj}
ref_VpartR <- function(ii, jj, xi) {
  n <- nrow(xi)
  k <- ncol(xi)
  V <- matrix(0, n, k)
  diag(V) <- 1
  for (j in min(k, n-1):1) {
    for (i in n:(j+1)) {
      if (j > jj || (j == jj && i > ii)) {
        V <- ref_R(i, j, n, xi[i, j]) %*% V
      }
    }
  }
  return(V)
}

# the rotation matrix for a set of angles
# xi is lower triangular
ref_Vxi <- function(xi) {
  n <- nrow(xi)
  if (n == 1) {
    return(matrix(c(1, rep(0, ncol(xi)-1)), 1, ncol(xi)))
  }
  k <- min(n, ncol(xi))
  V <- diag(n)[, 1:k, drop=FALSE]
  for (j in min(k, n-1):1) {
    for (i in n:(j+1)) {
      V <- ref_R(i, j, n, xi[i, j]) %*% V
    }
  }
  return(cbind(V, matrix(0, n, max(0, ncol(xi)-n))))
}

ref_unvectorize <- function(x, k, n, p) {
  i <- 1
  xi <- list()
  for (view in 1:n) {
    xi[[view]] <- matrix(x[i:(i+k*p[view]-1)], p[view], k)
    i <- i + k*p[view]
  }
  D <- matrix(x[i:(i+k*n-1)], k, n)
  return(list(xi=xi, D=D))
}

