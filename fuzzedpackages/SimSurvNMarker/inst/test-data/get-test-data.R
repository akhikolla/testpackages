#####
# input
d_m <- 2L
d_g <- 3L
d_b <- 2L
n_y <- 2L
d_z <- 4L
d_x <- 3L
K <- n_y * d_m

#####
# output
b_ks <- seq(log(1e-1), log(10), length.out = d_b)
m_ks <- seq(0     , 10     , length.out = d_m)
g_ks <- seq(0     , 10     , length.out = d_g)

if(d_b > 0){
  omega <- runif(d_b, -1, 1)
  omega <- omega + local({
    library(splines)
    Boundary_knots <- b_ks[c(1, d_b)]
    a_knots <- sort(c(rep(Boundary_knots, 4), b_ks[-c(1, d_b)]))

    points <- b_ks
    basis <- splineDesign(a_knots, points, ord = 4L, outer.ok = TRUE)
    const <- splineDesign(a_knots, Boundary_knots, ord = 4L,
                          derivs = c(2L, 2L))
    qr.const <- qr(t(const))
    basis <- (t(qr.qty(qr.const, t(basis))))[, -(1L:2L), drop = FALSE]

    omega_const <- solve(basis, rep(-2, d_b))
    drop(omega + omega_const)
  })
} else
  omega <- numeric()

dput(omega <- round(omega, 2))
dput(B <- matrix(round(runif(n_y * d_g, -1, 1), 2), nrow = d_g))
dput(Psi <- round(drop(rWishart(1L, 2 * K, diag(K) / K)), 2))
dput(sig <- diag(round(runif(n_y, .1, .4)^2, 2), n_y))
dput(delta <- round(runif(d_z, -1, 1) / d_z, 2))
dput(alpha <- round(runif(n_y, -1, 1) / n_y, 2))
dput(gamma <- matrix(round(runif(d_x * n_y, -1, 1) / d_x, 2),
                     nc = n_y))
n_obs <- 2000L
