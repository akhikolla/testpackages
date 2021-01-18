args_env <- new.env()
local(envir = args_env, {
  d_m <- 4L
  d_g <- 5L
  d_b <- 5L
  n_y <- 1L
  d_z <- 2L
  d_x <- 2L
  K <- n_y * d_m

  b_ks <- seq(log(1e-2), log(10), length.out = d_b)
  m_ks <- seq(0        , 10     , length.out = d_m)
  g_ks <- seq(0        , 10     , length.out = d_g)

  omega <- c(-0.84, -2.67, -0.37, -4.62, 0.63)
  B <- structure(c(0.26, 0.3, 0.48, 0.8, 0.83), .Dim = c(5L, 1L))
  Psi <- structure(c(1.02, 0.18, -1.27, -0.39, 0.18, 1.33, 0.47, -0.74,
                     -1.27, 0.47, 3.63, 0.2, -0.39, -0.74, 0.2, 1.56),
                   .Dim = c(4L, 4L))
  sig <- structure(0.02, .Dim = c(1L, 1L))
  delta <- c(-0.1, -0.3)
  alpha <- -0.62
  gamma <- structure(c(0.11, 0.33), .Dim = 2:1)

  n_obs <- 2000L
})
