context("Testing polynomial example")

# simulate the data
set.seed(1)
n_global <- 3L
n_clusters <- 10L

mu_global <- rnorm(n_global)
idx_start <- n_global

# cluster_dat <- replicate(n_clusters, {
#   n_members <- sample.int(n_global, 1L)
#   g_idx <- sort(sample.int(n_global, n_members))
#   mu_cluster <- rnorm(n_members)
#   Psi <- matrix(rnorm(n_members * n_members), n_members, n_members)
#
#   out <- list(idx = idx_start + 1:n_members, g_idx = g_idx,
#               mu_cluster = mu_cluster, Psi = Psi)
#   idx_start <<- idx_start + n_members
#   out
# }, simplify = FALSE)
#
# saveRDS(cluster_dat, file.path("tests", "testthat", "poly.RDS"))
f <- file.path("tests", "testthat", "poly.RDS")
if(!file.exists(f))
  f <- "poly.RDS"
cluster_dat <- readRDS(f)

test_that("Poly example gives the same", {
  skip_if_not_installed("Rcpp")
  skip_if_not_installed("RcppArmadillo")
  skip_on_macOS()

  Rcpp::sourceCpp(system.file("poly-ex.cpp", package = "psqn"))

  optimizer <- get_poly_optimizer(
    cluster_dat, max_threads = 2L, mu_global = mu_global)

  val <- c(1.77, 0.72, 0.91, 0.38, 1.68, -0.64, -0.46, 1.43, -0.65, -0.21,
           -0.39, -0.32, -0.28, 0.49, -0.18, -0.51, 1.34, -0.21, -0.18,
           -0.1, 0.71, -0.07)
  fv <- eval_poly(val = val, ptr = optimizer, n_threads = 1L)
  expect_equal(fv, 70.7249805284658)
  fv <- eval_poly(val = val, ptr = optimizer, n_threads = 2L)
  expect_equal(fv, 70.7249805284658)

  gr_res <- structure(
    c(40.002794468683, -4.69365742157453, 27.8942951540135,
      -2.24607774807025, 4.44258239006445, -4.88972063019875, -1.42070419832834,
      7.52804215023212, -5.41336909997731, -5.555925884566, -1.95691439753132,
      -0.259270763193206, -0.567523464319468, 0.990177191416215, 1.42494449144973,
      -1.94231181231703, 4.80032005237134, -0.672040114327287, -2.50984701206127,
      2.04148823243275, 4.21292066165825, -4.04007073762148),
    value = 70.7249805284658)
  gr <- grad_poly(val = val, ptr = optimizer, n_threads = 1L)
  expect_equal(gr, gr_res)
  gr <- grad_poly(val = val, ptr = optimizer, n_threads = 2L)
  expect_equal(gr, gr_res)

  rel_eps <- sqrt(.Machine$double.eps)
  opt <- optim_poly(
    val = val, ptr = optimizer, rel_eps = rel_eps, max_it = 100L,
    c1 = 1e-4, c2 = .9, n_threads = 1L)
  opt_res <- list(par = c(-0.626453810632069, 0.183643324345628, -0.835628612250485,
                          -2.04960116005313, 0.602201895204366, -0.601922171312659, 0.936438001760361,
                          1.7873595917768, -1.01160126624839, -0.125807962111816, 1.30646428460648,
                          0.823781885813167, 0.357668134623467, 0.196958422127988, -0.760960914726239,
                          0.762492431522127, 0.866996127883627, -0.234758322046668, 0.307564405437029,
                          1.43619906759021, 1.32150098518575, -0.873195703378549), value = 2.28888765400896e-19,
                  info = 0L, counts = c(`function` = 36, gradient = 35, n_cg = 244
                  ), convergence = TRUE)
  tol <- sqrt(rel_eps) * 10
  do_check <- !names(opt) %in% "counts"
  expect_equal(opt[do_check], opt_res[do_check], tolerance = tol)
  opt <- optim_poly(
    val = val, ptr = optimizer, rel_eps = rel_eps, max_it = 100L,
    c1 = 1e-4, c2 = .9, n_threads = 2L)
  expect_equal(opt[do_check], opt_res[do_check], tolerance = tol)

  opt <- optim_poly(
    val = val, ptr = optimizer, rel_eps = rel_eps, max_it = 100L,
    c1 = 1e-4, c2 = .9, n_threads = 1L, strong_wolfe = FALSE)
  opt_res <- list(par = c(-0.626453810632069, 0.183643324345628, -0.835628612250485,
                          -2.04960116005313, 0.602201895204366, -0.601922171312659, 0.936438001760361,
                          1.7873595917768, -1.01160126624839, -0.125807962111816, 1.30646428460648,
                          0.823781885813167, 0.357668134623467, 0.196958422127988, -0.760960914726239,
                          0.762492431522127, 0.866996127883627, -0.234758322046668, 0.307564405437029,
                          1.43619906759021, 1.32150098518575, -0.873195703378549), value = 2.28888765400896e-19,
                  info = 0L, counts = c(`function` = 36, gradient = 35, n_cg = 244
                  ), convergence = TRUE)
  expect_equal(opt[do_check], opt_res[do_check], tolerance = tol)
})
