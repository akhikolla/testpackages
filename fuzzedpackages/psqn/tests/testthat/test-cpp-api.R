context("Testing C++ headers from R")

# assign model parameters and number of random effects and fixed effects
q <- 4
p <- 5
beta <- sqrt((1:p) / sum(1:p))
Sigma <- diag(q)

# # simulate data set
# n_clusters <- 4L
# set.seed(1)
#
# sim_dat <- replicate(n_clusters, {
#   n_members <- sample.int(20L, 1L) + 2L
#   X <- matrix(runif(p * n_members, -sqrt(6 / 2), sqrt(6 / 2)),
#               p)
#   u <- drop(rnorm(q) %*% chol(Sigma))
#   Z <- matrix(runif(q * n_members, -sqrt(6 / 2 / q), sqrt(6 / 2 / q)),
#               q)
#   eta <- drop(beta %*% X + u %*% Z)
#   y <- as.numeric((1 + exp(-eta))^(-1) > runif(n_members))
#
#   list(X = X, Z = Z, y = y, u = u, Sigma_inv = solve(Sigma))
# }, simplify = FALSE)
# saveRDS(sim_dat, file.path("tests", "testthat", "mixed-logit.RDS"))
f <- file.path("tests", "testthat", "mixed-logit.RDS")
if(!file.exists(f))
  f <- "mixed-logit.RDS"
sim_dat <- readRDS(f)

test_that("mixed logit model gives the same", {
  skip_if_not_installed("Rcpp")
  skip_if_not_installed("RcppArmadillo")
  skip_on_macOS()

  Rcpp::sourceCpp(system.file("mlogit-ex.cpp", package = "psqn"))
  optimizer <- get_mlogit_optimizer(sim_dat, max_threads = 2L)

  val <- c(beta, sapply(sim_dat, function(x) x$u))
  fv <- eval_mlogit(val = val, ptr = optimizer, n_threads = 1L)
  expect_equal(fv, 35.248256794375)
  fv <- eval_mlogit(val = val, ptr = optimizer, n_threads = 2L)
  expect_equal(fv, 35.248256794375)

  gr_res <- structure(
    c(-3.09880891801443, -3.11919456739973, -1.26969269606313,
      5.50412085634227, -5.07396341960751, 0.407174373086653, 0.590879440022872,
      0.122784612825495, 0.222039129397368, -0.335740954704027, 0.967036607234169,
      -2.05362226641706, -0.517311694420387, -2.33859175587387, 0.911175257815629,
      -1.80538928384093, -0.81631781729293, 3.14749812772018, -0.894143195465396,
      0.647354085121652, 0.344172484580979), value = 35.248256794375)
  gr <- grad_mlogit(val = val, ptr = optimizer, n_threads = 1L)
  expect_equal(gr, gr_res)
  gr <- grad_mlogit(val = val, ptr = optimizer, n_threads = 2L)
  expect_equal(gr, gr_res)

  rel_eps <- sqrt(.Machine$double.eps)
  opt <- optim_mlogit(
    val = val, ptr = optimizer, rel_eps = rel_eps, max_it = 100L,
    c1 = 1e-4, c2 = .9, n_threads = 1L)
  opt_res <- list(par = c(0.647310437469711, 0.773205671911178, 0.656983834118283,
                          0.119681528415682, 0.855849653168437, -0.484311298187747, 0.412484752266496,
                          0.766994317977443, 0.422874453513229, 0.605169744625452, -0.783701817869368,
                          0.096194899967712, -0.922934757392211, -0.483542761956257, 0.438193243670319,
                          -0.493309514536335, 0.0274193762771719, 0.426816737611666, 0.0281704108821127,
                          -0.180877909582702, -0.150512019701299), value = 25.1215015278071,
                  info = 0L, counts = c(`function` = 12, gradient = 11, n_cg = 40
                  ), convergence = TRUE)
  tol <- sqrt(rel_eps) * 10
  do_check <- !names(opt) %in% "counts"
  expect_equal(opt[do_check], opt_res[do_check], tolerance = tol)
  opt <- optim_mlogit(
    val = val, ptr = optimizer, rel_eps = rel_eps, max_it = 100L,
    c1 = 1e-4, c2 = .9, n_threads = 2L)
  expect_equal(opt[do_check], opt_res[do_check], tolerance = tol)

  # check the function to optimize the private parameters
  start_priv <- opt$par
  start_priv[-seq_along(beta)] <- 0
  out <- optim_mlogit_private(
    val = start_priv, ptr = optimizer, rel_eps = rel_eps, max_it = 100L,
    c1 = 1e-4, c2 = .9, n_threads = 2L)
  expect_equal(opt$par, out, tolerance = tol, check.attributes = FALSE)
  expect_true(all(opt$par[seq_along(beta)] == out[seq_along(beta)]))
  expect_equal(opt$value, attr(out, "value"), tolerance = rel_eps)

  opt <- optim_mlogit(
    val = val, ptr = optimizer, rel_eps = rel_eps, max_it = 100L,
    c1 = 1e-4, c2 = .9, use_bfgs = FALSE,
    n_threads = 1L)
  opt_res <- list(par = c(0.647309984423475, 0.773099829472392, 0.657030091436329,
                          0.119542556826119, 0.856048420971233, -0.48432079434579, 0.412471124068345,
                          0.76694534092756, 0.422786612641562, 0.605124122567012, -0.783742636851702,
                          0.0960937729102007, -0.922965266682886, -0.483501533118553, 0.438327834070562,
                          -0.493318173514831, 0.0274749226734802, 0.426944683060686, 0.0287992526677307,
                          -0.180834889398274, -0.150315014163914), value = 25.121501682662,
                  info = 0L, counts = c(`function` = 14, gradient = 11, n_cg = 22
                  ), convergence = TRUE)
  expect_equal(opt[do_check], opt_res[do_check], tolerance = tol)
  opt <- optim_mlogit(
    val = val, ptr = optimizer, rel_eps = rel_eps, max_it = 100L,
    c1 = 1e-4, c2 = .9, use_bfgs = FALSE,
    n_threads = 2L)
  expect_equal(opt[do_check], opt_res[do_check], tolerance = tol)

  opt <- optim_mlogit(
    val = val, ptr = optimizer, rel_eps = rel_eps, max_it = 100L,
    c1 = 1e-4, c2 = .9, n_threads = 1L, strong_wolfe = FALSE)
  opt_res <- list(par = c(0.647310437469711, 0.773205671911178, 0.656983834118283,
                          0.119681528415682, 0.855849653168437, -0.484311298187747, 0.412484752266496,
                          0.766994317977443, 0.422874453513229, 0.605169744625452, -0.783701817869368,
                          0.096194899967712, -0.922934757392211, -0.483542761956257, 0.438193243670319,
                          -0.493309514536335, 0.0274193762771719, 0.426816737611666, 0.0281704108821127,
                          -0.180877909582702, -0.150512019701299), value = 25.1215015278071,
                  info = 0L, counts = c(`function` = 12, gradient = 11, n_cg = 40
                  ), convergence = TRUE)
  expect_equal(opt[do_check], opt_res[do_check], tolerance = tol)
})
