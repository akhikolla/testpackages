## ----setup, include=FALSE-------------------------------------------------------------------------
op <- options(width = 100)
knitr::opts_chunk$set(
  collapse = TRUE,
  error = TRUE, # ensures compilation even if testthat checks fail
  comment = "#>"
)

## ----validity-pkg, eval=TRUE----------------------------------------------------------------------
library("fddm")
library("rtdists")
library("RWiener")
source(system.file("extdata", "Gondan_et_al_density.R", package = "fddm", mustWork = TRUE))

## ----validity-run, eval=TRUE----------------------------------------------------------------------
# Define parameter space
RT <- c(0.001, 0.1, 1, 10)
A <- c(0.5, 1, 5)
V <- c(-5, 0, 5)
t0 <- 1e-4 # must be nonzero for RWiener
W <- c(0.2, 0.5, 0.8)
SV <- c(0, 0.5, 1.5)
SV_THRESH <- 1e-6
eps <- 1e-6 # this is the setting from rtdists

nRT <- length(RT)
nA <- length(A)
nV <- length(V)
nW <- length(W)
nSV <- length(SV)
resp <- rep("lower", nRT) # for RWiener

fnames <- c("fs_SWSE_17", "fs_SWSE_14", "fb_SWSE_17", "fb_SWSE_14",
            "fs_Gon_17", "fs_Gon_14", "fb_Gon_17", "fb_Gon_14",
            "fs_Nav_17", "fs_Nav_14", "fb_Nav_17", "fb_Nav_14",
            "fl_Nav_09", "RWiener", "Gondan", "rtdists")
nf <- length(fnames)

res <- data.frame(matrix(ncol = 9, nrow = nf*nRT*nA*nV*nW*nSV))
colnames(res) <- c('rt', 'a', 'v', 'w', 'sv', 'FuncName', 'res', 'dif',
                   'log_res')
start <- 1
stop <- nf

# Loop through each combination of parameters and record results
for (rt in 1:nRT) {
  for (a in 1:nA) {
    for (v in 1:nV) {
      for (w in 1:nW) {
        for (sv in 1:nSV) {
          # add the rt, v, a, w, and function names to the dataframe
          res[start:stop, 1] <- rep(RT[rt], nf)
          res[start:stop, 2] <- rep(A[a]  , nf)
          res[start:stop, 3] <- rep(V[v]  , nf)
          res[start:stop, 4] <- rep(W[w]  , nf)
          res[start:stop, 5] <- rep(SV[sv], nf)
          res[start:stop, 6] <- fnames

          # calculate "lower" density
          res[start,    7] <- dfddm(rt = RT[rt], response = resp[rt], a = A[a],
                                    v = V[v], t0 = t0, w = W[w], sv = SV[sv],
                                    log = FALSE, n_terms_small = "SWSE",
                                    summation_small = "2017", scale = "small",
                                    err_tol = eps)
          res[start+1,  7] <- dfddm(rt = RT[rt], response = resp[rt], a = A[a],
                                    v = V[v], t0 = t0, w = W[w], sv = SV[sv],
                                    log = FALSE, n_terms_small = "SWSE",
                                    summation_small = "2014", scale = "small",
                                    err_tol = eps)
          res[start+2,  7] <- dfddm(rt = RT[rt], response = resp[rt], a = A[a],
                                    v = V[v], t0 = t0, w = W[w], sv = SV[sv],
                                    log = FALSE, n_terms_small = "SWSE",
                                    summation_small = "2017", scale = "both",
                                    err_tol = eps)
          res[start+3,  7] <- dfddm(rt = RT[rt], response = resp[rt], a = A[a],
                                    v = V[v], t0 = t0, w = W[w], sv = SV[sv],
                                    log = FALSE, n_terms_small = "SWSE",
                                    summation_small = "2014", scale = "both",
                                    err_tol = eps)
          res[start+4,  7] <- dfddm(rt = RT[rt], response = resp[rt], a = A[a],
                                    v = V[v], t0 = t0, w = W[w], sv = SV[sv],
                                    log = FALSE, n_terms_small = "Gondan",
                                    summation_small = "2017", scale = "small",
                                    err_tol = eps)
          res[start+5,  7] <- dfddm(rt = RT[rt], response = resp[rt], a = A[a],
                                    v = V[v], t0 = t0, w = W[w], sv = SV[sv],
                                    log = FALSE, n_terms_small = "Gondan",
                                    summation_small = "2014", scale = "small",
                                    err_tol = eps)
          res[start+6,  7] <- dfddm(rt = RT[rt], response = resp[rt], a = A[a],
                                    v = V[v], t0 = t0, w = W[w], sv = SV[sv],
                                    log = FALSE, n_terms_small = "Gondan",
                                    summation_small = "2017", scale = "both",
                                    err_tol = eps)
          res[start+7,  7] <- dfddm(rt = RT[rt], response = resp[rt], a = A[a],
                                    v = V[v], t0 = t0, w = W[w], sv = SV[sv],
                                    log = FALSE, n_terms_small = "Gondan",
                                    summation_small = "2014", scale = "both",
                                    err_tol = eps)
          res[start+8,  7] <- dfddm(rt = RT[rt], response = resp[rt], a = A[a],
                                    v = V[v], t0 = t0, w = W[w], sv = SV[sv],
                                    log = FALSE, n_terms_small = "Navarro",
                                    summation_small = "2017", scale = "small",
                                    err_tol = eps)
          res[start+9,  7] <- dfddm(rt = RT[rt], response = resp[rt], a = A[a],
                                    v = V[v], t0 = t0, w = W[w], sv = SV[sv],
                                    log = FALSE, n_terms_small = "Navarro",
                                    summation_small = "2014", scale = "small",
                                    err_tol = eps)
          res[start+10, 7] <- dfddm(rt = RT[rt], response = resp[rt], a = A[a],
                                    v = V[v], t0 = t0, w = W[w], sv = SV[sv],
                                    log = FALSE, n_terms_small = "Navarro",
                                    summation_small = "2017", scale = "both",
                                    err_tol = eps)
          res[start+11, 7] <- dfddm(rt = RT[rt], response = resp[rt], a = A[a],
                                    v = V[v], t0 = t0, w = W[w], sv = SV[sv],
                                    log = FALSE, n_terms_small = "Navarro",
                                    summation_small = "2014", scale = "both",
                                    err_tol = eps)
          res[start+12,  7] <- dfddm(rt = RT[rt], response = resp[rt], a = A[a],
                                    v = V[v], t0 = t0, w = W[w], sv = SV[sv],
                                    log = FALSE, n_terms_small = "Navarro",
                                    scale = "large", err_tol = eps)
          res[start+13, 7] <- dwiener(RT[rt], resp = resp[rt], alpha = A[a],
                                      delta = V[v], tau = t0, beta = W[w],
                                      give_log = FALSE)
          res[start+14, 7] <- fs(t = RT[rt]-t0, a = A[a], v = V[v],
                                 w = W[w], eps = eps)
          res[start+15, 7] <- ddiffusion(RT[rt], resp[rt], a = A[a], v = V[v],
                                        t0 = t0, z = W[w]*A[a], sv = SV[sv])
          if (sv > SV_THRESH) { # multiply to get density with sv
            t <- RT[rt] - t0
            M <- exp(V[v] * A[a] * W[w] + V[v]*V[v] * t / 2 +
                     (SV[sv]*SV[sv] * A[a]*A[a] * W[w]*W[w] -
                       2 * V[v] * A[a] * W[w] - V[v]*V[v] * t) /
                     (2 + 2 * SV[sv]*SV[sv] * t)) / sqrt(1 + SV[sv]*SV[sv] * t)
            res[start+13, 7] <- M * res[start+11, 7] # RWiener
            res[start+14, 7] <- M * res[start+12, 7] # Gondan_R
          }

          # calculate differences
          ans <- res[start + 2, 7] # use fb_SWSE_17 as truth
          res[start,    8] <- abs(res[start,    7] - ans)
          res[start+1,  8] <- abs(res[start+1,  7] - ans)
          res[start+2,  8] <- abs(res[start+2,  7] - ans)
          res[start+3,  8] <- abs(res[start+3,  7] - ans)
          res[start+4,  8] <- abs(res[start+4,  7] - ans)
          res[start+5,  8] <- abs(res[start+1,  7] - ans)
          res[start+6,  8] <- abs(res[start+6,  7] - ans)
          res[start+7,  8] <- abs(res[start+7,  7] - ans)
          res[start+8,  8] <- abs(res[start+8,  7] - ans)
          res[start+9,  8] <- abs(res[start+9,  7] - ans)
          res[start+10, 8] <- abs(res[start+10, 7] - ans)
          res[start+11, 8] <- abs(res[start+11, 7] - ans)
          res[start+12, 8] <- abs(res[start+12, 7] - ans)
          res[start+13, 8] <- abs(res[start+13, 7] - ans)
          res[start+14, 8] <- abs(res[start+14, 7] - ans)
          res[start+15, 8] <- abs(res[start+15, 7] - ans)

          # calculate log of "lower" density
          res[start,    9] <- dfddm(rt = RT[rt], response = resp[rt], a = A[a],
                                    v = V[v], t0 = t0, w = W[w], sv = SV[sv],
                                    log = TRUE, n_terms_small = "SWSE",
                                    summation_small = "2017", scale = "small",
                                    err_tol = eps)
          res[start+1,  9] <- dfddm(rt = RT[rt], response = resp[rt], a = A[a],
                                    v = V[v], t0 = t0, w = W[w], sv = SV[sv],
                                    log = TRUE, n_terms_small = "SWSE",
                                    summation_small = "2014", scale = "small",
                                    err_tol = eps)
          res[start+2,  9] <- dfddm(rt = RT[rt], response = resp[rt], a = A[a],
                                    v = V[v], t0 = t0, w = W[w], sv = SV[sv],
                                    log = TRUE, n_terms_small = "SWSE",
                                    summation_small = "2017", scale = "both",
                                    err_tol = eps)
          res[start+3,  9] <- dfddm(rt = RT[rt], response = resp[rt], a = A[a],
                                    v = V[v], t0 = t0, w = W[w], sv = SV[sv],
                                    log = TRUE, n_terms_small = "SWSE",
                                    summation_small = "2014", scale = "both",
                                    err_tol = eps)
          res[start+4,  9] <- dfddm(rt = RT[rt], response = resp[rt], a = A[a],
                                    v = V[v], t0 = t0, w = W[w], sv = SV[sv],
                                    log = TRUE, n_terms_small = "Gondan",
                                    summation_small = "2017", scale = "small",
                                    err_tol = eps)
          res[start+5,  9] <- dfddm(rt = RT[rt], response = resp[rt], a = A[a],
                                    v = V[v], t0 = t0, w = W[w], sv = SV[sv],
                                    log = TRUE, n_terms_small = "Gondan",
                                    summation_small = "2014", scale = "small",
                                    err_tol = eps)
          res[start+6,  9] <- dfddm(rt = RT[rt], response = resp[rt], a = A[a],
                                    v = V[v], t0 = t0, w = W[w], sv = SV[sv],
                                    log = TRUE, n_terms_small = "Gondan",
                                    summation_small = "2017", scale = "both",
                                    err_tol = eps)
          res[start+7,  9] <- dfddm(rt = RT[rt], response = resp[rt], a = A[a],
                                    v = V[v], t0 = t0, w = W[w], sv = SV[sv],
                                    log = TRUE, n_terms_small = "Gondan",
                                    summation_small = "2014", scale = "both",
                                    err_tol = eps)
          res[start+8,  9] <- dfddm(rt = RT[rt], response = resp[rt], a = A[a],
                                    v = V[v], t0 = t0, w = W[w], sv = SV[sv],
                                    log = TRUE, n_terms_small = "Navarro",
                                    summation_small = "2017", scale = "small",
                                    err_tol = eps)
          res[start+9,  9] <- dfddm(rt = RT[rt], response = resp[rt], a = A[a],
                                    v = V[v], t0 = t0, w = W[w], sv = SV[sv],
                                    log = TRUE, n_terms_small = "Navarro",
                                    summation_small = "2014", scale = "small",
                                    err_tol = eps)
          res[start+10, 9] <- dfddm(rt = RT[rt], response = resp[rt], a = A[a],
                                    v = V[v], t0 = t0, w = W[w], sv = SV[sv],
                                    log = TRUE, n_terms_small = "Navarro",
                                    summation_small = "2017", scale = "both",
                                    err_tol = eps)
          res[start+11, 9] <- dfddm(rt = RT[rt], response = resp[rt], a = A[a],
                                    v = V[v], t0 = t0, w = W[w], sv = SV[sv],
                                    log = TRUE, n_terms_small = "Navarro",
                                    summation_small = "2014", scale = "both",
                                    err_tol = eps)
          res[start+12,  9] <- dfddm(rt = RT[rt], response = resp[rt], a = A[a],
                                    v = V[v], t0 = t0, w = W[w], sv = SV[sv],
                                    log = TRUE, n_terms_small = "Navarro",
                                    scale = "large", err_tol = eps)
          res[start+13, 9] <- dwiener(RT[rt], resp = resp[rt], alpha = A[a],
                                      delta = V[v], tau = t0, beta = W[w],
                                      give_log = TRUE)
          res[start+14, 9] <- log(fs(t = RT[rt]-t0, a = A[a], v = V[v],
                                     w = W[w], eps = eps))
          res[start+15, 9] <- log(ddiffusion(RT[rt], resp[rt], a = A[a],
                                             v = V[v], t0 = t0, z = W[w]*A[a],
                                             sv = SV[sv]))
          if (sv > SV_THRESH) { # add to get log of density with sv
            t <- RT[rt] - t0
            M <- V[v] * A[a] * W[w] + V[v]*V[v] * t / 2 +
                 (SV[sv]*SV[sv] * A[a]*A[a] * W[w]*W[w] -
                  2 * V[v] * A[a] * W[w] - V[v]*V[v] * t) /
                 (2 + 2 * SV[sv]*SV[sv] * t) - 0.5 * log(1 + SV[sv]*SV[sv] * t)
            res[start+13, 9] <- M + res[start+11, 9] # RWiener
            res[start+14, 9] <- M + res[start+12, 9] # Gondan_R
          }

          # iterate start and stop values
          start = start + nf
          stop = stop + nf
        }
      }
    }
  }
}

## ----validity-test, eval=TRUE---------------------------------------------------------------------
library("testthat")

# Subset results
SWSE_s <- res[res[["FuncName"]] %in% fnames[c(1, 2)], ]
SWSE_b <- res[res[["FuncName"]] %in% fnames[c(3, 4)], ]
Gondan_s <- res[res[["FuncName"]] %in% fnames[c(5, 6)], ]
Gondan_b <- res[res[["FuncName"]] %in% fnames[c(7, 8)], ]
Navarro_s <- res[res[["FuncName"]] %in% fnames[c(9, 10)], ]
Navarro_b <- res[res[["FuncName"]] %in% fnames[c(11, 12)], ]
Navarro_l <- res[res[["FuncName"]] %in% fnames[13], ]
RWiener <- res[res[["FuncName"]] %in% fnames[14], ]
Gondan_R <- res[res[["FuncName"]] %in% fnames[15], ]
rtdists <- res[res[["FuncName"]] %in% fnames[16], ]


# Ensure all densities are non-negative
test_that("Non-negativity of densities", {
  expect_true(all(SWSE_s[["res"]] >= 0))
  expect_true(all(SWSE_b[["res"]] >= 0))
  expect_true(all(Gondan_s[["res"]] >= 0))
  expect_true(all(Gondan_b[["res"]] >= 0))
  expect_true(all(Navarro_s[["res"]] >= 0))
  expect_true(all(Navarro_b[["res"]] >= 0))
  expect_true(all(Navarro_l[["res"]] >= 0))
  expect_true(all(RWiener[["res"]] >= 0))
  expect_true(all(Gondan_R[["res"]] >= 0))
  expect_true(all(rtdists[["res"]] >= 0))
})

# Test accuracy within 2*eps (allows for convergence from above and below)
test_that("Consistency among internal methods", {
  expect_true(all(SWSE_s[["dif"]] < 2*eps))
  expect_true(all(SWSE_b[["dif"]] < 2*eps))
  expect_true(all(Gondan_s[["dif"]] < 2*eps))
  expect_true(all(Gondan_b[["dif"]] < 2*eps))
  expect_true(all(Navarro_s[["dif"]] < 2*eps))
  expect_true(all(Navarro_b[["dif"]] < 2*eps))
  expect_true(all(Navarro_l[["dif"]] < 2*eps))
})

test_that("Accuracy relative to established packages", {
  expect_true(all(RWiener[RWiener[["sv"]] < SV_THRESH, "dif"] < 2*eps)) # see KE 1
  expect_true(all(Gondan_R[Gondan_R[["sv"]] < SV_THRESH, "dif"] < 2*eps)) # see KE 1
  expect_true(all(rtdists[["dif"]] < 2*eps))
})

# Test consistency in log vs non-log (see KE 2)
test_that("Log-Consistency among internal methods", {
  expect_equal(SWSE_s[SWSE_s[["res"]] > eps*eps, "log_res"],
               log(SWSE_s[SWSE_s[["res"]] > eps*eps, "res"]))
  expect_equal(SWSE_b[SWSE_b[["res"]] > eps*eps, "log_res"],
               log(SWSE_b[SWSE_b[["res"]] > eps*eps, "res"]))
  expect_equal(Gondan_s[Gondan_s[["res"]] > eps*eps, "log_res"],
               log(Gondan_s[Gondan_s[["res"]] > eps*eps, "res"]))
  expect_equal(Gondan_b[Gondan_b[["res"]] > eps*eps, "log_res"],
               log(Gondan_b[Gondan_b[["res"]] > eps*eps, "res"]))
  expect_equal(Navarro_s[Navarro_s[["res"]] > eps*eps, "log_res"],
               log(Navarro_s[Navarro_s[["res"]] > eps*eps, "res"]))
  expect_equal(Navarro_b[Navarro_b[["res"]] > eps*eps, "log_res"],
               log(Navarro_b[Navarro_b[["res"]] > eps*eps, "res"]))
  expect_equal(Navarro_l[Navarro_l[["res"]] > eps*eps, "log_res"],
               log(Navarro_l[Navarro_l[["res"]] > eps*eps, "res"]))
})

test_that("Log-Consistency of established packages", {
  expect_equal(RWiener[RWiener[["res"]] > eps*eps, "log_res"],
               log(RWiener[RWiener[["res"]] > eps*eps, "res"]))
  expect_equal(Gondan_R[Gondan_R[["res"]] > eps*eps, "log_res"],
               log(Gondan_R[Gondan_R[["res"]] > eps*eps, "res"]))
  expect_equal(rtdists[rtdists[["res"]] > eps*eps, "log_res"],
               log(rtdists[rtdists[["res"]] > eps*eps, "res"]))
})

## ----known-errors, eval=TRUE----------------------------------------------------------------------
rt <- 1.5
t <- rt - 1e-4
a <- 0.5
v <- 4.5
w <- 0.5
eps <- 1e-6
sv <- 0.9
sv0 <- exp(-v*a*w - v*v*t/2) / (a*a) # for constant drift rate
sv0_9 <- exp((-2*v*a*w - v*v*t + sv*sv*a*a*w*w)/(2 + 2*sv*sv*t)) /
         (a*a*sqrt(1+sv*sv*t)) # for variable drift rate
ks_0 <- ks(t/(a*a), w, eps/sv0) # = 2; the summation will only calculate 2 terms
ks_9 <- ks(t/(a*a), w, eps/sv0_9) # = 5; but the summation needs 5 terms

cat("the summation will only calculate", ks_0, "terms, but it needs", ks_9, "terms.")

## ----fitting-pkg, eval=TRUE-----------------------------------------------------------------------
library("fddm")
library("rtdists")

## ----loglik-fun, eval=TRUE------------------------------------------------------------------------
ll_fb_SWSE_17 <- function(pars, rt, resp, truth, err_tol) {
  v <- numeric(length(rt))
  v[truth == "upper"] <- pars[[1]]
  v[truth == "lower"] <- pars[[2]]
  dens <- dfddm(rt = rt, response = resp, a = pars[[3]], v = v,
                t0 = pars[[4]], w = pars[[5]], sv = pars[[6]], log = TRUE,
                n_terms_small = "SWSE", summation_small = "2017",
                scale = "both", err_tol = err_tol)
  return( ifelse(any(!is.finite(dens)), 1e6, -sum(dens)) )
}

ll_fb_Gon_17 <- function(pars, rt, resp, truth, err_tol) {
  v <- numeric(length(rt))
  v[truth == "upper"] <- pars[[1]]
  v[truth == "lower"] <- pars[[2]]
  dens <- dfddm(rt = rt, response = resp, a = pars[[3]], v = v,
                t0 = pars[[4]], w = pars[[5]], sv = pars[[6]], log = TRUE,
                n_terms_small = "Gondan", summation_small = "2017",
                scale = "both", err_tol = err_tol)
  return( ifelse(any(!is.finite(dens)), 1e6, -sum(dens)) )
}

ll_fb_Nav_17 <- function(pars, rt, resp, truth, err_tol) {
  v <- numeric(length(rt))
  v[truth == "upper"] <- pars[[1]]
  v[truth == "lower"] <- pars[[2]]
  dens <- dfddm(rt = rt, response = resp, a = pars[[3]], v = v,
                t0 = pars[[4]], w = pars[[5]], sv = pars[[6]], log = TRUE,
                n_terms_small = "Navarro", summation_small = "2017",
                scale = "both", err_tol = err_tol)
  return( ifelse(any(!is.finite(dens)), 1e6, -sum(dens)) )
}

ll_RTDists <- function(pars, rt, resp, truth) {
  rtu <- rt[truth == "upper"]
  rtl <- rt[truth == "lower"]
  respu <- resp[truth == "upper"]
  respl <- resp[truth == "lower"]

  densu <- ddiffusion(rtu, respu, a = pars[[3]], v = pars[[1]],
                      z = pars[[5]]*pars[[3]], t0 = pars[[4]], sv = pars[[6]])
  densl <- ddiffusion(rtl, respl, a = pars[[3]], v = pars[[2]],
                      z = pars[[5]]*pars[[3]], t0 = pars[[4]], sv = pars[[6]])

  densities <- c(densu, densl)
  if (any(densities <= 0)) return(1e6)
  return(-sum(log(densities)))
}

## ----fitting-fun, eval=TRUE-----------------------------------------------------------------------
rt_fit <- function(data, id_idx = NULL, rt_idx = NULL, response_idx = NULL,
                   truth_idx = NULL, response_upper = NULL, err_tol = 1e-6) {

  # Format data for fitting
  if (all(is.null(id_idx), is.null(rt_idx), is.null(response_idx),
      is.null(truth_idx), is.null(response_upper))) {
    df <- data # assume input data is already formatted
  } else {
    if(any(data[,rt_idx] < 0)) {
      stop("Input data contains negative response times; fit will not be run.")
    }
    if(any(is.na(data[,response_idx]))) {
      stop("Input data contains invalid responses (NA); fit will not be run.")
    }

    nr <- nrow(data)
    df <- data.frame(id = character(nr),
                     rt = double(nr),
                     response = character(nr),
                     truth = character(nr),
                     stringsAsFactors = FALSE)

    if (!is.null(id_idx)) { # relabel identification tags
      for (i in 1:length(id_idx)) {
        idi <- unique(data[,id_idx[i]])
        for (j in 1:length(idi)) {
          df[["id"]][data[,id_idx[i]] == idi[j]] <- paste(
            df[["id"]][data[,id_idx[i]] == idi[j]], idi[j], sep = " ")
        }
      }
      df[["id"]] <- trimws(df[["id"]], which = "left")
    }

    df[["rt"]] <- as.double(data[,rt_idx])

    df[["response"]] <- "lower"
    df[["response"]][data[,response_idx] == response_upper] <- "upper"

    df[["truth"]] <- "lower"
    df[["truth"]][data[,truth_idx] == response_upper] <- "upper"
  }

  # Preliminaries
  ids <- unique(df[["id"]])
  nids <- max(length(ids), 1) # if inds is null, there is only one individual

  init_vals <- data.frame(vu = c( 0,  10, -.5,  0,  0,  0,  0,  0,  0,   0,  0),
                          vl = c( 0, -10,  .5,  0,  0,  0,  0,  0,  0,   0,  0),
                          a  = c( 1,   1,   1, .5,  5,  1,  1,  1,  1,   1,  1),
                          t0 = c( 0,   0,   0,  0,  0,  0,  0,  0,  0,   0,  0),
                          w  = c(.5,  .5,  .5, .5, .5, .5, .5, .2, .8,  .5, .5),
                          sv = c( 1,   1,   1,  1,  1,  1,  1,  1,  1, .05,  5))
  ninit_vals <- nrow(init_vals)

  algo_names <- c("fb_SWSE_17", "fb_Gon_17", "fb_Nav_17", "rtdists")
  nalgos <- length(algo_names)
  ni <- nalgos*ninit_vals

  # Initilize the result dataframe
  cnames <- c("ID", "Algorithm", "Convergence", "Objective",
              "vu_init", "vl_init", "a_init", "t0_init", "w_init", "sv_init",
              "vu_fit", "vl_fit", "a_fit", "t0_fit", "w_fit", "sv_fit")
  res <- data.frame(matrix(ncol = length(cnames), nrow = nids*ninit_vals*nalgos))
  colnames(res) <- cnames

  # label the result dataframe
  res[["ID"]] <- rep(ids, each = ni) # label individuals
  res[["Algorithm"]] <- rep(algo_names, each = ninit_vals) # label algorithms
  res[["vu_init"]] <- init_vals[["vu"]] # label initial vu
  res[["vl_init"]] <- init_vals[["vl"]] # label initial vl
  res[["a_init"]]  <- init_vals[["a"]]  # label initial a
  res[["w_init"]]  <- init_vals[["w"]]  # label initial w
  res[["sv_init"]] <- init_vals[["sv"]] # label initial sv

  # Loop through each individual and starting values
  for (i in 1:nids) {
    # extract data for id i
    dfi <- df[df[["id"]] == ids[i], ]
    rti <- dfi[["rt"]]
    respi <- dfi[["response"]]
    truthi <- dfi[["truth"]]

    # starting value for t0 must be smaller than the smallest rt
    min_rti <- min(rti)
    t0_lo <- 0.01*min_rti
    t0_me <- 0.50*min_rti
    t0_hi <- 0.99*min_rti
    init_vals[["t0"]] <- c(rep(t0_me, 5), t0_lo, t0_hi, rep(t0_me, 4))

    # label the result dataframe
    res[["t0_init"]][((i-1)*ni+1):(i*ni)] <- init_vals[["t0"]] # label initial t0

    # loop through all of the starting values
    for (j in 1:ninit_vals) {
      temp <- nlminb(init_vals[j, ], ll_fb_SWSE_17,
                     rt = rti, resp = respi, truth = truthi, err_tol = err_tol,
                     # limits:   vu,   vl,   a,      t0, w,  sv
                     lower = c(-Inf, -Inf, .01,       0, 0,   0),
                     upper = c( Inf,  Inf, Inf, min_rti, 1, Inf))
      res[["Convergence"]][(i-1)*ni+0*ninit_vals+j] <- temp[["convergence"]]
      res[["Objective"]][(i-1)*ni+0*ninit_vals+j] <- temp[["objective"]]
      res[(i-1)*ni+0*ninit_vals+j, 11:16] <- temp[["par"]]

      temp <- nlminb(init_vals[j, ], ll_fb_Gon_17,
                     rt = rti, resp = respi, truth = truthi, err_tol = err_tol,
                     # limits:   vu,   vl,   a,      t0, w,  sv
                     lower = c(-Inf, -Inf, .01,       0, 0,   0),
                     upper = c( Inf,  Inf, Inf, min_rti, 1, Inf))
      res[["Convergence"]][(i-1)*ni+1*ninit_vals+j] <- temp[["convergence"]]
      res[["Objective"]][(i-1)*ni+1*ninit_vals+j] <- temp[["objective"]]
      res[(i-1)*ni+1*ninit_vals+j, 11:16] <- temp[["par"]]

      temp <- nlminb(init_vals[j, ], ll_fb_Nav_17,
                     rt = rti, resp = respi, truth = truthi, err_tol = err_tol,
                     # limits:   vu,   vl,   a,      t0, w,  sv
                     lower = c(-Inf, -Inf, .01,       0, 0,   0),
                     upper = c( Inf,  Inf, Inf, min_rti, 1, Inf))
      res[["Convergence"]][(i-1)*ni+2*ninit_vals+j] <- temp[["convergence"]]
      res[["Objective"]][(i-1)*ni+2*ninit_vals+j] <- temp[["objective"]]
      res[(i-1)*ni+2*ninit_vals+j, 11:16] <- temp[["par"]]

      temp <- nlminb(init_vals[j, ], ll_RTDists,
                     rt = rti, resp = respi, truth = truthi,
                     # limits:   vu,   vl,   a,      t0, w,  sv
                     lower = c(-Inf, -Inf, .01,       0, 0,   0),
                     upper = c( Inf,  Inf, Inf, min_rti, 1, Inf))
      res[["Convergence"]][(i-1)*ni+3*ninit_vals+j] <- temp[["convergence"]]
      res[["Objective"]][(i-1)*ni+3*ninit_vals+j] <- temp[["objective"]]
      res[(i-1)*ni+3*ninit_vals+j, 11:16] <- temp[["par"]]
    }
  }
  return(res)
}

## ----fitting-run, eval=FALSE----------------------------------------------------------------------
#  data(med_dec, package = "fddm")
#  med_dec <- med_dec[which(med_dec[["rt"]] >= 0), ]
#  fit <- rt_fit(med_dec, id_idx = c(2,1), rt_idx = 8, response_idx = 7,
#                truth_idx = 5, response_upper = "blast", err_tol = 1e-6)

## ----fitting-run-internal, eval=FALSE, include=FALSE----------------------------------------------
#  save(fit, compress = "xz", compression_level = 9,
#       file = "inst/extdata/valid_fit.Rds")

## ----fitting-prep, eval=TRUE----------------------------------------------------------------------
fit_prep <- function(fit) {
  nr <- nrow(fit)
  fit[["Obj_diff"]] <- rep(0, nr)
  fit[["vu_diff"]]  <- rep(0, nr)
  fit[["vl_diff"]]  <- rep(0, nr)
  fit[["a_diff"]]   <- rep(0, nr)
  fit[["t0_diff"]]  <- rep(0, nr)
  fit[["w_diff"]]   <- rep(0, nr)
  fit[["sv_diff"]]  <- rep(0, nr)

  ids <- unique(fit[["ID"]])
  nids <- length(ids)
  algos <- unique(fit[["Algorithm"]])
  nalgos <- length(algos)

  fit_idx <- c(4, 11:16)
  dif_idx <- 17:23
  ninit <- nrow(fit[fit[["ID"]] == ids[1] & fit[["Algorithm"]] == algos[1], ])
  for (i in 1:nids) {
    for (j in 1:ninit) {
      actual_idx <- seq((i-1)*ninit*nalgos+j, i*ninit*nalgos, by = ninit)
      min_obj_idx <- actual_idx[which.min(fit[actual_idx, 4])]
      best_fit <- fit[min_obj_idx, fit_idx]
      for (k in 0:(nalgos-1)) {
        fit[(i-1)*(ninit*nalgos) + k*ninit + j, dif_idx] <-
          fit[(i-1)*(ninit*nalgos) + k*ninit + j, fit_idx] - best_fit

      }
    }
  }
  return(fit)
}

## ----fit-load, eval=TRUE--------------------------------------------------------------------------
# load data, will be in the variable 'fit'
load(system.file("extdata", "valid_fit.Rds", package = "fddm", mustWork = TRUE))
fit <- fit_prep(fit)

cat("Results for ID = experienced 2")
fit[(0:3)*11+1, ]

## ----fit-estimates, eval=TRUE---------------------------------------------------------------------
# Define error tolerance
eps <- 1e-4

out <- fit[unique(which(abs(fit[, c(3, 17:23)]) > eps, arr.ind = TRUE)[, 1]), ]
out[, -c(1:2)] <- zapsmall(out[, -c(1:2)])
out

## ----session-info, collapse=TRUE------------------------------------------------------------------
sessionInfo()

## ----reset-options, include=FALSE---------------------------------------------
options(op)  # reset options

