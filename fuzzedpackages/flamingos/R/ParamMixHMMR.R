#' A Reference Class which contains parameters of a mixture of HMMR models.
#'
#' ParamMixHMMR contains all the parameters of a mixture of HMMR models.
#'
#' @field fData [FData][FData] object representing the sample (covariates/inputs
#'   `X` and observed responses/outputs `Y`).
#' @field K The number of clusters (Number of HMMR models).
#' @field R The number of regimes (HMMR components) for each cluster.
#' @field p The order of the polynomial regression.
#' @field variance_type Character indicating if the model is homoskedastic
#'   (`variance_type = "homoskedastic"`) or heteroskedastic (`variance_type =
#'   "heteroskedastic"`). By default the model is heteroskedastic.
#' @field order_constraint A logical indicating whether or not a mask of order
#'   one should be applied to the transition matrix of the Markov chain to
#'   provide ordered states. For the purpose of segmentation, it must be set to
#'   `TRUE` (which is the default value).
#' @field alpha Cluster weights. Matrix of dimension \eqn{(K, 1)}.
#' @field prior The prior probabilities of the Markov chains. `prior` is a
#'   matrix of dimension \eqn{(R, K)}. The k-th column represents the prior
#'   distribution of the Markov chain asociated to the cluster k.
#' @field trans_mat The transition matrices of the Markov chains. `trans_mat` is
#'   an array of dimension \eqn{(R, R, K)}.
#' @field mask Mask applied to the transition matrices `trans_mat`. By default,
#'   a mask of order one is applied.
#' @field beta Parameters of the polynomial regressions. `beta` is an array of
#'   dimension \eqn{(p + 1, R, K)}, with `p` the order of the polynomial
#'   regression. `p` is fixed to 3 by default.
#' @field sigma2 The variances for the `K` clusters. If MixHMMR model is
#'   heteroskedastic (`variance_type = "heteroskedastic"`) then `sigma2` is a
#'   matrix of size \eqn{(R, K)} (otherwise MixHMMR model is homoskedastic
#'   (`variance_type = "homoskedastic"`) and `sigma2` is a matrix of size
#' @field nu The degree of freedom of the MixHMMR model representing the
#'   complexity of the model.
#' @field phi A list giving the regression design matrix for the polynomial regressions.
#' @export
ParamMixHMMR <- setRefClass(
  "ParamMixHMMR",
  fields = list(
    fData = "FData",
    phi = "matrix",

    K = "numeric", # Number of clusters
    R = "numeric", # Number of regimes (HMM states)
    p = "numeric", # Dimension of beta (order of polynomial regression)
    variance_type = "character",
    order_constraint = "logical",
    nu = "numeric", # Degrees of freedom

    alpha = "matrix", # Cluster weights
    prior = "matrix", # Initial distributions
    trans_mat = "array", # Transition matrices
    beta = "array", # Polynomial regression coefficient vectors
    sigma2 = "matrix", # Variances
    mask = "matrix"
  ),
  methods = list(
    initialize = function(fData = FData(numeric(1), matrix(1)), K = 2, R = 1, p = 3, variance_type = "heteroskedastic", order_constraint = TRUE) {

      fData <<- fData

      phi <<- designmatrix(x = fData$X, p = p)$XBeta

      K <<- K
      R <<- R
      p <<- p

      variance_type <<- variance_type

      order_constraint <<- order_constraint

      if (order_constraint) {
        if (variance_type == "homoskedastic") {
          nu <<- (K - 1) + K * ((R - 1) + R + (R - 1) + R * (p + 1) + 1)
        } else {
          nu <<- (K - 1) + K * ((R - 1) + R + (R - 1) + R * (p + 1) + R)
        }
      } else {
        if (variance_type == "homoskedastic") {
          nu <<- (K - 1) + K * ((R - 1) + R * (R - 1) + R * (p + 1) + 1)
        } else {
          nu <<- (K - 1) + K * ((R - 1) + R * (R - 1) + R * (p + 1) + R)
        }
      }

      alpha <<- matrix(NA, nrow = K)
      prior <<- matrix(NA, nrow = R, ncol = K)
      trans_mat <<- array(NA, dim = c(R, R, K))
      beta <<- array(NA, dim = c(p + 1, R, K))

      if (variance_type == "homoskedastic") {
        sigma2 <<- matrix(NA, ncol = K)
      } else {
        sigma2 <<- matrix(NA, nrow = R, ncol = K)
      }
      mask <<- matrix(NA, R, R)

    },

    initParam = function(init_kmeans = TRUE, try_algo = 1) {
      "Method to initialize parameters \\code{alpha}, \\code{prior},
      \\code{trans_mat}, \\code{beta} and \\code{sigma2}.

      If \\code{init_kmeans = TRUE} then the curve partition is initialized by
      the K-means algorithm. Otherwise the curve partition is initialized
      randomly.

      If \\code{try_algo = 1} then \\code{beta} and \\code{sigma2} are
      initialized by segmenting  the time series \\code{Y} uniformly into
      \\code{R} contiguous segments. Otherwise, \\code{beta} and
      \\code{sigma2} are initialized by segmenting randomly the time series
      \\code{Y} into \\code{R} segments."

      # 1. Initialization of cluster weights
      alpha <<- 1 / K * matrix(1, K, 1)

      # Initialization of the initial distributions and the transition matrices
      if (order_constraint) {

        # Initialization taking into account the constraint:

        # Initialization of the transition matrix
        maskM <- diag(R) # Mask of order 1
        if (R > 1) {
          for (r in 1:(R - 1)) {
            ind <- which(maskM[r,] != 0)
            maskM[r, ind + 1] <- 1
          }
        }

        # Initialization of the initial distribution
        for (k in 1:K) {
          prior[, k] <<- c(1, matrix(0, R - 1, 1))
          trans_mat[, , k] <<- normalize(maskM, 2)$M
        }

        mask <<- maskM

      } else {

        for (k in 1:K) {
          prior[, k] <<- c(1, matrix(0, R - 1, 1))
          trans_mat[, , k] <<- mkStochastic(matrix(runif(R), R, R))
        }
      }

      # 2. Initialisation of regression coefficients and variances
      if (init_kmeans) {
        max_iter_kmeans <- 400
        n_tries_kmeans <- 20
        verbose_kmeans <- 0
        solution <- kmeans(fData$Y, K, n_tries_kmeans, max_iter_kmeans, verbose_kmeans)

        for (k in 1:K) {
          Yk <- fData$Y[solution$klas == k ,] #if kmeans
          initRegressionParam(Yk, k, R, phi, variance_type, try_algo)
        }

      } else {
        ind <- sample(1:fData$n, fData$n)
        for (k in 1:K) {
          if (k < K) {
            Yk <- fData$Y[ind[((k - 1) * round(fData$n / K) + 1):(k * round(fData$n / K))],]
          } else {
            Yk <- fData$Y[ind[((k - 1) * round(fData$n / K) + 1):fData$n],]
          }

          initRegressionParam(Yk, k, R, phi, variance_type, try_algo)

        }
      }
    },

    initRegressionParam = function(Y, k, R, phi, variance_type, try_algo) {
      "Initialize \\code{beta} and \\code{sigma2} for the cluster \\code{k}."

      n <- nrow(Y)
      m <- ncol(Y)

      if (variance_type == "homoskedastic") {
        s <- 0
      }

      if (try_algo == 1) {
        zi <- round(m / R) - 1
        for (r in 1:R) {
          i <- (r - 1) * zi + 1
          j <- r * zi

          Yij <- Y[, i:j, drop = FALSE]
          Yij <- matrix(t(Yij), ncol = 1, byrow = T)

          phi_ij <- phi[i:j,]
          Phi_ij <- repmat(phi_ij, n, 1)

          bk <- solve(t(Phi_ij) %*% Phi_ij) %*% t(Phi_ij) %*% Yij
          beta[, r, k] <<- bk

          if (variance_type == "homoskedastic") {
            s <- s + sum((Yij - Phi_ij %*% bk) ^ 2)
            sigma2[, k] <<- s / (n * m)
          } else {
            mk <- j - i + 1
            z <- Yij - Phi_ij %*% bk
            sk <- t(z) %*% z / (n * mk)
            sigma2[r, k] <<- sk
          }
        }

      } else {
        Lmin <- round(m / (R + 1))
        tr_init <- matrix(0, 1, R + 1)
        tr_init[1] <- 0
        R_1 <- R
        for (r in 2:R) {
          R_1 <- R_1 - 1
          temp <- seq(tr_init[r - 1] + Lmin, m - R_1 * Lmin)
          ind <- sample(length(temp))
          tr_init[r] <- temp[ind[1]]
        }

        tr_init[R + 1] <- m
        for (r in 1:R) {
          i <- tr_init[r] + 1
          j <- tr_init[r + 1]
          Yij <- Y[, i:j, drop = FALSE]
          Yij <- matrix(t(Yij), ncol = 1, byrow = T)

          phi_ij <- phi[i:j,]
          Phi_ij <- repmat(phi_ij, n, 1)

          bk <- solve(t(Phi_ij) %*% Phi_ij) %*% t(Phi_ij) %*% Yij
          beta[, r, k] <<- bk

          if (variance_type == "homoskedastic") {
            s <- s + sum((Yij - Phi_ij %*% bk) ^ 2)
            sigma2[, k] <<- s / (n * m)
          } else {
            mk <- j - i + 1
            z <- Yij - Phi_ij %*% bk
            sk <- t(z) %*% z / (n * mk)
            sigma2[r, k] <<- sk
          }
        }
      }
    },

    MStep = function(statMixHMMR) {
      "Method which implements the M-step of the EM algorithm to learn the
      parameters of the MixHMMR model based on statistics provided by the
      object \\code{statMixHMMR} of class \\link{StatMixHMMR} (which contains
      the E-step)."

      # Maximization of Q1 w.r.t alpha
      alpha <<- matrix(apply(statMixHMMR$tau_ik, 2, sum)) / fData$n

      exp_num_trans_k <- array(0, dim = c(R, R, fData$n))

      for (k in 1:K) {
        if (variance_type == "homoskedastic") {
          s <- 0
        }

        weights_cluster_k <- statMixHMMR$tau_ik[, k]

        # Maximization of Q2 w.r.t \pi_k
        exp_num_trans_k_from_l <- (matrix(1, R, 1) %*% t(weights_cluster_k)) * statMixHMMR$exp_num_trans_from_l[, , k] # [R x n]

        prior[, k] <<- (1 / sum(statMixHMMR$tau_ik[, k])) * apply(exp_num_trans_k_from_l, 1, sum) # sum over i

        # Maximization of Q3 w.r.t A_k
        for (r in 1:R) {
          if (fData$n == 1) {
            exp_num_trans_k[r, ,] <- t(matrix(1, R, 1) %*% weights_cluster_k) * drop(statMixHMMR$exp_num_trans[r, , , k])
          } else {
            exp_num_trans_k[r, ,] <- (matrix(1, R, 1) %*% t(weights_cluster_k)) * drop(statMixHMMR$exp_num_trans[r, , , k])
          }
        }

        if (fData$n == 1) {
          temp <- exp_num_trans_k
        } else {
          temp <- apply(exp_num_trans_k, MARGIN = c(1, 2), sum) # sum over i
        }

        trans_mat[, , k] <<- mkStochastic(temp)

        # If HMM with order constraints
        if (order_constraint) {
          trans_mat[, , k] <<- mkStochastic(mask * trans_mat[, , k])
        }

        # Maximisation of Q4 w.r.t with betak et sigma2k

        Nk <- apply(statMixHMMR$tau_ik, 2, sum) # Nbr of individuals within the cluster k ,k=1...K estimated at iteration q
        nk <- Nk # Cardinal nbr of the cluster k
        # Each sequence i (m observations) is first weighted by the cluster weights
        weights_cluster_k <- matrix(t(statMixHMMR$tau_ik[, k]), nrow = fData$m, ncol = ncol(t(statMixHMMR$tau_ik)), byrow = T)
        weights_cluster_k <- matrix(as.vector(weights_cluster_k), length(as.vector(weights_cluster_k)), 1)

        # Secondly, the m observations of each sequence are weighted by the
        # weights of each segment k (post prob of the segments for each
        # cluster k)

        gamma_ijk <- statMixHMMR$gamma_ikjr[, , k] # [n*m R]

        for (r in 1:R) {

          weights_seg_k <- matrix(as.matrix(gamma_ijk)[, r])

          Xkr <- (sqrt(weights_cluster_k * weights_seg_k) %*% matrix(1, 1, p + 1)) * repmat(phi, fData$n, 1) # [n*m x (p+1)]
          Ykr <- (sqrt(weights_cluster_k * weights_seg_k)) * fData$vecY # [n*m x 1]

          # Weighted least squares: maximization w.r.t beta
          beta[, r, k] <<- solve(t(Xkr) %*% Xkr) %*% t(Xkr) %*% Ykr # Maximization w.r.t beta

          # Maximization w.r.t sigma2k :
          z <- sqrt(weights_cluster_k * weights_seg_k) * (fData$vecY - repmat(phi, fData$n, 1) %*% beta[, r, k])

          if (variance_type == "homoskedastic") {
            s <- s + (t(z) %*% z)
            nkm <- sum((weights_cluster_k %*% matrix(1, 1, R)) * gamma_ijk)
            sigma2[k] <<- s / nkm
          } else {
            nkmr <- sum(weights_cluster_k * weights_seg_k)

            sigma2[r, k] <<- (t(z) %*% z) / nkmr
          }
        }
      }
    }
  )
)
