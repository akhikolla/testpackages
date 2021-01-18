#' A Reference Class which contains parameters of a mixture of HMM models.
#'
#' ParamMixHMM contains all the parameters of a mixture of HMM models.
#'
#' @field fData [FData][FData] object representing the sample (covariates/inputs
#'   `X` and observed responses/outputs `Y`).
#' @field K The number of clusters (Number of HMM models).
#' @field R The number of regimes (HMM components) for each cluster.
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
#' @field mu Means. Matrix of dimension \eqn{(R, K)}. The k-th column gives
#'   represents the k-th cluster and gives the means for the `R` regimes.
#' @field sigma2 The variances for the `K` clusters. If MixHMM model is
#'   heteroskedastic (`variance_type = "heteroskedastic"`) then `sigma2` is a
#'   matrix of size \eqn{(R, K)} (otherwise MixHMM model is homoskedastic
#'   (`variance_type = "homoskedastic"`) and `sigma2` is a matrix of size
#'   \eqn{(1, K)}).
#' @field nu The degrees of freedom of the MixHMM model representing the
#'   complexity of the model.
#' @export
ParamMixHMM <- setRefClass(
  "ParamMixHMM",
  fields = list(
    fData = "FData",

    K = "numeric", # Number of clusters
    R = "numeric", # Number of regimes (HMM states)
    variance_type = "character",
    order_constraint = "logical",
    nu = "numeric", # Degree of freedom

    alpha = "matrix", # Cluster weights
    prior = "matrix", # Initial distributions
    trans_mat = "array", # Transition matrices
    mu = "matrix", # Means
    sigma2 = "matrix", # Variances
    mask = "matrix"
  ),
  methods = list(
    initialize = function(fData = FData(numeric(1), matrix(1)), K = 2, R = 1, p = 3, variance_type = "heteroskedastic", order_constraint = TRUE) {
      fData <<- fData

      K <<- K
      R <<- R

      variance_type <<- variance_type

      order_constraint <<- order_constraint

      if (order_constraint) {
        if (variance_type == "homoskedastic") {
          nu <<- (K - 1) + K * ((R - 1) + R + (R - 1) + R + 1)
        } else {
          nu <<- (K - 1) + K * ((R - 1) + R + (R - 1) + R + R)
        }
      } else {
        if (variance_type == "homoskedastic") {
          nu <<- (K - 1) + K * ((R - 1) + R * (R - 1) + R + 1)
        } else {
          nu <<- (K - 1) + K * ((R - 1) + R * (R - 1) + R + R)
        }
      }

      alpha <<- matrix(NA, nrow = K)
      prior <<- matrix(NA, nrow = R, ncol = K)
      trans_mat <<- array(NA, dim = c(R, R, K))
      mu <<- matrix(NA, nrow = R, ncol = K)

      if (variance_type == "homoskedastic") {
        sigma2 <<- matrix(NA, ncol = K)
      } else {
        sigma2 <<- matrix(NA, nrow = R, ncol = K)
      }
      mask <<- matrix(1, R, R)

    },

    initParam = function(init_kmeans = TRUE, try_algo = 1) {
      "Method to initialize parameters \\code{alpha}, \\code{prior}, \\code{trans_mat},
      \\code{mu} and \\code{sigma2}.

      If \\code{init_kmeans = TRUE} then the curve partition is initialized by
      the K-means algorithm. Otherwise the curve partition is initialized
      randomly.

      If \\code{try_algo = 1} then \\code{mu} and \\code{sigma2} are
      initialized by segmenting  the time series \\code{Y} uniformly into
      \\code{R} contiguous segments. Otherwise, \\code{mu} and
      \\code{sigma2} are initialized by segmenting randomly the time series
      \\code{Y} into \\code{R} segments."

      # 1. Initialization of cluster weights
      alpha <<- 1 / K * matrix(1, K, 1)

      # Initialization of the initial distributions and the transition matrices
      if (order_constraint) { # Initialization taking into account the constraint:

        maskM <- diag(R) # Mask of order 1
        if (R > 1) {
          for (r in 1:(R - 1)) {
            ind <- which(maskM[r,] != 0)
            maskM[r, ind + 1] <- 1
          }
        }

        for (k in 1:K) {
          prior[, k] <<- c(1, matrix(0, R - 1, 1))
          trans_mat[, , k] <<- normalize(maskM, 2)$M
        }

        mask <<- maskM

      } else {

        for (k in 1:K) {
          prior[, k] <<- 1 / R * matrix(1, R, 1)
          trans_mat[, , k] <<- mkStochastic(matrix(runif(R), R, R))
        }
      }

      # Initialization of the means and variances
      if (init_kmeans) {
        max_iter_kmeans <- 400
        n_tries_kmeans <- 20
        verbose_kmeans <- 0
        solution <- kmeans(fData$Y, K, n_tries_kmeans, max_iter_kmeans, verbose_kmeans)

        for (k in 1:K) {
          Yk <- fData$Y[solution$klas == k ,] # If kmeans
          initGaussParamHmm(Yk, k, R, variance_type, try_algo)
        }

      } else {
        ind <- sample(1:fData$n, fData$n)
        for (k in 1:K) {
          if (k < K) {
            Yk <- fData$Y[ind[((k - 1) * round(fData$n / K) + 1):(k * round(fData$n / K))],]
          } else {
            Yk <- fData$Y[ind[((k - 1) * round(fData$n / K) + 1):fData$n],]
          }

          initGaussParamHmm(Yk, k, R, variance_type, try_algo)

        }
      }
    },

    initGaussParamHmm = function(Y, k, R, variance_type, try_algo) {
      "Initialize the means \\code{mu} and \\code{sigma2} for the cluster
      \\code{k}."

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
          Yij <- matrix(t(Yij), 1, byrow = T)
          mu[r, k] <<- mean(Yij)

          if (variance_type == "homoskedastic") {
            s <- s + sum((Yij - mu[r, k]) ^ 2)
            sigma2[, k] <<- s / (n * m)
          } else {
            m_r <- j - i + 1
            sigma2[r, k] <<- sum((Yij - mu[r, k]) ^ 2) / (n * m_r)
          }
        }

      } else {
        Lmin <- 2
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

          mu[r, k] <<- mean(Yij)

          if (variance_type == "homoskedastic") {
            s <- s + sum((Yij - mu[r, k]) ^ 2)
            sigma2[, k] <<- s / (n * m)
          } else {
            m_r <- j - i + 1
            sigma2[r, k] <<- sum((Yij - mu[r, k]) ^ 2) / (n * m_r)
          }
        }
      }
    },

    MStep = function(statMixHMM) {
      "Method which implements the M-step of the EM algorithm to learn the
      parameters of the MixHMM model based on statistics provided by the object
      \\code{statMixHMM} of class \\link{StatMixHMM} (which contains the
      E-step)."

      # Maximization of Q1 w.r.t alpha
      alpha <<- matrix(apply(statMixHMM$tau_ik, 2, sum)) / fData$n

      exp_num_trans_k <- array(0, dim = c(R, R, fData$n))

      for (k in 1:K) {
        if (variance_type == "homoskedastic") {
          s <- 0
        }

        weights_cluster_k <- statMixHMM$tau_ik[, k]

        # Maximization of Q2 w.r.t \pi_k
        exp_num_trans_k_from_l <- (matrix(1, R, 1) %*% t(weights_cluster_k)) * statMixHMM$exp_num_trans_from_l[, , k] # [K x n]

        prior[, k] <<- (1 / sum(statMixHMM$tau_ik[, k])) * apply(exp_num_trans_k_from_l, 1, sum) # Sum over i

        # Maximization of Q3 w.r.t A_k
        for (r in 1:R) {
          if (fData$n == 1) {
            exp_num_trans_k[r, ,] <- t(matrix(1, R, 1) %*% weights_cluster_k) * drop(statMixHMM$exp_num_trans[r, , , k])
          } else {
            exp_num_trans_k[r, ,] <- (matrix(1, R, 1) %*% t(weights_cluster_k)) * drop(statMixHMM$exp_num_trans[r, , , k])
          }
        }

        if (fData$n == 1) {
          temp <- exp_num_trans_k
        } else{
          temp <- apply(exp_num_trans_k, MARGIN = c(1, 2), sum) # Sum over i
        }

        trans_mat[, , k] <<- mkStochastic(temp)

        # If HMM with order constraints
        if (order_constraint) {
          trans_mat[, , k] <<- mkStochastic(mask * trans_mat[, , k])
        }

        # Maximization of Q4 w.r.t muk and sigmak
        # each sequence i (m observations) is first weighted by the cluster weights
        weights_cluster_k <- matrix(t(statMixHMM$tau_ik[, k]), nrow = fData$m, ncol = ncol(t(statMixHMM$tau_ik)), byrow = T)
        weights_cluster_k <- matrix(as.vector(weights_cluster_k), length(as.vector(weights_cluster_k)), 1)

        # Secondly, the m observations of each sequence are weighted by the
        # weights of each segment k (post prob of the segments for each
        # cluster k)
        gamma_ijk <- statMixHMM$gamma_ikjr[, , k] # [n*m K]

        for (r in 1:R) {

          # Maximization w.r.t the Gaussian means muk
          if (R == 1) {
            weights_seg_k <- matrix(gamma_ijk)
          } else {
            weights_seg_k <- matrix(gamma_ijk[, r])
          }

          mu[r, k] <<- (1 / sum(weights_cluster_k * weights_seg_k)) %*% sum((weights_cluster_k * weights_seg_k) * fData$vecY)

          # Maximization w.r.t the Gaussian variances sigma2k
          z <- sqrt(weights_cluster_k * weights_seg_k) * (fData$vecY - matrix(1, fData$n * fData$m, 1) * mu[r, k])

          if (variance_type == "homoskedastic") {
            s <- s + (t(z) %*% z)
            nkr <- sum((weights_cluster_k %*% matrix(1, 1, R)) * gamma_ijk)
            sigma2[k] <<- s / nkr
          } else{
            nkmr <- sum(weights_cluster_k * weights_seg_k)
            sigma2[r, k] <<- (t(z) %*% z) / (nkmr)
          }
        }

      }
    }
  )
)
