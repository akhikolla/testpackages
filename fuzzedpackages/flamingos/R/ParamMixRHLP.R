#' A Reference Class which contains parameters of a mixture of RHLP models.
#'
#' ParamMixRHLP contains all the parameters of a mixture of RHLP models.
#'
#' @field fData [FData][FData] object representing the sample (covariates/inputs
#'   `X` and observed responses/outputs `Y`).
#' @field K The number of clusters (Number of RHLP models).
#' @field R The number of regimes (RHLP components) for each cluster.
#' @field p The order of the polynomial regression.
#' @field q The dimension of the logistic regression. For the purpose of
#'   segmentation, it must be set to 1.
#' @field variance_type Character indicating if the model is homoskedastic
#'   (`variance_type = "homoskedastic"`) or heteroskedastic (`variance_type =
#'   "heteroskedastic"`). By default the model is heteroskedastic.
#' @field alpha Cluster weights. Matrix of dimension \eqn{(1, K)}.
#' @field W Parameters of the logistic process. \eqn{\boldsymbol{W} =
#'   (\boldsymbol{w}_{1},\dots,\boldsymbol{w}_{K})}{W = (w_{1},\dots,w_{K})} is
#'   an array of dimension \eqn{(q + 1, R - 1, K)}, with \eqn{\boldsymbol{w}_{k}
#'   = (\boldsymbol{w}_{k,1},\dots,\boldsymbol{w}_{k,R-1})}{w_{k} =
#'   (w_{k,1},\dots,w_{k,R-1})}, \eqn{k = 1,\dots,K}, and `q` the order of the
#'   logistic regression. `q` is fixed to 1 by default.
#' @field beta Parameters of the polynomial regressions. \eqn{\boldsymbol{\beta}
#'   = (\boldsymbol{\beta}_{1},\dots,\boldsymbol{\beta}_{K})}{\beta =
#'   (\beta_{1},\dots,\beta_{K})} is an array of dimension \eqn{(p + 1, R, K)},
#'   with \eqn{\boldsymbol{\beta}_{k} =
#'   (\boldsymbol{\beta}_{k,1},\dots,\boldsymbol{\beta}_{k,R})}{\beta_{k} =
#'   (\beta_{k,1},\dots,\beta_{k,R})}, \eqn{k = 1,\dots,K}, `p` the order of the
#'   polynomial regression. `p` is fixed to 3 by default.
#' @field sigma2 The variances for the `K` clusters. If MixRHLP model is
#'   heteroskedastic (`variance_type = "heteroskedastic"`) then `sigma2` is a
#'   matrix of size \eqn{(R, K)} (otherwise MixRHLP model is homoskedastic
#'   (`variance_type = "homoskedastic"`) and `sigma2` is a matrix of size
#'   \eqn{(K, 1)}).
#' @field nu The degree of freedom of the MixRHLP model representing the
#'   complexity of the model.
#' @field phi A list giving the regression design matrices for the polynomial
#'   and the logistic regressions.
#' @export
ParamMixRHLP <- setRefClass(
  "ParamMixRHLP",
  fields = list(
    fData = "FData",
    phi = "list",

    K = "numeric", # Number of clusters
    R = "numeric", # Number of regimes
    p = "numeric", # Dimension of beta (order of polynomial regression)
    q = "numeric", # Dimension of w (order of logistic regression)
    variance_type = "character",
    nu = "numeric", # Degree of freedom

    alpha = "matrix", # Cluster weights
    W = "array", # W = (w_1,...w_K), w_k = (W_{k1},...,w_{k,R-1}) parameters of the logistic process: matrix of dimension [(q+1)x(R-1)] with q the order of logistic regression.
    beta = "array", # beta = (beta_1,...,beta_K), beta_k = (beta_{k1},...,beta_{kR}) polynomial regression coefficients: matrix of dimension [(p+1)xR] p being the polynomial degree.
    sigma2 = "matrix" # sigma2 = (sigma^2_k1,...,sigma^2_KR) : the variances for the R regmies for each cluster $k$.
  ),
  methods = list(
    initialize = function(fData = FData(numeric(1), matrix(1)), K = 1, R = 1, p = 3, q = 1, variance_type = "heteroskedastic") {

      fData <<- fData

      phi <<- designmatrix(x = fData$X, p = p, q = q, n = fData$n)

      K <<- K
      R <<- R
      p <<- p
      q <<- q
      variance_type <<- variance_type

      if (variance_type == "homoskedastic") {
        nu <<- (K - 1) + K * ((q + 1) * (R - 1) + R * (p + 1) + 1)
      } else {
        nu <<- (K - 1) + K * ((q + 1) * (R - 1) + R * (p +  1) + R)
      }

      W <<- array(0, dim = c(q + 1, R - 1, K))
      beta <<- array(NA, dim = c(p + 1, R, K))
      if (variance_type == "homoskedastic") {
        sigma2 <<- matrix(NA, K)
      } else {
        sigma2 <<- matrix(NA, R, K)
      }
      alpha <<- matrix(NA, K)
    },

    initParam = function(init_kmeans = TRUE, try_algo = 1) {
      "Method to initialize parameters \\code{alpha}, \\code{W}, \\code{beta}
      and \\code{sigma2}.

      If \\code{init_kmeans = TRUE} then the curve partition is initialized by
      the R-means algorithm. Otherwise the curve partition is initialized
      randomly.

      If \\code{try_algo = 1} then \\code{beta} and \\code{sigma2} are
      initialized by segmenting  the time series \\code{Y} uniformly into
      \\code{R} contiguous segments. Otherwise, \\code{W}, \\code{beta} and
      \\code{sigma2} are initialized by segmenting randomly the time series
      \\code{Y} into \\code{R} segments."

      # 1. Initialization of cluster weights
      alpha <<- 1 / (K * ones(K, 1))

      # 2. Initialization of the model parameters for each cluster: W, betak and sigmak

      # Setting W
      if (try_algo == 1) {
        for (k in (1:K)) { # Random initialization of parameter vector for IRLS
          W[, , k] <<- zeros(q + 1, R - 1)
        }
      } else {
        for (k in (1:K)) { # Random initialization of parameter vector for IRLS
          W[, , k] <<- rand(q + 1, R - 1)
        }
      }

      # beta_kr and sigma2_kr
      if (init_kmeans) {

        kmeans_res <- kmeans(fData$Y, K, nbr_runs = 20, nbr_iter_max = 400, verbose = FALSE)
        klas <- kmeans_res$klas
        for (k in 1:K) {
          Yk <- fData$Y[klas == k,]
          initRegressionParam(Yk, k, try_algo)
        }
      } else {
        ind <- sample(fData$n)
        for (k in 1:K) {
          if (k < K) {
            Yk <- fData$Y[ind[((k - 1) * round(fData$n / K) + 1):(k * round(fData$n / K))],]
          } else {
            Yk <- fData$Y[ind[((k - 1) * round(fData$n / K) + 1):length(ind)],]
          }
          initRegressionParam(Yk, k, try_algo)
        }
      }
    },

    initRegressionParam = function(Yk, k, try_algo  = 1) {
      "Initialize the matrix of polynomial regression coefficients beta_k for
      the cluster \\code{k}."
      n <- nrow(Yk)
      m <- ncol(Yk)
      if (try_algo == 1) { # Uniform segmentation into R contiguous segments, and then a regression
        zi <- round(m / R) - 1

        beta_k <- matrix(NA, p + 1, R)
        sig2 <- c()

        for (r in 1:R) {
          i <- (r - 1) * zi + 1
          j <- r * zi
          Yij <- Yk[, i:j, drop = FALSE]
          Yij <- matrix(t(Yij), ncol = 1)
          phi_ij <- phi$XBeta[i:j, , drop = FALSE]
          Phi_ij <- repmat(phi_ij, n, 1)

          br <- solve(t(Phi_ij) %*% Phi_ij) %*% t(Phi_ij) %*% Yij
          beta_k[, r] <- br

          if (variance_type == "homoskedastic") {
            sig2 <- var(Yij)
          } else {
            mr <- j - i + 1 # length(Xij);
            z <- Yij - Phi_ij %*% br

            sr <- t(z) %*% z / (n * mr)

            sig2[r] <- sr

          }
        }
      } else {# Random segmentation into R contiguous segments, and then a regression
        Lmin <- round(m / R) # nbr pts min into one segment
        tr_init <- zeros(1, R + 1)
        R_1 <- R
        for (r in 2:R) {
          R_1 <- R_1 - 1

          temp <- tr_init[r - 1] + Lmin:(m - (R_1 * Lmin) - tr_init[r - 1])

          ind <- sample(length(temp))

          tr_init[r] <- temp[ind[1]]
        }
        tr_init[R + 1] <- m
        beta_r <- matrix(NA, p + 1, R)
        sig2 <- c()
        for (r in 1:R) {
          i <- tr_init[r] + 1
          j <- tr_init[r + 1]
          Yij <- Xg[, i:j]
          Yij <- matrix(t(Yij), ncol = 1)
          phi_ij <- phi$XBeta[i:j, , drop = FALSE]
          Phi_ij <- repmat(phi_ij, n, 1)

          br <- solve(t(Phi_ij) %*% Phi_ij) %*% t(Phi_ij) %*% Yij
          beta_r[, r] <- br

          if (variance_type == "homoskedastic") {
            sig2 <- var(Yij)
          } else {
            mr <- j - i + 1 #length(Xij);
            z <- Yij - Phi_ij %*% br
            sr <- t(z) %*% z / (n * mr)
            sig2[r] <- sr
          }
        }
      }

      beta[, , k] <<- beta_k
      if (variance_type == "homoskedastic") {
        sigma2[k] <<- sig2
      } else {
        sigma2[, k] <<- sig2
      }
    },

    CMStep = function(statMixRHLP, verbose_IRLS = FALSE) {
      "Method which implements the M-step of the CEM algorithm to learn the
      parameters of the MixRHLP model based on statistics provided by the
      object \\code{statMixRHLP} of class \\link{StatMixRHLP} (which contains
      the E-step and the C-step)."

      good_segmentation = TRUE

      alpha <<- t(colSums(statMixRHLP$z_ik)) / fData$n

      # Maximization w.r.t beta_kr et sigma2_kr
      cluster_labels <- t(repmat(statMixRHLP$klas, 1, fData$m)) # [m x n]
      cluster_labels <- as.vector(cluster_labels)

      for (k in 1:K) {
        Yk = fData$vecY[cluster_labels == k,] # Cluster k (found from a hard clustering)
        gamma_ijk <- as.matrix(statMixRHLP$gamma_ijkr[cluster_labels == k, , k]) #[(nk xm) x R]

        if (variance_type == "homoskedastic") {
          s <- 0
        } else {
          sigma2_kr <- zeros(R, 1)
        }

        beta_kr <- matrix(NA, p + 1, R)

        for (r in 1:R) {

          segments_weights <- gamma_ijk[, r, drop = F]
          phikr <- (sqrt(segments_weights) %*% ones(1, p + 1)) * phi$XBeta[cluster_labels == k,] # [(nk*m)*(p+1)]
          Ykr <- sqrt(segments_weights) * Yk

          # maximization w.r.t beta_kr: Weighted least squares
          beta_kr[, r] <- solve(t(phikr) %*% phikr + .Machine$double.eps * diag(p + 1)) %*% t(phikr) %*% Ykr # Maximization w.r.t beta_kr

          #    the same as
          #                 W_kr = diag(cluster_weights.*segment_weights);
          #                 beta_kr(:,r) = inv(phiBeta'*W_kr*phiBeta)*phiBeta'*W_kr*X;
          #   Maximization w.r.t au sigma2_kr :
          if (variance_type == "homoskedastic") {
            sr <- colSums((Ykr - phikr %*% beta_kr[, r]) ^ 2)
            s <- s + sr
            sigma2_kr <- s / sum(gamma_ijk)
          } else {
            sigma2_kr[r] <-
              colSums((Ykr - phikr %*% beta_kr[, r]) ^ 2) / (sum(segments_weights))
            if ((sum(segments_weights) == 0)) {
              good_segmentation = FALSE
              return(list(0, good_segmentation))
            }
          }
        }

        beta[, , k] <<- beta_kr
        if (variance_type == "homoskedastic") {
          sigma2[k] <<- sigma2_kr
        } else {
          sigma2[, k] <<- sigma2_kr

        }

        # Maximization w.r.t W

        # Setting of W[,,k]
        Wk_init <- matrix(W[, , k], nrow = q + 1)

        res_irls <- IRLS(phi$Xw[cluster_labels == k,], gamma_ijk, ones(nrow(gamma_ijk), 1), Wk_init, verbose_IRLS)

        W[, , k] <<- res_irls$W
        piik <- res_irls$piik
        reg_irls <- res_irls$reg_irls
      }
      return(list(reg_irls, good_segmentation))
    },

    MStep = function(statMixRHLP, verbose_IRLS = FALSE) {
      "Method which implements the M-step of the EM algorithm to learn the
      parameters of the MixRHLP model based on statistics provided by the
      object \\code{statMixRHLP} of class \\link{StatMixRHLP} (which contains
      the E-step)."

      alpha <<- t(colSums(statMixRHLP$tau_ik)) / fData$n
      for (k in 1:K) {
        temp <- repmat(statMixRHLP$tau_ik[, k], 1, fData$m) # [m x n]
        cluster_weights <-
          matrix(t(temp), fData$m * fData$n, 1) # Cluster_weights(:) [mn x 1]
        gamma_ijk <- as.matrix(statMixRHLP$gamma_ijkr[, , k]) # [(nxm) x R]

        if (variance_type == "homoskedastic") {
          s <- 0
        } else {
          sigma2_kr <- zeros(R, 1)
        }

        beta_kr <- matrix(NA, p + 1, R)

        for (r in 1:R) {

          segments_weights <- gamma_ijk[, r, drop = F]
          phikr <- (sqrt(cluster_weights * segments_weights) %*% ones(1, p + 1)) * phi$XBeta #[(n*m)*(p+1)]
          Ykr <- sqrt(cluster_weights * segments_weights) * fData$vecY

          # Maximization w.r.t beta_kr: Weighted least squares
          beta_kr[, r] <- solve(t(phikr) %*% phikr + .Machine$double.eps * diag(p + 1)) %*% t(phikr) %*% Ykr # Maximization w.r.t beta_kr

          #    the same as
          #                 W_kr = diag(cluster_weights.*segment_weights);
          #                 beta_kr(:,r) = inv(phiBeta'*W_kr*phiBeta)*phiBeta'*W_kr*X;
          #   Maximization w.r.t au sigma2_kr :
          if (variance_type == "homoskedastic") {
            sr <- colSums((Ykr - phikr %*% beta_kr[, r]) ^ 2)
            s <- s + sr
            sigma2_kr <- s / sum(colSums((cluster_weights %*% ones(1, R)) * gamma_ijk))
          } else {
            sigma2_kr[r] <- colSums((Ykr - phikr %*% beta_kr[, r]) ^ 2) / (colSums(cluster_weights * segments_weights))
          }
        }

        beta[, , k] <<- beta_kr
        if (variance_type == "homoskedastic") {
          sigma2[k] <<- sigma2_kr
        } else {
          sigma2[, k] <<- sigma2_kr
        }

        # Maximization w.r.t W
        # Setting of W[,,k]
        Wk_init <- matrix(W[, , k], nrow = q + 1)

        res_irls <- IRLS(phi$Xw, gamma_ijk, cluster_weights, Wk_init, verbose_IRLS)

        W[, , k] <<- res_irls$W
        piik <- res_irls$piik

      }
    }
  )
)
