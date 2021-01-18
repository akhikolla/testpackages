#' A Reference Class which contains statistics of a mixture of HMM model.
#'
#' StatMixHMM contains all the statistics associated to a [MixHMM][ParamMixHMM]
#' model, in particular the E-Step of the EM algorithm.
#'
#' @field tau_ik Matrix of size \eqn{(n, K)} giving the posterior probabilities
#'   that the curve \eqn{\boldsymbol{y}_{i}}{y_{i}} originates from the
#'   \eqn{k}-th HMM model.
#' @field gamma_ikjr Array of size \eqn{(nm, R, K)} giving the posterior
#'   probabilities that the observation \eqn{\boldsymbol{y}_{ij}}{y_{ij}}
#'   originates from the \eqn{r}-th regime of the \eqn{k}-th HMM model.
#' @field loglik Numeric. Log-likelihood of the MixHMM model.
#' @field stored_loglik Numeric vector. Stored values of the log-likelihood at
#'   each iteration of the EM algorithm.
#' @field klas Row matrix of the labels issued from `tau_ik`. Its elements are
#'   \eqn{klas[i] = z\_i}{klas[i] = z_i}, \eqn{i = 1,\dots,n}.
#' @field z_ik Hard segmentation logical matrix of dimension \eqn{(n, K)}
#'   obtained by the Maximum a posteriori (MAP) rule: \eqn{z\_ik = 1 \
#'   \textrm{if} \ z\_i = \textrm{arg} \ \textrm{max}_{k} \ P(z_{ik} = 1 |
#'   \boldsymbol{y}_{i}; \boldsymbol{\Psi}) = tau\_tk;\ 0 \
#'   \textrm{otherwise}}{z_ik = 1 if z_i = arg max_k P(z_{ik} = 1 | y_{i};
#'   \Psi) = tau_ik; 0 otherwise}.
#' @field smoothed Matrix of size \eqn{(m, K)} giving the smoothed time series.
#'   The smoothed time series are computed by combining the time series
#'   \eqn{\boldsymbol{y}_{i}}{y_{i}} with both the estimated posterior regime
#'   probabilities `gamma_ikjr` and the corresponding estimated posterior
#'   cluster probability `tau_ik`. The k-th column gives the estimated mean
#'   series of cluster k.
#' @field BIC Numeric. Value of BIC (Bayesian Information Criterion).
#' @field AIC Numeric. Value of AIC (Akaike Information Criterion).
#' @field ICL1 Numeric. Value of ICL (Integrated Completed Likelihood
#'   Criterion).
#' @field log_alpha_k_fyi Private. Only defined for calculations.
#' @field exp_num_trans Private. Only defined for calculations.
#' @field exp_num_trans_from_l Private. Only defined for calculations.
#' @seealso [ParamMixHMM]
#' @export
StatMixHMM <- setRefClass(
  "StatMixHMM",
  fields = list(
    tau_ik = "matrix",
    gamma_ikjr = "array",
    log_alpha_k_fyi = "matrix",
    exp_num_trans = "array",
    exp_num_trans_from_l = "array",
    loglik = "numeric",
    stored_loglik = "numeric",
    klas = "matrix",
    z_ik = "matrix",
    smoothed = "matrix",
    mean_curves = "array",
    BIC = "numeric",
    AIC = "numeric",
    ICL1 = "numeric"
  ),
  methods = list(
    initialize = function(paramMixHMM = ParamMixHMM()) {
      tau_ik <<- matrix(NA, paramMixHMM$fData$n, paramMixHMM$K)
      gamma_ikjr <<- array(NA, dim = c(paramMixHMM$fData$n * paramMixHMM$fData$m, paramMixHMM$R, paramMixHMM$K))
      log_alpha_k_fyi <<- matrix(NA, paramMixHMM$fData$n, paramMixHMM$K)
      exp_num_trans <<- array(NA, dim = c(paramMixHMM$R, paramMixHMM$R, paramMixHMM$fData$n, paramMixHMM$K))
      exp_num_trans_from_l <<- array(NA, dim = c(paramMixHMM$R, paramMixHMM$fData$n, paramMixHMM$K))
      loglik <<- -Inf
      stored_loglik <<- numeric()
      klas <<- matrix(NA, paramMixHMM$fData$n, 1) # klas: [nx1 double]
      z_ik <<- matrix(NA, paramMixHMM$fData$n, paramMixHMM$K) # z_ik: [nxK]
      smoothed <<- matrix(NA, paramMixHMM$fData$m, paramMixHMM$K)
      mean_curves <<- array(NA, dim = c(paramMixHMM$fData$m, paramMixHMM$R, paramMixHMM$K))
      BIC <<- -Inf
      AIC <<- -Inf
      ICL1 <<- -Inf

    },

    MAP = function() {
      "MAP calculates values of the fields \\code{z_ik} and \\code{klas}
      by applying the Maximum A Posteriori Bayes allocation rule.

      \\eqn{z\\_ik = 1 \\ \\textrm{if} \\ z\\_i = \\textrm{arg} \\
      \\textrm{max}_{k} \\ P(z_{ik} = 1 | \\boldsymbol{y}_{i};
      \\boldsymbol{\\Psi}) = tau\\_tk;\\ 0 \\ \\textrm{otherwise}}{z_ik = 1 if
      z_i = arg max_k P(z_{ik} = 1 | y_{i};\\Psi) = tau_ik; 0 otherwise}."

      N <- nrow(tau_ik)
      K <- ncol(tau_ik)
      ikmax <- max.col(tau_ik)
      ikmax <- matrix(ikmax, ncol = 1)
      z_ik <<- ikmax %*% ones(1, K) == ones(N, 1) %*% (1:K) # partition_MAP
      klas <<- ones(N, 1)
      for (k in 1:K) {
        klas[z_ik[, k] == 1] <<- k
      }
    },

    computeStats = function(paramMixHMM) {
      "Method used in the EM algorithm to compute statistics based on
      parameters provided by the object \\code{paramMixHMM} of class
      \\link{ParamMixHMM}."

      for (k in 1:paramMixHMM$K) {
        weighted_segments <- apply(gamma_ikjr[, , k] * (paramMixHMM$fData$vecY %*% matrix(1, 1, paramMixHMM$R)), 1, sum)

        dim(weighted_segments) <- c(paramMixHMM$fData$m, paramMixHMM$fData$n)
        weighted_clusters <- (matrix(1, paramMixHMM$fData$m, 1) %*% t(tau_ik[, k])) * weighted_segments
        smoothed[, k] <<- (apply(weighted_clusters, 1, sum)) / sum(tau_ik[, k])
      }

      # BIC AIC and ICL*
      BIC <<- loglik - (paramMixHMM$nu * log(paramMixHMM$fData$n) / 2)
      AIC <<- loglik - paramMixHMM$nu
      # ICL*
      # Compute the comp-log-lik
      cik_log_alpha_k_fyi <- (z_ik) * (log_alpha_k_fyi)
      comp_loglik <- sum(cik_log_alpha_k_fyi)
      ICL1 <<- comp_loglik - paramMixHMM$nu * log(paramMixHMM$fData$n) / 2 # n*m/2

    },

    EStep = function(paramMixHMM) {
      "Method used in the EM algorithm to update statistics based on parameters
      provided by the object \\code{paramMixHMM} of class \\link{ParamMixHMM}
      (prior and posterior probabilities)."

      exp_num_trans_ck  <- array(0, dim = c(paramMixHMM$R, paramMixHMM$R, paramMixHMM$fData$n))
      exp_num_trans_from_l_ck <- matrix(0, paramMixHMM$R, paramMixHMM$fData$n)

      alpha_k_fyi <- matrix(0, paramMixHMM$fData$n, paramMixHMM$K)

      for (k in 1:paramMixHMM$K) {
        # Run a hmm for each sequence
        log_fkr_yij <- matrix(0, paramMixHMM$R, paramMixHMM$fData$m)
        fkr_yij <- matrix(0, paramMixHMM$R, paramMixHMM$fData$m)

        Li <- matrix(0, paramMixHMM$fData$n, 1) # To store the loglik for each example

        mu_kr <- paramMixHMM$mu[, k]
        num_log_post_prob <- matrix(0, paramMixHMM$fData$n, paramMixHMM$K)

        for (i in 1:paramMixHMM$fData$n) {
          y_i <- paramMixHMM$fData$Y[i,]

          for (r in 1:paramMixHMM$R) {
            mukr <- mu_kr[r]

            if (paramMixHMM$variance_type == "homoskedastic") {
              sigma2 <- paramMixHMM$sigma2[k]
              sig2kr <- sigma2
            } else {
              sigma2 <- paramMixHMM$sigma2[, k]
              sig2kr <- sigma2[r]
            }
            z <- ((y_i - mukr * matrix(1, 1, paramMixHMM$fData$m)) ^ 2) / sig2kr
            log_fkr_yij[r,] <- -0.5 * matrix(1, 1, paramMixHMM$fData$m) * (log(2 * pi) + log(sig2kr)) - 0.5 * z# log pdf yij |z_i = k, h_i = r
            fkr_yij[r,] <- dnorm(y_i, mukr * matrix(1, 1, paramMixHMM$fData$m), sqrt(sig2kr))

          }

          # Calcul of p(y) : forwards backwards
          fb <- forwardsBackwards(as.matrix(paramMixHMM$prior[, k]), as.matrix(paramMixHMM$trans_mat[, , k]), fkr_yij)

          gamma_ik <- fb$tau_tk
          xi_ik <- fb$xi_tk
          fwd_ik <- fb$alpha_tk
          backw_ik <- fb$beta_tk
          loglik_i <- fb$loglik

          Li[i] <- loglik_i # loglik for the ith curve  ( logProb(Yi)

          gamma_ikjr[(((i - 1) * paramMixHMM$fData$m + 1):(i * paramMixHMM$fData$m)), , k] <<- t(gamma_ik) # [n*m R K]

          exp_num_trans_ck[, , i] <- apply(xi_ik, MARGIN = c(1, 2), sum) # [R R n]
          exp_num_trans_from_l_ck[, i] <- gamma_ik[, 1] # [R x n]
        }

        exp_num_trans_from_l[, , k] <<- exp_num_trans_from_l_ck # [R n K]
        exp_num_trans[, , , k] <<- exp_num_trans_ck # [R R n K]

        # For the MAP partition:  the numerator of the cluster post probabilities
        # num_log_post_prob[,k] <- log(param$alpha[k]) + Li

        # For computing the global loglik
        alpha_k_fyi[, k] <- paramMixHMM$alpha[k] * exp(Li) # [nx1]

        log_alpha_k_fyi[, k] <<- log(paramMixHMM$alpha[k]) + Li
      }

      log_alpha_k_fyi <<- pmin(log_alpha_k_fyi, log(.Machine$double.xmax))
      log_alpha_k_fyi <<- pmax(log_alpha_k_fyi, log(.Machine$double.xmin))

      tau_ik <<- exp(log_alpha_k_fyi) / (apply(exp(log_alpha_k_fyi), 1, sum) %*% matrix(1, 1, paramMixHMM$K))

      # Log-likelihood
      loglik <<- sum(log(apply(exp(log_alpha_k_fyi), 1, sum)))

    }
  )
)
