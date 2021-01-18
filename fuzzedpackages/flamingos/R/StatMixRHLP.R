#' A Reference Class which contains statistics of a mixture of RHLP models.
#'
#' StatMixRHLP contains all the statistics associated to a
#' [MixRHLP][ParamMixRHLP] model, in particular the E-Step (and C-Step) of the
#' (C)EM algorithm.
#'
#' @field pi_jkr Array of size \eqn{(nm, R, K)} representing the logistic
#'   proportion for cluster k.
#' @field tau_ik Matrix of size \eqn{(n, K)} giving the posterior probabilities
#'   (fuzzy segmentation matrix) that the curve \eqn{\boldsymbol{y}_{i}}{y_{i}}
#'   originates from the \eqn{k}-th RHLP model.
#' @field z_ik Hard segmentation logical matrix of dimension \eqn{(n, K)}
#'   obtained by the Maximum a posteriori (MAP) rule: \eqn{z\_ik = 1 \
#'   \textrm{if} \ z\_i = \textrm{arg} \ \textrm{max}_{k} \ tau\_ik;\ 0 \
#'   \textrm{otherwise}}{z_ik = 1 if z_i = arg max_k tau_ik; 0 otherwise}.
#' @field klas Column matrix of the labels issued from `z_ik`. Its elements are
#'   \eqn{klas[i] = z\_i}{klas[i] = z_i}, \eqn{i = 1,\dots,n}.
#' @field gamma_ijkr Array of size \eqn{(nm, R, K)} giving the posterior
#'   probabilities that the observation \eqn{\boldsymbol{y}_{ij}}{y_{ij}}
#'   originates from the \eqn{r}-th regime of the \eqn{k}-th RHLP model.
#' @field polynomials Array of size \eqn{(m, R, K)} giving the values of the
#'   estimated polynomial regression components.
#' @field weighted_polynomials Array of size \eqn{(m, R, K)} giving the values
#'   of the estimated polynomial regression components weighted by the prior
#'   probabilities `pi_jkr`.
#' @field Ey Matrix of size \emph{(m, K)}. `Ey` is the curve expectation
#'   (estimated signal): sum of the polynomial components weighted by the
#'   logistic probabilities `pi_jkr`.
#' @field loglik Numeric. Observed-data log-likelihood of the MixRHLP model.
#' @field com_loglik Numeric. Complete-data log-likelihood of the MixRHLP model.
#' @field stored_loglik Numeric vector. Stored values of the log-likelihood at
#'   each EM iteration.
#' @field stored_com_loglik Numeric vector. Stored values of the Complete
#'   log-likelihood at each EM iteration.
#' @field BIC Numeric. Value of BIC (Bayesian Information Criterion).
#' @field ICL Numeric. Value of ICL (Integrated Completed Likelihood).
#' @field AIC Numeric. Value of AIC (Akaike Information Criterion).
#' @field log_fk_yij Matrix of size \eqn{(n, K)} giving the values of the
#'   probability density function \eqn{f(\boldsymbol{y}_{i} | z_i = k,
#'   \boldsymbol{x}, \boldsymbol{\Psi})}{f(y_{i} | z_i = k, x, \Psi)}, \eqn{i =
#'   1,\dots,n}.
#' @field log_alphak_fk_yij Matrix of size \eqn{(n, K)} giving the values of the
#'   logarithm of the joint probability density function
#'   \eqn{f(\boldsymbol{y}_{i}, \ z_{i} = k | \boldsymbol{x},
#'   \boldsymbol{\Psi})}{f(y_{i}, z_{i} = k | x, \Psi)}, \eqn{i = 1,\dots,n}.
#' @field log_gamma_ijkr Array of size \eqn{(nm, R, K)} giving the logarithm of
#'   `gamma_ijkr`.
#' @seealso [ParamMixRHLP]
#' @export
StatMixRHLP <- setRefClass(
  "StatMixRHLP",
  fields = list(
    pi_jkr = "array", # pi_jkr :logistic proportions for cluster k
    tau_ik = "matrix", # tau_ik = prob(curve|cluster_k) : post prob (fuzzy segmentation matrix of dim [nxK])
    z_ik = "matrix", # z_ik : Hard partition obtained by the MAP rule:  z_{ik} = 1
    # if and only z_i = arg max_k tau_ik (k=1,...,K)
    klas = "matrix", # klas : column vector of cluster labels
    Ey = "matrix",
    # Ey: curve expectation: sum of the polynomial components beta_kr weighted by
    # the logitic probabilities pi_jkr: Ey(j) = sum_{r=1}^R pi_jkr beta_kr rj, j=1,...,m. Ey
    # is a column vector of dimension m for each k.
    loglik = "numeric", # the loglikelihood of the EM or CEM algorithm
    com_loglik = "numeric", # the complete loglikelihood of the EM (computed at the convergence) or CEM algorithm
    stored_loglik = "numeric", # vector of stored valued of the comp-log-lik at each EM teration
    stored_com_loglik = "numeric",
    gamma_ijkr = "array",
    # gamma_ijkr prob(y_{ij}|kth_segment,cluster_g), fuzzy
    # segmentation for the cluster k. matrix of dimension
    # [nmxR] for each k  (k=1,...,K).
    log_gamma_ijkr = "array",
    BIC = "numeric", # BIC value = loglik - nu*log(nm)/2.
    ICL = "numeric", # ICL value = comp-loglik_star - nu*log(nm)/2.
    AIC = "numeric", # AIC value = loglik - nu.
    log_fk_yij = "matrix",
    log_alphak_fk_yij = "matrix",
    polynomials = "array",
    weighted_polynomials = "array"
  ),
  methods = list(
    initialize = function(paramMixRHLP = ParamMixRHLP()) {

      pi_jkr <<- array(0, dim = c(paramMixRHLP$fData$m * paramMixRHLP$fData$n, paramMixRHLP$R, paramMixRHLP$K))
      tau_ik <<- matrix(NA, paramMixRHLP$fData$n, paramMixRHLP$K)
      z_ik <<- matrix(NA, paramMixRHLP$fData$n, paramMixRHLP$K)
      klas <<- matrix(NA, paramMixRHLP$fData$n, 1)
      Ey <<- matrix(NA, paramMixRHLP$fData$m, paramMixRHLP$K)
      loglik <<- -Inf
      com_loglik <<- -Inf
      stored_loglik <<- numeric()
      stored_com_loglik <<- numeric()
      BIC <<- -Inf
      ICL <<- -Inf
      AIC <<- -Inf
      log_fk_yij <<- matrix(0, paramMixRHLP$fData$n, paramMixRHLP$K)
      log_alphak_fk_yij <<- matrix(0, paramMixRHLP$fData$n, paramMixRHLP$K)
      polynomials <<- array(NA, dim = c(paramMixRHLP$fData$m, paramMixRHLP$R, paramMixRHLP$K))
      weighted_polynomials <<- array(NA, dim = c(paramMixRHLP$fData$m, paramMixRHLP$R, paramMixRHLP$K))
      gamma_ijkr <<- array(0, dim = c(paramMixRHLP$fData$n * paramMixRHLP$fData$m, paramMixRHLP$R, paramMixRHLP$K))
      log_gamma_ijkr <<- array(0, dim = c(paramMixRHLP$fData$n * paramMixRHLP$fData$m, paramMixRHLP$R, paramMixRHLP$K))
    },

    MAP = function() {
      "MAP calculates values of the fields \\code{z_ik} and \\code{klas}
      by applying the Maximum A Posteriori Bayes allocation rule.

      \\eqn{z\\_ik = 1 \\ \\textrm{if} \\ z\\_i = \\textrm{arg} \\
      \\textrm{max}_{k} \\ tau\\_ik;\\ 0 \\ \\textrm{otherwise}
      }{z_ik = 1 if z_i = arg max_k tau_ik; 0 otherwise}."

      N <- nrow(tau_ik)
      K <- ncol(tau_ik)
      ikmax <- max.col(tau_ik)
      ikmax <- matrix(ikmax, ncol = 1)
      z_ik <<- ikmax %*% ones(1, K) == ones(N, 1) %*% (1:K)
      klas <<- ones(N, 1)
      for (k in 1:K) {
        klas[z_ik[, k] == 1] <<- k
      }
    },

    computeStats = function(paramMixRHLP) {
      "Method used in the EM algorithm to compute statistics based on
      parameters provided by the object \\code{paramMixRHLP} of class
      \\link{ParamMixRHLP}."

      for (k in 1:paramMixRHLP$K) {

        polynomials[, , k] <<- paramMixRHLP$phi$XBeta[1:paramMixRHLP$fData$m, ] %*% matrix(paramMixRHLP$beta[, , k], nrow = paramMixRHLP$p + 1)

        weighted_polynomials[, , k] <<- pi_jkr[1:paramMixRHLP$fData$m, , k] * polynomials[, , k]
        Ey[, k] <<- rowSums(as.matrix(weighted_polynomials[, , k]))

      }

      Ey <<- matrix(Ey, nrow = paramMixRHLP$fData$m)

      BIC <<- loglik - (paramMixRHLP$nu * log(paramMixRHLP$fData$n) / 2)
      AIC <<- loglik - paramMixRHLP$nu

      cig_log_alphak_fk_yij <- (z_ik) * (log_alphak_fk_yij)

      com_loglik <<- sum(rowSums(cig_log_alphak_fk_yij))

      ICL <<- com_loglik - paramMixRHLP$nu * log(paramMixRHLP$fData$n) / 2
    },

    CStep = function(reg_irls) {
      "Method used in the CEM algorithm to update statistics."

      tau_ik <<- exp(lognormalize(log_alphak_fk_yij))

      MAP() # Setting klas and z_ik

      # Compute the optimized criterion
      cig_log_alphak_fk_yij <- (z_ik) * log_alphak_fk_yij
      com_loglik <<- sum(cig_log_alphak_fk_yij) +  reg_irls
    },

    EStep = function(paramMixRHLP) {
      "Method used in the EM algorithm to update statistics based on parameters
      provided by the object \\code{paramMixRHLP} of class \\link{ParamMixRHLP}
      (prior and posterior probabilities)."

      for (k in 1:paramMixRHLP$K) {

        alpha_k <- paramMixRHLP$alpha[k]
        beta_k <- matrix(paramMixRHLP$beta[, , k], nrow = paramMixRHLP$p + 1)
        Wk <- matrix(paramMixRHLP$W[, , k], nrow = paramMixRHLP$q + 1)
        piik <- multinomialLogit(Wk, paramMixRHLP$phi$Xw, ones(nrow(paramMixRHLP$phi$Xw), ncol(Wk) + 1), ones(nrow(paramMixRHLP$phi$Xw), 1))$piik
        pi_jkr[, , k] <<- as.matrix(repmat(piik[1:paramMixRHLP$fData$m,], paramMixRHLP$fData$n, 1))

        log_pijkr_fkr_yij <- zeros(paramMixRHLP$fData$n * paramMixRHLP$fData$m, paramMixRHLP$R)

        for (r in 1:paramMixRHLP$R) {

          beta_kr <- as.matrix(beta_k[, r])
          if (paramMixRHLP$variance_type == "homoskedastic") {
            s2kr <- paramMixRHLP$sigma2[k]
          } else {
            s2kr <- paramMixRHLP$sigma2[r, k]
          }
          z <- ((paramMixRHLP$fData$vecY - paramMixRHLP$phi$XBeta %*% beta_kr) ^ 2) / s2kr
          log_pijkr_fkr_yij[, r] <- log(pi_jkr[, r, k]) - 0.5 * (log(2 * pi) + log(s2kr)) - 0.5 * z # cond pdf of yij given z_i = k and h_i = r
        }

        log_pijkr_fkr_yij <- pmin(log_pijkr_fkr_yij, log(.Machine$double.xmax))
        log_pijkr_fkr_yij <- pmax(log_pijkr_fkr_yij, log(.Machine$double.xmin))

        pijkr_fkr_yij <- exp(log_pijkr_fkr_yij)
        sumk_pijkr_fkr_yij <- rowSums(pijkr_fkr_yij) # sum over k
        log_sumk_pijkr_fkr_yij <- log(sumk_pijkr_fkr_yij) # [n*m, 1]

        log_gamma_ijkr[, , k] <<- log_pijkr_fkr_yij - log_sumk_pijkr_fkr_yij %*% ones(1, paramMixRHLP$R)
        gamma_ijkr[, , k] <<- exp(lognormalize(log_gamma_ijkr[, , k]))

        log_fk_yij[, k] <<- rowSums(t(matrix(log_sumk_pijkr_fkr_yij, paramMixRHLP$fData$m, paramMixRHLP$fData$n))) # [n x 1]:  sum over j=1,...,m: fk_yij = prod_j sum_r pi_{jkr} N(x_{ij},mu_{kr},s_{kr))
        log_alphak_fk_yij[, k] <<- log(alpha_k) + log_fk_yij[, k] # [nxK]
      }
      log_alphak_fk_yij <<- pmin(log_alphak_fk_yij, log(.Machine$double.xmax))
      log_alphak_fk_yij <<- pmax(log_alphak_fk_yij, log(.Machine$double.xmin))

      tau_ik <<- exp(lognormalize(log_alphak_fk_yij))

      loglik <<- sum(log(rowSums(exp(log_alphak_fk_yij))))
    }

  )
)
