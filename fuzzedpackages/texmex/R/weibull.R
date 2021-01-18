#' @include texmexFamily.R
#' @include weibull.info.R
#' @export weibull
NULL

weibull <- texmexFamily(name = 'Weibull',
          log.lik = function(data, th, ...) {
            y <- data$y
            X.phi <- data$D$phi
            X.gamma <- data$D$gamma

            n.phi <- ncol(X.phi)
            n.end <- n.phi + ncol(X.gamma)

            function(param) {
              stopifnot(length(param) == n.end)
              phi <- X.phi %*% param[1:n.phi]
              gamma <- X.gamma %*% param[(1 + n.phi):n.end]
              sum(dweibull(y-th, shape=gamma, scale=exp(phi), log=TRUE))
            }
          }, # Close log.lik
          param = c(phi=0, gamma=1),
          info = weibull.info,
          sandwich = NULL,
          start = function(data){
            y <- data$y
            X.phi <- data$D$phi
            X.gamma <- data$D$gamma
            c(log(mean(y) - min(y)), rep(1e-05, -1 + ncol(X.phi) + ncol(X.gamma)))
          }, # Close start

          resid = function(o){
            p <- texmexMakeParams(coef(o), o$data$D)
            ((o$data$y - o$threshold) / exp(p[,1]))^p[,2]  # Standard exponential
          }, # Close resid

          endpoint = function(param, model){
            Inf
          },
          rl = function(m, param, model){
            ## write in terms of qweibull; let's not reinvent the wheel
            qweibull(1/(m * model$rate),
                 scale=exp(param[,1]), shape = param[,2],
                 lower.tail=FALSE) + model$threshold
          },
          delta = function(param,m,model){ # follows argument in Coles eqn (4.15) for GPD
            out <- rep(1,2) # can have vector output as call only ever assumes m length 1
            lmr <- log(m * model$rate)

            out[1] <- exp(param[1]) * lmr^(1/param[2])
            out[2] <- -exp(param[1])/param[2]^2 * (lmr)^(1/param[2]) * log(lmr)
            out
          }, # Close delta
          density = function(n, param, model, log.d=FALSE){
            dweibull(n-model$threshold, scale = c(exp(param[, 1])), shape = c(param[, 2]), log=log.d)
          },

          rng = function(n, param, model){
            rweibull(n, scale = c(exp(param[, 1])), shape = c(param[, 2])) + model$threshold
          },
          prob = function(n, param, model){
            pweibull(n-model$threshold, scale = c(exp(param[, 1])), shape = c(param[, 2]))
          },
          quant = function(n, param, model){
            qweibull(n, scale = c(exp(param[, 1])), shape = c(param[, 2])) + model$threshold
          }
)
