#' @include texmexFamily.R
#' @include gpd.info.R
#' @include gpd.sandwich.R
#' @export cgpd
NULL

cgpd <- texmexFamily(name = 'CGPD',
                     log.lik = function(data, th, ...) {
                       y <- data$y
                       X.phi <- data$D$phi
                       X.eta <- data$D$eta
                       n.phi <- ncol(X.phi)
                       n.end <- n.phi + ncol(X.eta)
                       function(param) {
                         stopifnot(length(param) == n.end)
                         phi <- X.phi %*% param[1:n.phi]
                         eta <- X.eta %*% param[(1 + n.phi):n.end]
                         sum(dgpd(y, exp(phi), exp(eta) - 1/2, u=th, log.d=TRUE))
                       }
                     }, # Close log.lik
                     param = c(phi=0, eta=exp(.5)),
                     info = NULL,
                     sandwich = NULL,
                     start = function(data){
                       y <- data$y
                       X.phi <- data$D[[1]]
                       X.eta <- data$D[[2]]
                       c(log(mean(y)), rep(1e-05, -1 + ncol(X.phi) + ncol(X.eta)))
                     }, # Close start

                     resid = function(o){
                       p <- texmexMakeParams(coef(o), o$data$D)
                       delta <- (o$data$y - o$threshold) / exp(p[,1])
                       .log1prel(delta * (exp(p[,2]) - 1/2)) * delta # Standard exponential
                     }, # Close resid

                     endpoint = function(param, model){
                       res <- model$threshold - exp(param[, 1]) / (exp(param[, 2]) - 0.5)
                       res[param[, 2] >= 0] <- Inf
                       res
                     },
                     rl = function(m, param, model){
                       ## write in terms of qgpd; let's not reinvent the wheel
                       qgpd(1/(m * model$rate),
                            exp(param[,1]), exp(param[,2]) - 0.5, u=model$threshold,
                            lower.tail=FALSE)
                     },
                     delta = function(param, m, model){
                       # This is not exact if a prior (penalty) function is used, but
                       # the CI is approximate anyway.
                       param <- c(model$rate, param)
                       out <- matrix(0, nrow=2, ncol=length(m))
                       if (param[3] == 0){ # exponential case
                         out[1,] <- exp(param[2]) * log(m * param[1])
                       } else {
                         # Next line contains exp(param[2]) because the derivative is of log(sigma), unlike in Coles page 82
                         out[1,] <- exp(param[2]) / (exp(param[3]) -1/2) * (m * param[1])^((exp(param[3]) - 1/2) - 1)

                         v <- exp(param[2]) / (exp(param[3]) - 1/2)
                         u <- ( (m * param[1])^(exp(param[3]) - 1/2) - 1)

                         dv <- -exp(param[2] + param[3]) / (exp(param[3]) - 1/2)^2
                         du <- (m * param[1])^(exp(param[3]) - 1/2) * log(m * param[1]) * exp(param[3])

                         out[2,] <- v * du + u * dv
                       }
                       out
                     }, # Close delta
                     density = function(n, param, model, log.d = FALSE){
                       dgpd(n, exp(c(param[, 1])), c(exp(param[, 2]) - 1/2), u=model$threshold, log.d = log.d)
                     },

                     rng = function(n, param, model){
                       rgpd(n, exp(c(param[, 1])), c(exp(param[, 2]) - 1/2), u=model$threshold)
                     },
                     prob = function(n, param, model){
                       pgpd(n, exp(c(param[, 1])), c(exp(param[, 2]) - 1/2), u=model$threshold)
                     },
                     quant = function(n, param, model){
                       qgpd(n, exp(c(param[, 1])), c(exp(param[, 2]) - 1/2), u=model$threshold)
                     }
)


