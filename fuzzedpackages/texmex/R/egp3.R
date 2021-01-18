#' @include texmexFamily.R
#' @export egp3
NULL

egp3 <- texmexFamily(name="EGP3",
                     param = c(lambda=0, phi=0, xi=0), # Standard exponential
                     density = function(n, param, model, log.d=FALSE){
                                 degp3(n, exp(c(param[, 1])), exp(c(param[, 2])), c(param[, 3]), u=model$threshold, log.d=log.d)
                     },
                     rng = function(n, param, model){
                             regp3(n, exp(c(param[, 1])), exp(c(param[, 2])), c(param[, 3]), u=model$threshold)
                     },
                     prob = function(n, param, model){
                              pegp3(n, exp(c(param[, 1])), exp(c(param[, 2])), c(param[, 3]), u=model$threshold)
                     },
                     quant = function(n, param, model){
                               qegp3(n, exp(c(param[, 1])), exp(c(param[, 2])), c(param[, 3]), u=model$threshold)
                     },
                     log.lik = function(data, th, ...){
                                 y <- data$y

                                 X.lambda <- data$D$lambda
                                 X.phi <- data$D$phi
                                 X.xi <- data$D$xi

                                 n.lambda <- ncol(X.lambda)
                                 n.phi <- n.lambda + ncol(X.phi)
                                 n.end <- n.phi + ncol(X.xi)

                                 function(param) {
                                   stopifnot(length(param) == n.end)
                                   lambda <- X.lambda %*% param[1:n.lambda]
                                   phi <- X.phi %*% param[(1 + n.lambda):n.phi]
                                   xi <- X.xi %*% param[(1 + n.phi):n.end]
                                   sum(degp3(y, exp(lambda), exp(phi), xi, u=th, log.d=TRUE))
                                  }
                    }, # Close log.lik
                    start = function(data){
                              c(0, log(sd(data$y)), rep(1e-05, sum(sapply(data$D, ncol)) - 2))
                    }, # Close start
                    rl = function(m, param, model){
                       qegp3(1/(m * model$rate),
                             exp(param[, 1]), exp(param[, 2]), param[, 3], u=model$threshold,
                             lower.tail=FALSE)
                    },
                    endpoint = function(param, model){
                       res <- model$threshold - exp(param[, 2]) / param[, 3]
                       res[param[, 3] >= 0] <- Inf
                       res
                     },
                     delta = function(param, m, model){
                               kappa <- exp(param[1])
                               sigma <- exp(param[2])
                               xi <- param[3]

                               out <- matrix(0, nrow=3, ncol=length(m))

                               z <- 1 - 1/(m*model$rate)

                               dk <- -(z^(1/kappa) * sigma * (1 - z^(1/kappa))^(-xi - 1) * log(z)) / kappa^2
                               ds <- ((1 - z^(1/kappa))^(-xi) -1 ) / xi

                               dz1 <- -sigma * (1 - z^(1/kappa))^(-xi) * log(1 - z^(1/kappa)) / xi
                               dz2 <- -sigma * ((1 - z^(1/kappa))^(-xi) -1) / xi^2
                               dz <- dz1 + dz2

                               out[1, ] <- dk * kappa # To account for working with exp(lambda)
                               out[2, ] <- ds * sigma # To account for working with exp(phi)
                               out[3, ] <- dz
                               out
                     },
                     resid = function(o){
                       p <- texmexMakeParams(coef(o), o$data$D)
                       delta <- (o$data$y - o$threshold) / exp(p[, 2]) * p[, 3]
                       #(texmex:::.log1prel(delta * p[, 3]) * delta)^{exp(p[, 1])} # Standard exponential
                       r <- (1 - (1 + delta)^(-1/p[, 3]))^exp(p[, 1])
                       -log(1-r)
                     } # Close resid
)

