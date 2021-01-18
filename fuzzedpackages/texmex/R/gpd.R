#' @include texmexFamily.R
#' @include gpd.info.R
#' @include gpd.sandwich.R
#' @export gpd
#' @aliases evm
gpd <- texmexFamily(name = 'GPD',
                   log.lik = function(data, th, ...) {
                               y <- data$y
                               X.phi <- data$D$phi
                               X.xi <- data$D$xi
                               n.phi <- ncol(X.phi)
                               n.end <- n.phi + ncol(X.xi)
                               function(param) {
                                 stopifnot(length(param) == n.end)
                                 phi <- X.phi %*% param[1:n.phi]
                                 xi <- X.xi %*% param[(1 + n.phi):n.end]
                                 sum(dgpd(y, exp(phi), xi, u=th, log.d=TRUE))
                               }
                   }, # Close log.lik
                   param = c(phi=0, xi=0),
                   info = gpd.info,
                   sandwich = gpd.sandwich,
                   start = function(data){
                             y <- data$y
                             X.phi <- data$D[[1]]
                             X.xi <- data$D[[2]]
                             c(log(mean(y)), rep(1e-05, -1 + ncol(X.phi) + ncol(X.xi)))
                    }, # Close start

                    resid = function(o){
                              p <- texmexMakeParams(coef(o), o$data$D)
                              delta <- (o$data$y - o$threshold) / exp(p[,1])
                              .log1prel(delta * p[,2]) * delta # Standard exponential
                    }, # Close resid

                    endpoint = function(param, model){
                        res <- model$threshold - exp(param[, 1]) / param[, 2]
                        res[param[, 2] >= 0] <- Inf
                        res
                    },
                    rl = function(m, param, model){
                      ## write in terms of qgpd; let's not reinvent the wheel
                      qgpd(1/(m * model$rate),
                           exp(param[,1]), param[,2], u=model$threshold,
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
                             out[1,] <- exp(param[2]) / param[3] * ((m * param[1])^param[3] - 1)
                             out[2,] <- -exp(param[2]) / (param[3] * param[3]) * ( (m * param[1] )^param[3] - 1 ) +
                                        exp(param[2]) / param[3] * (m * param[1])^param[3] * log(m * param[1])
                             }
                             out
                    }, # Close delta
                    density = function(n, param, model, log.d=FALSE){
                                dgpd(n, exp(c(param[, 1])), c(param[, 2]), u=model$threshold, log.d=log.d)
                    },

                    rng = function(n, param, model){
                            rgpd(n, exp(c(param[, 1])), c(param[, 2]), u=model$threshold)
                    },
                    prob = function(n, param, model){
                             pgpd(n, exp(c(param[, 1])), c(param[, 2]), u=model$threshold)
                    },
                    quant = function(n, param, model){
                              qgpd(n, exp(c(param[, 1])), c(param[, 2]), u=model$threshold)
                    }
)

#' @export gpdIntCensored
#' @aliases evm
gpdIntCensored <- texmexFamily(name = 'gpdIntCensored',
							   log.lik= function (data, th, dp=2, ...) 
							   {
							   	y <- data$y
							   	X.phi <- data$D$phi
							   	X.xi <- data$D$xi
							   	n.phi <- ncol(X.phi)
							   	n.end <- n.phi + ncol(X.xi)
							   	function(param) {
							   		stopifnot(length(param) == n.end)
							   		phi <- X.phi %*% param[1:n.phi]
							   		xi <- X.xi %*% param[(1 + n.phi):n.end]
							   		sum(log(pgpd(y+(10^(-dp))/2,exp(phi),xi,u=th) - pgpd(y-(10^(-dp))/2,exp(phi),xi,u=th)))
							   	}
							   },
							   start=gpd$start,
							   resid=gpd$resid,
							   param=gpd$param,
							   rl=gpd$rl,
							   delta=gpd$delta,
							   endpoint=gpd$endpoint,
							   density=gpd$density,
							   rng=gpd$rng,
							   prob=gpd$prob,
							   quant=gpd$quant)

