asymptotic_distributions <- list(
  MOI = function(k) {
    sigma_k <- choose(2*k, k)^2 * (sqrt(pi) * gamma(k+1.5)*gamma(2*k+1) +
                                   2*gamma(2*k + 1.5) * (gamma(k+1.5) -
                                           sqrt(pi)*gamma(k+1))) /
      (2^(4*k-1)*(2*k+1)^2*gamma(k+1.5)*gamma(2*k+1.5))
    sigma <- sqrt((2 * k + 1)^2 * sigma_k)
    function(q) pnorm(q, 0, sigma)
  },
  M = pnorm,
  CM = pnorm,
  MGG = pnorm,
  WCX = function(x) pnorm(x, sd = sqrt(1/3)),
  SGN = function(x) pnorm(x, sd = 1/2),
  BHI = function(q) pnorm(q, sd = sqrt(1/20)),
  NAI = function(k) {
    k <- k + 1
    sigma_k <- (2^(3 - 2*k)*(2*k*(-1 + 2*k)*sqrt(pi)*gamma(k) +
                               (2^(3 + k) + (-8 - 2^(4 + k) + 4^k)*k +
                                  16*k^2)*gamma(0.5 + k)))/
      (k^3*(-1 + 2*k)*gamma(0.5 + k))
    sigma <- sqrt(k^2 * sigma_k)
    function(q) pnorm(q, 0, sigma)
  }
)
