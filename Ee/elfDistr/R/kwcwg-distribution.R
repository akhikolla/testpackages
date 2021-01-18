

#' Kumaraswamy Complementary Weibull Geometric Probability Distribution
#'
#' Density, distribution function, quantile function and random generation
#' for the Kumaraswamy Complementary Weibull Geometric (Kw-CWG) probability
#' distribution.
#'
#' @param x,q	          vector of quantiles.
#' @param p	              vector of probabilities.
#' @param n	              number of observations. If \code{length(n) > 1},
#'                        the length is taken to be the number required.
#' @param alpha,beta,gamma,a,b	      Parameters of the distribution. 0 < alpha < 1, and the other parameters mustb e positive.
#' @param log,log.p	      logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail	    logical; if TRUE (default), probabilities are \eqn{P[X \le x]}
#'                        otherwise, \eqn{P[X > x]}.
#'
#' @details
#'
#' Probability density function
#' \deqn{
#'    f(x) = \alpha^a \beta \gamma a b (\gamma x)^{\beta - 1} \exp[-(\gamma x)^\beta] \cdot
#'    \frac{\{1 - \exp[-(\gamma x)^\beta]\}^{a-1}}{\{ \alpha + (1 - \alpha) \exp[-(\gamma x)^\beta] \}^{a+1}} \cdot
#' }
#' \deqn{
#'    \cdot \bigg\{ 1 - \frac{\alpha^a[1 - \exp[-(\gamma x)^\beta]]^a}{\{ \alpha + (1 - \alpha) \exp[-(\gamma x)^\beta] \}^a} \bigg\}
#' }
#'
#' Cumulative density function
#' \deqn{
#'     F(x) = 1 - \bigg\{  1 - \bigg[ \frac{\alpha (1 - \exp[-(\gamma x)^\beta]) }{ \alpha + (1 - \alpha) \exp[-(\gamma x)^\beta] } \bigg]^a \bigg\}^b
#' }
#'
#' Quantile function
#' \deqn{
#'    Q(u) = \gamma^{-1} \bigg\{
#'        \log\bigg[\frac{
#'            \alpha + (1 - \alpha) \sqrt[a]{1 - \sqrt[b]{1 - u} }
#'        }{
#'            \alpha (1 - \sqrt[a]{1 - \sqrt[b]{1 - u} } )
#'        }\bigg]
#'    \bigg\}^{1/\beta}, 0 < u < 1
#' }
#'
#' @references
#' Afify, A.Z., Cordeiro, G.M., Butt, N.S., Ortega, E.M. and
#' Suzuki, A.K. (2017). A new lifetime model with variable shapes for
#' the hazard rate. Brazilian Journal of Probability and Statistics
#' 
#' @name Kw-CWG
#' @aliases Kw-CWG
#' @aliases kwcwg
#'
#' @keywords distribution
#' @keywords univar
#' @keywords models
#' @keywords survival
#' @concept Univariate
#' @concept Continuous
#' @concept Lifetime
#'
#' @export

dkwcwg <- function(x, alpha, beta, gamma, a, b, log = FALSE) {
  cpp_dkwcwg(x, alpha, beta, gamma, a, b, log[1L])
}


#' @rdname Kw-CWG
#' @export

pkwcwg <- function(q, alpha, beta, gamma, a, b, lower.tail = TRUE, log.p = FALSE) {
  cpp_pkwcwg(q, alpha, beta, gamma, a, b, lower.tail[1L], log.p[1L])
}


#' @rdname Kw-CWG
#' @export

qkwcwg <- function(p, alpha, beta, gamma, a, b, lower.tail = TRUE, log.p = FALSE) {
  cpp_qkwcwg(p, alpha, beta, gamma, a, b, lower.tail[1L], log.p[1L])
}


#' @rdname Kw-CWG
#' @export

rkwcwg <- function(n, alpha, beta, gamma, a, b) {
  if (length(n) > 1) n <- length(n)
  cpp_rkwcwg(n, alpha, beta, gamma, a, b)
}

