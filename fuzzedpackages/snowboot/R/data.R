#' 10 Simulated Networks of Order 2000 with Polylogarithmic (0.1, 2)
#' Degree Distributions
#'
#' A list called "artificial_networks". The length of the list is 10, and each element is a
#' network object of order 2000. These networks were simulated using the
#' polylogarithmic (aka Gutenberg--Richter law) degree distribution with parameters
#' \eqn{\delta = 0.1} and \eqn{\lambda = 2} as shown in the following equations:
#' \deqn{f(k) = k^{-{\delta}}e^{-{k/{\lambda}}}/Li_{\delta}(e^{-{1/\lambda}})}{f(k)=k^-\delta exp(-k/\lambda )/Li[\delta](exp(-1/\lambda)),}
#' \deqn{Li_{\delta}(z)=\sum_{j=1}^{\infty} z^{-j}/{j^{\delta}},}{Li[\delta](z)=\sum_{j=1}^{\infty} z^{-j}/{j^{\delta}},}
#' where \eqn{\lambda > 0}
#' \insertCite{@see @newman_etal_2001, @gel_etal_2017, and @chen_etal_2018_snowboot for details}{snowboot}.
#'
#' @format A list containing 10 network objects. Each network object is a list
#' with three elements:
#' \describe{
#'   \item{\code{degree}}{the degree sequence of the network, which is
#'      an integer vector of length \eqn{n};}
#'   \item{\code{edges}}{the edgelist, which is a two-column
#'      matrix, where each row is an edge of the network;}
#'   \item{\code{n}}{the network order (number of nodes in the network).
#'   The order is 2000.}
#' }
#'
#' @references
#' \insertAllCited{}

"artificial_networks"
