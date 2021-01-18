# Copyright (c) 2020 Chaim Even-Zohar
#
# This file is part of the R package "independence".
#
# "independence" is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# "independence" is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with "independence".  If not, see <https://www.gnu.org/licenses/>.

#' independence
#'
#' Fast Rank-Based Independence Testing
#'
#' @details
#' The package \code{independence} provides three ranking-based nonparametric tests
#' for the independence of two continuous variables X and Y:
#' \enumerate{
#' \item the classical Hoeffding's D test:
#' \code{\link{hoeffding.D.test}}
#' \item a refined variant of it, named R:
#' \code{\link{hoeffding.refined.test}}
#' \item the Bergsma-Dassios T* sign covariance:
#' \code{\link{tau.star.test}}
#' }
#' The first test is consistent assuming an absolutely continuous joint distribution,
#' i.e., the population coefficient D=0 iff the variables are independent.
#' The latter two are consistent under no restriction on the distribution.
#'
#' Given an iid sample (X1,Y1),...,(Xn,Yn),
#' all three statistics are computed in time O(n log n)
#' improving upon previous implementations.
#' The statistics R and T* are computed by a new algorithm,
#' following work of Even-Zohar and Leng.
#' It is based on the fast counting of certain patterns
#' in the permutation that relates the ranks of X and Y.
#' See [arxiv:2010.09712] and references therein.
#'
#' @seealso
#' \code{\link{tau.star.test}},
#' \code{\link{hoeffding.D.test}},
#' \code{\link{hoeffding.refined.test}}
#' \code{\link{relative.order}}
#'
#' @examples
#' library(independence)
#'
#' ## independent
#' set.seed(123)
#' xs = rnorm(10000)
#' ys = rnorm(10000)
#' hoeffding.D.test(xs,ys)
#' hoeffding.refined.test(xs,ys)
#' tau.star.test(xs,ys)
#'
#' ## dependent, even though uncorrelated
#' set.seed(123)
#' xs = rnorm(10000,0,3001:13000)
#' ys = rnorm(10000,0,3001:13000)
#' hoeffding.D.test(xs,ys)
#' hoeffding.refined.test(xs,ys)
#' tau.star.test(xs,ys)
#'
#' ## dependent but not absolutely continuous, fools Hoeffding's D
#' set.seed(123)
#' xs = runif(200)
#' f = function(x,y) ifelse(x>y, pmin(y,x/2), pmax(y,(x+1)/2))
#' ys = f(xs,runif(200))
#' hoeffding.D.test(xs,ys)
#' hoeffding.refined.test(xs,ys)
#' tau.star.test(xs,ys)
#'
#' @aliases
#' independence-test dependence-test
#'
#' @docType package
#' @name independence
#' @author Chaim Even-Zohar <\email{chaim@@ucdavis.edu}>
#' @import Rcpp
#' @importFrom Rcpp evalCpp
#' @useDynLib independence, .registration = TRUE
"_PACKAGE"
#> [1] "_PACKAGE"
