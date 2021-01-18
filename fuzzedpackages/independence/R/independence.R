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


#' The relative order of two vectors
#'
#' Given (X1,Y1),...,(Xn,Yn), many nonparametric statistics depend only
#' on the permutation P that satisfies rank Yi = P[rank Xi].
#' The function \code{relative_order} computes such P given X1,...,Xn and Y1,...,Yn.
#'
#' By default, the function removes missing values, and warns of
#' repeating values.
#' Then it computes the relative order by calling the base R function
#' \code{\link{order}} twice: \code{order(xs[order(ys)])}.
#'
#' Ties may be broken arbitrarily, depending on the behavior of the function
#' \code{order}.
#'
#' @param xs,ys Numeric vectors of same length.
#' @param na.rm Logical: Should missing values, \code{NaN}, and \code{Inf} be removed?
#' @param collisions Logical: Warn of repeating values in \code{xs} or \code{ys}.
#'
#' @return
#' An integer vector which describes the ordering of the second argument
#' \code{ys}, in terms of the ordering of the corresponding values in the
#' first argument \code{xs}.
#'
#' For example, if \code{xs[3]} is the \code{i}th smallest
#' and \code{ys[3]} is the \code{j}th smallest,
#' then the returned value in position \code{i} is \code{j}.
#'
#' @export
#' @examples
#' relative.order(1:5, c(10,30,50,40,20))
#' ## [1] 1 3 5 4 2
#'
#' relative.order(c(1,2,5,3,4), c(10,30,50,40,20))
#' ## [1] 1 3 4 2 5
#'
#' set.seed(123)
#' relative.order(runif(8), runif(8))
#' ## [1] 5 4 8 1 3 2 7 6

relative.order <- function(xs, ys, na.rm = TRUE, collisions = TRUE)
{
  stopifnot(is.vector(xs), is.vector(ys), is.numeric(xs), is.numeric(ys),
            length(xs) == length(ys), is.logical(na.rm), is.logical(collisions))
  if (na.rm)
  {
    stopifnot(any(good <- is.finite(xs) & is.finite(ys)))
    xs <- xs[good]
    ys <- ys[good]
  }
  if (collisions)
  {
    dx = sum(duplicated(xs))
    dy = sum(duplicated(ys))
    if (dx) warning(sprintf("%d collisions in 1st variable", dx))
    if (dy) warning(sprintf("%d collisions in 2nd variable", dy))
  }
  return(order(xs[order(ys)]))
}


#' Generic rank test for paired samples
#'
#' An internal function unifying several nonparametric tests for
#' paired samples.
#'
#' The function \code{.generic.rank.test} first calls
#' \code{relative.ordering} with \code{xs} and \code{ys}.
#' Then it uses the given function to compute the test statistic
#' from the resulting permutation.
#' The statistic is rescaled by multiplication with
#' \code{(n-1)*limit_law_coef}, where \code{n} is the sample size.
#' Finally, it computes the p-value by calling
#' \code{\link[TauStar]{pHoeffInd}} from the package \code{TauStar}.
#'
#' @inheritParams relative.order
#'
#' @param xs,ys Same-length numeric vectors, containing paired samples.
#' @param test Function computing the test statistic given a relative order.
#' @param letter Notation for the test statistic, e.g., "D" for Hoeffding's D.
#' @param description Full name of test.
#' @param limit_law_coef Scaling of test statistic for standard null distribution.
#' @param min_samples,max_samples Data size limits.
#' @param precision of p-value, between 0 and 1. Otherwise p-value=\code{NA}.
#'
#' @return A list, of class \code{"indtest"}:
#' \tabular{ll}{
#' \code{method} \tab
#' the test's name \cr
#' \code{n} \tab
#' number of data points used \cr
#' \code{Tn}/\code{Dn}/\code{Rn}/\dots{} \tab
#' the test statistic, measure of dependence \cr
#' \code{scaled} \tab
#' the test statistic rescaled for a standard null distribution \cr
#' \code{p.value} \tab
#' the asymptotic p-value, by TauStar::\code{\link[TauStar]{pHoeffInd}} \cr}
#'
#' @seealso
#' \code{\link{independence}},
#' \code{\link{relative.order}},
#' \code{\link{tau.star.test}},
#' \code{\link{hoeffding.D.test}},
#' \code{\link{hoeffding.refined.test}}
#'
#' @section P-value:
#' The null distribution of the test statistic was described by Hoeffding.
#' The p-value is approximated by calling the function
#' \code{\link[TauStar]{pHoeffInd}} from the package \code{TauStar} by
#' Luca Weihs.
#'
#' By default, the p-value's \code{precision} parameter is set to \code{1e-5}.
#' It seems that better precision would cost a considerable amount of time,
#' especially for large values of the test statistic.
#' It is therefore recommended to modify this parameter only upon need.
#'
#' In case that \code{TauStar} is unavailable, or to save time in repeated use,
#' set \code{precision = 1} to avoid computing p-values altogether.
#' The \code{scaled} test statistic may be used instead.
#' Its asymptotic distribution does not depend on any parameter.
#' Also the raw test statistic may be used, descriptively,
#' as a measure of dependence.
#' Only its accuracy depends on the sample size.
#'
#' @section Ties:
#' This package currently assumes that the variables under consideration are
#' non-atomic, so that ties are not expected, other than by occasional effects of
#' numerical precision. Addressing ties rigorously is left for future versions.
#'
#' The flag \code{collisions = TRUE} invokes checking for ties in \code{xs}
#' and in \code{ys}, and produces an appropriate warning if they exist.
#' The current implementation breaks such ties arbitrarily, not randomly.
#'
#' By the averaging nature of the test statistic,
#' it seems that a handful of ties should not be of much concern.
#' In case of more than a handful of ties, our current advice
#' to the user is to break them uniformly at random beforehand.
#'
#' @section Big Data:
#' The test statistic is computed in almost linear time, O(n log n),
#' given a sample of size n. Its computation involves integer arithmetics
#' of order n^4 or n^5, which should fit into an integer data type
#' supported by the compiler.
#'
#' Most 64-bit compilers emulate 128-bit arithmetics.
#' Otherwise we use the standard 64-bit arithmetics.
#' Find the upper limits of your environment using
#' \enumerate{
#' \item \code{\link[independence]{max_taustar}()}
#' \item \code{\link[independence]{max_hoeffding}()}}
#'
#' Another limitation is 2^31-1, the maximum size and value of
#' an integer vector in a 32-bit build of R.
#' This is only relevant for the tau star statistic in 128-bit mode,
#' which could otherwise afford about three times that size.
#' If your sample size falls in this range, try recompiling the
#' function \code{.\link{calc.taustar}}
#' according to the instructions in the cpp source file.
#'
#' @export

.generic.rank.test <- function(xs, ys, test, letter, description,
                               na.rm = TRUE,
                               collisions = TRUE,
                               precision = 1e-5,
                               limit_law_coef = 1.0,
                               min_samples = 1,
                               max_samples = Inf)
{
  v = list(method = description)
  class(v) <- "indtest"

  ordering = relative.order(xs, ys, na.rm = na.rm, collisions = collisions) - 1
  v$n <- n <- length(ordering)

  if (n < min_samples) stop(sprintf(
    "too few samples: %d are given, at least %d are needed", n, min_samples))
  if (n > max_samples) stop(sprintf(
    "too many samples: %d are given, compiler limit is %d", n, max_samples))

  v[[sprintf("%sn", letter)]] <- Tn <- test(ordering)
  if (Tn == -1) stop("computation failed")
  v[["scaled"]] <- scaled <- Tn*(n-1)*limit_law_coef

  if ((precision > 0) && (precision < 1) && requireNamespace("TauStar"))
    p.value <- TauStar::pHoeffInd(scaled, lower.tail = FALSE, precision)[1]
  else p.value <- NA
  v$p.value <- p.value

  return(v)
}


#' The Bergsma--Dassios tau* independence test
#'
#' The function \code{tau.star.test} provides an
#' independence test for two continuous numeric variables,
#' that is consistent for all alternatives.
#' It is based on the Yanagimoto-Bergsma-Dassios coefficient,
#' which compares the frequencies of concordant and discordant
#' quadruples of data points.
#' The test statistic is efficiently computed in O(n log n) time,
#' following work of Even-Zohar and Leng.
#'
#' @param xs,ys Same-length numeric vectors, containing paired samples.
#' @inherit .generic.rank.test return
#' @inheritParams .generic.rank.test
#' @inheritParams relative.order
#' @inheritSection .generic.rank.test P-value
#' @inheritSection .generic.rank.test Ties
#' @inheritSection .generic.rank.test Big Data
#' @aliases taustar Bergsma-Dassios Yanagimoto-Bergsma-Dassios
#'
#' @seealso
#' \code{\link{independence}},
#' \code{\link{hoeffding.D.test}},
#' \code{\link{hoeffding.refined.test}}
#'
#' @examples
#' ## independent, $p.value is 0.2585027
#' set.seed(123)
#' tau.star.test(rnorm(10000),rnorm(10000))
#'
#' ## dependent, even though uncorrelated, $p.value is 0.000297492
#' set.seed(123)
#' tau.star.test(rnorm(10000,0,3001:13000), rnorm(10000,0,3001:13000))
#'
#' @references
#' Bergsma, Wicher; Dassios, Angelos. A consistent test of independence based
#' on a sign covariance related to Kendall's tau. \emph{Bernoulli} 20 (2014), no.
#' 2, 1006--1028.
#' \cr\cr
#' Yanagimoto, Takemi. "On measures of association and a related problem."
#' Annals of the Institute of Statistical Mathematics 22.1 (1970): 57-63.
#' \cr\cr
#' Luca Weihs (2019). TauStar: Efficient Computation and Testing of the
#' Bergsma-Dassios Sign Covariance. R package version 1.1.4.
#' https://CRAN.R-project.org/package=TauStar
#' \cr\cr
#' Even-Zohar, Chaim, and Calvin Leng. "Counting Small Permutation Patterns."
#' arXiv preprint arXiv:1911.01414 (2019).
#' \cr\cr
#' Even-Zohar, Chaim. "independence: Fast Rank Tests."
#' arXiv preprint arXiv:2010.09712 (2020).
#'
#' @export
tau.star.test <- function(xs, ys, na.rm = TRUE, collisions = TRUE, precision = 1e-5)
{
  return(.generic.rank.test(xs, ys, .calc.taustar, "T",
                            "Bergsma-Dassios t* Independence Test",
                            na.rm = na.rm,
                            collisions = collisions,
                            precision = precision,
                            limit_law_coef = 1.0,
                            min_samples = 4,
                            max_samples = max_taustar()))
}


#' Hoeffding's D independence test
#'
#' The function \code{hoeffding.D.test} provides
#' independence testing for two continuous numeric variables,
#' that is consistent for absolutely-continuous alternative
#' bivariate distributions.
#' It implements the classical D statistic by Hoeffding,
#' which in terms of CDFs estimates the integral of (Fxy-Fx*Fy)^2 dFxy.
#' It may also be expressed in terms of the ordering types of quintuples
#' of data points.
#' Its efficient O(n log n) computation seems to be new in R.
#'
#' @param xs,ys Same-length numeric vectors, containing paired samples.
#' @inherit .generic.rank.test return
#' @inheritParams .generic.rank.test
#' @inheritParams relative.order
#' @inheritSection .generic.rank.test P-value
#' @inheritSection .generic.rank.test Ties
#' @inheritSection .generic.rank.test Big Data
#' @aliases Hoeffding-test
#'
#' @seealso
#' \code{\link{independence}},
#' \code{\link{tau.star.test}},
#' \code{\link{hoeffding.refined.test}}
#'
#' @examples
#' ## independent, $p.value is 0.2582363
#' set.seed(123)
#' hoeffding.D.test(rnorm(10000),rnorm(10000))
#'
#' ## dependent, even though uncorrelated, $p.value is 0.0002891223
#' set.seed(123)
#' hoeffding.D.test(rnorm(10000,0,3001:13000), rnorm(10000,0,3001:13000))
#'
#' @references
#' Hoeffding, Wassily. "A non-parametric test of independence."
#' The annals of mathematical statistics (1948): 546-557.
#' \cr\cr
#' Luca Weihs (2019). TauStar: Efficient Computation and Testing of the
#' Bergsma-Dassios Sign Covariance. R package version 1.1.4.
#' https://CRAN.R-project.org/package=TauStar
#' \cr\cr
#' Frank E Harrell Jr, with contributions from Charles Dupont and many
#' others. (2020). Hmisc: Harrell Miscellaneous. R package version
#' 4.4-0. https://CRAN.R-project.org/package=Hmisc
#' \cr\cr
#' Even-Zohar, Chaim. "independence: Fast Rank Tests."
#' arXiv preprint arXiv:2010.09712 (2020).
#'
#' @export
hoeffding.D.test <- function(xs, ys, na.rm = TRUE, collisions = TRUE, precision = 1e-5)
{
  return(.generic.rank.test(xs, ys, .calc.hoeffding, "D",
                            "Hoeffding's D Independence Test",
                            na.rm = na.rm,
                            collisions = collisions,
                            precision = precision,
                            limit_law_coef = 36.0,
                            min_samples = 5,
                            max_samples = max_hoeffding()))
}


#' The refined Hoeffding independence test
#'
#' The function \code{hoeffding.refined.test} provides
#' independence testing for two continuous numeric variables,
#' that is consistent for all alternatives.
#' The test statistic is a variant of the classical Hoeffding's D statistic.
#' In terms of CDFs, it estimates the integral of (Fxy-Fx*Fy)^2 dFx dFy,
#' based on the ordering types of quintuples of data points.
#' This test statistic is efficiently computed via a new O(n log n)-time
#' algorithm, following work of Even-Zohar and Leng.
#'
#' @param xs,ys Same-length numeric vectors, containing paired samples.
#' @inherit .generic.rank.test return
#' @inheritParams .generic.rank.test
#' @inheritParams relative.order
#' @inheritSection .generic.rank.test P-value
#' @inheritSection .generic.rank.test Ties
#' @inheritSection .generic.rank.test Big Data
#' @aliases refined-Hoeffding-test
#'
#' @seealso
#' \code{\link{independence}},
#' \code{\link{tau.star.test}},
#' \code{\link{hoeffding.D.test}},
#'
#' @examples
#' ## independent, $p.value is 0.258636
#' set.seed(123)
#' hoeffding.refined.test(rnorm(10000),rnorm(10000))
#'
#' ## dependent, even though uncorrelated, $p.value is 0.0003017679
#' set.seed(123)
#' hoeffding.refined.test(rnorm(10000,0,3001:13000), rnorm(10000,0,3001:13000))
#'
#' @references
#' Hoeffding, Wassily. "A non-parametric test of independence."
#' The annals of mathematical statistics (1948): 546-557.
#' \cr\cr
#' Yanagimoto, Takemi. "On measures of association and a related problem."
#' Annals of the Institute of Statistical Mathematics 22.1 (1970): 57-63.
#' \cr\cr
#' Luca Weihs (2019). TauStar: Efficient Computation and Testing of the
#' Bergsma-Dassios Sign Covariance. R package version 1.1.4.
#' https://CRAN.R-project.org/package=TauStar
#' \cr\cr
#' Even-Zohar, Chaim. "Patterns in Random Permutations."
#' arXiv preprint arXiv:1811.07883 (2018).
#' \cr\cr
#' Even-Zohar, Chaim, and Calvin Leng. "Counting Small Permutation Patterns."
#' arXiv preprint arXiv:1911.01414 (2019).
#' \cr\cr
#' Even-Zohar, Chaim. "independence: Fast Rank Tests."
#' arXiv preprint arXiv:2010.09712 (2020).
#'
#' @export
hoeffding.refined.test <- function(xs, ys, na.rm = TRUE, collisions = TRUE, precision = 1e-5)
{
  return(.generic.rank.test(xs, ys, .calc.refined, "R",
                            "Refined Hoeffding Independence Test",
                            na.rm = na.rm,
                            collisions = collisions,
                            precision = precision,
                            limit_law_coef = 36.0,
                            min_samples = 5,
                            max_samples = max_hoeffding()))
}
