#' @useDynLib fastpos
#' @importFrom Rcpp sourceCpp
NULL

#' @details
#' In most cases you will just need the function \code{\link{find_critical_pos}}
#' which will you give you the critical point of stability for your specific
#' parameters. If you are interested in more complicated analysis you might want
#' to look at the function \code{\link{simulate_pos}}, which is a C++
#' functions to calculate correlations and return points of stability.
#'
#' @references Schönbrodt, F. D. & Perugini, M. (2013). At what sample size do
#' correlations stabilize? \emph{Journal of Research in Personality, 47},
#' 609-612. \url{https://doi.org/10.1016/j.jrp.2013.05.009}
#'
#' Schönbrodt, F. D. & Perugini, M. (2018) Corrigendum to “At what sample size
#' do correlations stabilize?” [J. Res. Pers. 47 (2013) 609–612.
#' https://doi.org/10.1016/j.jrp.2013.05.009]. \emph{Journal of Research in
#' Personality, 74}, 194. \url{https://doi.org/10.1016/j.jrp.2018.02.010}
#' @keywords internal
"_PACKAGE"

#' Stops execution without giving error
#'
#' This is useful to have a consistent behavior when the user interrupts
#' function execution but this interruption is not catched by C++. If this
#' happens nothing is returned. But if C++ catches the interrupt, we need to
#' stop execution ourselves (and also return nothing).
#' @noRd
stop_quietly <- function() {
  blankMsg <- sprintf("\r%s\r", paste(rep("", getOption("width")-1L),
                                      collapse=" "));
  stop(simpleError(blankMsg));
}

#' Creates a population with a specified correlation.
#'
#' @param rho Population correlation.
#' @param size Population size.
#' @return Two-dimensional population matrix with a specific correlation.
#' @examples
#' pop <- create_pop(0.5, 100000)
#' cor(pop)
#' @noRd
create_pop_inexact <- function(rho, size){
  mu <- c(1, 2)
  s1 <- 2
  s2 <- 8
  sigma <- matrix(c(s1^2, s1 * s2 * rho, s1 * s2 * rho, s2^2), 2)
  pop <- MASS::mvrnorm(n = size, mu = mu, Sigma = sigma)
  pop
}

#' Creates a population with a specified correlation.
#'
#' The correlation will be exactly the one specified. The used method is
#' described here:
#' https://stats.stackexchange.com/questions/15011/generate-a-random-variable-with-a-defined-correlation-to-an-existing-variables/15040#15040
#'
#' @param rho Population correlation.
#' @param size Population size.
#' @return Two-dimensional population matrix with a specific correlation.
#' @examples
#' pop <- create_pop(0.5, 100000)
#' cor(pop)
#' @export
create_pop <- function(rho, size) {
  y <- stats::rnorm(size)
  x <- stats::rnorm(size)
  y.perp <- stats::residuals(stats::lm.fit(cbind(1, y), x))
  x <- rho * stats::sd(y.perp) * y + y.perp * stats::sd(y) * sqrt(1 - rho^2)
  matrix(c(x, y), ncol = 2)
}

#' Run simulation for one specific correlation.
#'
#' @param rho Population correlation.
#' @inheritParams find_critical_pos
#' @return A list with two elements, (1) a data frame called "summary"
#'   containing all the above information as well as the critical sample sizes
#'   (points of stability) for the confidence-levels specified and (2) vector
#'   "n" with the sample size from each study (e.g. for plotting the
#'   distribution)
#' @examples
#' find_one_critical_pos(rho = 0.5)
#' @noRd
find_one_critical_pos <- function(rho, sample_size_min = 20,
                                  sample_size_max = 1000,
                                  replace = TRUE, n_studies = 1000,
                                  pop_size = 1e6, precision = .1,
                                  precision_rel = T,
                                  confidence_levels = c(.8, .9, .95),
                                  n_cores = 1) {
  # create corridor of stability
  if (precision_rel) {
    limits <- rho * (1 + c(-1, 1) * precision)
  } else {
    limits <- rho + c(-1, 1) * precision
  }
  lower_limit <- limits[1]
  upper_limit <- limits[2]
  # create bivariate population distribution
  pop <- create_pop(rho, pop_size)

  x <- pop[,1]
  y <- pop[,2]
  rho_pop <- stats::cor(x, y)

  # create dist of pos

  # we will use 1 normal invocation and n_cores - 1 futures, this way we have
  # a progress bar
  if (n_cores > 1){
    f <- list()
    for (ii in seq(n_cores - 1)) {
      f[[ii]] <- future::future({
        simulate_pos(x, y, ceiling(n_studies/(n_cores)), sample_size_min,
                     sample_size_max, replace, lower_limit, upper_limit)
      }, seed = TRUE)
    }
    res <- simulate_pos(x, y, ceiling(n_studies/n_cores), sample_size_min,
                        sample_size_max, replace, lower_limit, upper_limit)
    v <- unlist(lapply(f, FUN = future::value))
    res <- c(res, v)
  } else {
    res <- simulate_pos(x, y, n_studies, sample_size_min, sample_size_max, T,
                        lower_limit, upper_limit)
  }

  # on interruption, C++ will return -1 (if R interrupts by itself, nothing
  # will be returned, it just stops)
  if (length(res) == 1) {
    if (res == -1) stop_quietly()
  }
  names(res) <- unlist(paste("study ", 1:length(res)))
  n_not_breached <- sum(is.na(res))

  # calc critical pos
  # use max sample size for those studies where the corridor was not breached
  thequantiles <- stats::quantile(c(res, rep(sample_size_max, n_not_breached)),
                                    confidence_levels, na.rm = T)
  return(list(summary = c(rho_pop = rho_pop, thequantiles,
                          sample_size_min = sample_size_min,
                          sample_size_max = sample_size_max,
                          lower_limit=lower_limit,
                          upper_limit=upper_limit,
                          n_studies = length(res),
                          n_not_breached = n_not_breached),
              n = res
              ))
}

#' Find the critical point of stability
#'
#' Run simulations for one or several population correlations and return the
#' critical points of stability (POS). The critical point of stability is the
#' sample size at which a certain percentage of studies will fall into an a
#' priori specified interval and stay in this interval if the sample size is
#' increased further.
#'
#' @param rhos Vector of population correlations (can also be a single
#'   correlation).
#' @param precision Precision around the correlation which is acceptable
#'   (defaults to 0.1). The precision will determine the corridor of stability
#'   which is just rho+-precision.
#' @param precision_rel Whether the precision is absolute (rho+-precision or
#'   relative rho+-rho*precision), boolean (defaults to FALSE).
#' @param n_studies Number of studies to run for each rho (defaults to 10e3).
#' @param sample_size_min Minimum sample size for each study (defaults to 20).
#' @param sample_size_max Maximum sample size for each study (defaults to 1e3).
#' @param replace Whether drawing samples is with replacement or not.
#' @param confidence_levels Confidence levels for point of stability. This
#'   corresponds to the quantile of the distribution of all found critical
#'   sample sizes (defaults to c(.8, .9, .95)).
#' @param pop_size Population size (defaults to 1e6).
#' @param n_cores Number of cores to use for simulation.
#' @return A data frame containing all the above information, as well as the
#'   points of stability.
#' @examples
#' find_critical_pos(rhos = 0.5)
#' find_critical_pos(rhos = c(0.4, 0.5), n_studies = 1e3)
#' @export
find_critical_pos <- function(rhos, precision = 0.1, precision_rel = FALSE,
                              sample_size_min = 20, sample_size_max = 1000,
                              replace = TRUE,
                              n_studies = 10000,
                              confidence_levels = c(.8, .9, .95),
                              pop_size = 1e6,
                              n_cores = 1) {
  defaultplan <- options("future.plan")[[1]]
  options("future.plan" = "multisession")
  result <- mapply(find_one_critical_pos, rhos,
                   sample_size_max = sample_size_max,
                   sample_size_min = sample_size_min,
                   n_studies = n_studies,
                   precision = precision,
                   precision_rel = precision_rel,
                   n_cores = n_cores,
                   MoreArgs = list(confidence_levels = confidence_levels),
                   SIMPLIFY = F)
  summary <- lapply(result, function(x) x[[1]])
  summary <- plyr::ldply(summary)
  summary <- cbind(summary, precision, precision_rel)
  sum_n_not_breached = sum(summary$n_not_breached)
  if (sum_n_not_breached > 0){
    warning(sum_n_not_breached, " simulation[s] did not reach the corridor of
            stability", ".\nIncrease sample_size_max and rerun the simulation.",
            sep = "")
  }
  options("future.plan" = defaultplan)
  summary
}

.onUnload <- function(libpath) {
  library.dynam.unload("fastpos", libpath)
}