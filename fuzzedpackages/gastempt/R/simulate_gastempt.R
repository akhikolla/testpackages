#' Simulate gastric emptying data following a linexp or powexp function
#'
#' @param n_records Number of records
#' @param v0_mean,v0_std Mean and between record standard deviation of initial volume, typically in ml.
#' @param tempt_mean,tempt_std Mean and between record standard deviation of parameter $t_{empt}$, typically in minutes.
#' @param kappa_mean,kappa_std For linexp only: Mean and between-record standard deviation of overshoot parameter \code{kappa}. For values of \code{kappa} above 1, curve has an overshoot that can be used to follow volume time series with secretion.
#' @param beta_mean,beta_std For powexp only: Mean and between-record standard deviation of the so called lag parameter.
#' @param noise Standard deviation of normal noise when \code{student_t_df = NULL}; scaling of noise when student_t_df >= 2.
#' @param student_t_df When NULL (default), Gaussian noise is added; when >= 2, Student_t distributed noise is added, which generates more realistic outliers. Values from 2 to 5 are useful, when higher values are used the result comes close to that of Gaussian noise. Values below 2 are rounded to 2.
#' @param missing When 0 (default), all curves have the same number of data points. When > 0, this is the fraction of points that were removed randomly to simulate missing points. Maximum value is 0.5.
#' @param model linexp(default) or powexp
#' @param seed optional seed; not set if seed = NULL (default)
#' @param max_minute Maximal time in minutes; if NULL, a sensible
#' default rounded to hours is used
#'
#' @return A list with 3 elements:
#' \describe{
#'   \item{record}{Data frame with columns
#'     \code{record(chr), v0, tempt, kappa/beta} giving the effective
#'     \code{linexp} or \code{powexp} parameters for the individual record.
#'     \code{v0} is rounded to nearest integer.}
#'   \item{data}{Data frame with columns
#'     \code{record(chr), minute(dbl), vol(dbl)} giving the
#'      time series and grouping parameters. \code{vol} is rounded
#'      to nearest integer.}
#'   \item{stan_data}{A list for use as \code{data} in Stan-based fits
#'   with elements \code{prior_v0, n, n_record, record, minute, volume}.}
#'  }
#'  A comment is attached to the return value that can be used as a title
#' @examples
#' suppressWarnings(RNGversion("3.5.0"))
#' set.seed(4711)
#' library(ggplot2)
#' vol_linexp = simulate_gastempt(n_records = 4, noise = 20)
#' ggplot(vol_linexp$data, aes(x = minute, y = vol)) + geom_point() +
#'   facet_wrap(~record) + ggtitle("linexp, noise = 0, no missing")
#'
#' vol_powexp = simulate_gastempt(n_records = 4, missing = 0.2, student_t_df = 2)
#' ggplot(vol_powexp$data, aes(x = minute, y = vol)) + geom_point() +
#'   facet_wrap(~record) + ggtitle("powexp, noise = 10 (default), 20% missing,
#'   Student-t (df = 2) noise")
#' @import dplyr
#' @import assertthat
#' @export
#'

simulate_gastempt = function(
  n_records = 10,
  v0_mean = 400, v0_std = 50,
  tempt_mean = ifelse(identical(model, linexp), 60, 120),
  tempt_std = tempt_mean/3,
  kappa_mean = 0.7, # linexp only
  kappa_std = kappa_mean/3,
  beta_mean = 0.7, # powexp only
  beta_std = beta_mean/3,
  noise = 20,
  student_t_df = NULL,
  missing = 0,
  model = linexp,
  seed = NULL,
  max_minute = NULL) {

  # Hack to avoid notes
  vol = . = NULL
  if (!is.null(seed)) {
    suppressWarnings(RNGversion("3.5.0"))
    set.seed(seed)
  }
  # Only linexp and powexp are supported
  assert_that(identical(model, linexp) || identical(model, powexp))
  model_name = ifelse(identical(model, linexp), "linexp", "powexp")
  assert_that(v0_std >= 0)
  assert_that(tempt_std >= 0)
  assert_that(kappa_std >= 0)
  assert_that(beta_std >= 0)

  assert_that(v0_mean > 2*v0_std)
  assert_that(tempt_mean > 2*tempt_std)
  assert_that(kappa_mean > 2*kappa_std)
  assert_that(beta_mean > 2*beta_std)


  # Time range as 3* t50
  if (is.null(max_minute)) {
    if (model_name == "linexp") {
      max_minute = 3*as.numeric(t50(c(tempt = tempt_mean, kappa = kappa_mean)))
    } else  {
      max_minute = 3*as.numeric(t50(c(tempt = tempt_mean, beta = beta_mean)))
    }
    # Round to hours
    max_minute = ((max_minute %/% 60) + 1) * 60
  }
  minute = c(seq(0,20, by = 5), seq(30, max_minute, by = 30))

  # Record
  rec = tibble(
    record = sprintf("rec_%02d", 1:n_records),
    v0 = round(rnorm(n_records, v0_mean, v0_std)),
    tempt = pmax(rnorm(n_records, tempt_mean, tempt_std), tempt_mean/3))
  if (model_name == "linexp") {
    rec$kappa = pmax(rnorm(n_records, kappa_mean, kappa_std), kappa_mean/3)
  } else {
    rec$beta = pmax(rnorm(n_records, beta_mean, beta_std), beta_mean/3)
  }

  # Noise term
  if (is.null(student_t_df) || student_t_df == 0 ) {
    if (noise == 0)
      warning("With noise == 0, non-linear fits might fail.")
    noise_d = rnorm(nrow(rec)*length(minute), 0, noise)
  } else {
    if (student_t_df < 2) {
      student_t_df = 2
      warning("Degrees of freedom of Student-t noise was set to 2")
    }
    noise_d = noise*rt(nrow(rec)*length(minute), student_t_df)
  }
  # Data
  data = rec %>%
    rowwise() %>%
    do(
      tibble(
        record = .$record,
        minute = minute,
        vol = model(minute, .[[2]], .[[3]], .[[4]])
        ))  %>%
    ungroup() %>%
    mutate(
      vol = vol + noise_d,
      vol = round(pmax(vol, abs(sample(noise_d))))
    )
  # Remove missing
  if (missing != 0) {
    missing1 = min(max(abs(missing), 0),0.5)
    if (missing1 != missing)
      warning("Fraction of missing was set to ", missing1)
    data = data[-sample.int(nrow(data), nrow(data)*missing1),]
  }

  # Add descriptor as comment
  kb = ifelse(model_name == "linexp",
              sprintf("kappa = %.2f\U00B1%.2f", kappa_mean,  kappa_std),
              sprintf("beta = %.2f\U00B1%.2f", beta_mean,  beta_std))
  sf = ifelse(is.null(student_t_df) || student_t_df == 0, "Gaussian",
              paste0("Student-t ", student_t_df," df"))
  comment = sprintf("%.0f records, %s, v0 = %.0f\U00B1%.0f, tempt =  %.0f\U00B1%.0f, %s,\n%s noise amplitude = %.0f, %.0f%% missing",
          n_records, model_name, v0_mean, v0_std, tempt_mean, tempt_std, kb,
          sf, noise, 100*missing)
  comment(data) = comment
  # Add list format for use with stan
  #
  data1 = data %>%
    mutate(
      record_i = as.integer(as.factor(data$record))
    )

  stan_data = list(
    prior_v0 = median(data1$vol[data1$minute < 10]),
    n = nrow(data1),
    n_record = max(data1$record_i),
    record = data1$record_i,
    minute = data1$minute,
    volume = data1$vol)

  list(record = rec, data = data, stan_data = stan_data)
}


if (FALSE) {
  n_records = 10
  n_points = 10
  v0_mean = 400
  v0_std = 50
  tempt_mean = 60
  tempt_std = 15
  kappa_mean = 0.7
  kappa_std = 0.3
  beta_mean = 0.7
  beta_std = 0.3
  noise = 10
  missing = 0.1
  student_t_df = 2
  model = linexp
  model_name = "linexp"
  library(dplyr)
}
