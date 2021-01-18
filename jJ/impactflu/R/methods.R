# Methods from Tokars (2018)
# Arseniy Khvorov
# Created 2020/01/03
# Last edit 2020/01/03

#' Analysis methods from Tokars (2018)
#'
#' Method 1 was said to be as current. Method 3 was determined to be the least
#' biased.
#'
#' @param init_pop_size Integer initial population size
#' @param vaccinations Integer vector counts of vaccinations
#' @param cases Integer vector counts of cases
#' @param ve Vector vaccine effectiveness. If length 1, assumed to not vary with
#'   time.
#'
#' @return
#' A \link[tibble]{tibble} with the following columns (method-dependent):
#'   \item{cases}{Observed cases}
#'   \item{vaccinations}{Observed vaccinations}
#'   \item{ve}{Assumed vaccine effectiveness}
#'   \item{pvac}{Proportion of the starting population vaccinated}
#'   \item{vc_lag}{Vaccine coverage lagged}
#'   \item{pops}{Susceptible population}
#'   \item{pflu}{Infection risk}
#'   \item{popn}{Non-cases is absence of vaccination}
#'   \item{cases_novac}{Cases in absence of vaccination}
#'   \item{avert}{Expected number of vaccinations}
#'
#' @references Tokars JI, Rolfes MA, Foppa IM, Reed C. An evaluation and update
#'   of methods for estimating the number of influenza cases averted by
#'   vaccination in the United States. Vaccine. 2018;36(48):7331–7337.
#'   doi:10.1016/j.vaccine.2018.10.026
#'
#' @importFrom tibble as_tibble
#' @importFrom dplyr mutate group_by ungroup summarise
#' @importFrom lubridate year month ymd
#'
#' @export
#'
#' @examples
#' library(dplyr)
#'
#' # Simulate a population
#' nsam <- 1e6L
#' ndays <- 304L
#' pop_tok <- sim_reference(
#'   init_pop_size = nsam,
#'   vaccinations = generate_counts(nsam, ndays, 0.55, mean = 100, sd = 50),
#'   cases_novac = generate_counts(nsam, ndays, 0.12, mean = 190, sd = 35),
#'   ve = 0.48,
#'   lag = 14,
#'   deterministic = TRUE
#' )
#'
#' # Summarise by month
#' pop_tok_month <- pop_tok %>%
#'   mutate(
#'     datestamp = generate_dates(
#'       timepoint, lubridate::ymd("2017-08-01"), "day"
#'     ),
#'     year = lubridate::year(datestamp),
#'     month = lubridate::month(datestamp)
#'  ) %>%
#'  group_by(year, month) %>%
#'  summarise(
#'    vaccinations = sum(vaccinations), cases = sum(cases), ve = mean(ve)
#'  ) %>%
#'  ungroup()
#'
#' # Estimate averted cases using the two different methods
#' m1 <- method1(
#'   nsam, pop_tok_month$vaccinations, pop_tok_month$cases, pop_tok_month$ve
#' )
#' m3 <- method3(
#'   nsam, pop_tok_month$vaccinations, pop_tok_month$cases, pop_tok_month$ve
#' )
#' sum(m1$avert)
#' sum(m3$avert)
method1 <- function(init_pop_size, vaccinations, cases, ve) {
  check_counts(vaccinations, cases, ve)
  as_tibble(method1_cpp(init_pop_size, vaccinations, cases, ve))
}

#' @rdname method1
#' @export
method3 <- function(init_pop_size, vaccinations, cases, ve) {
  check_counts(vaccinations, cases, ve)
  as_tibble(method3_cpp(init_pop_size, vaccinations, cases, ve))
}
