#' Simulated data: Human cases of campylobacteriosis and numbers of source samples positive for \emph{Campylobacter}.
#'
#' A simulated dataset containing the number of human cases of campylobacteriosis, the numbers of source samples
#' positive for \emph{Campylobacter} for each bacterial subtype, and the overall source prevalence.
#'
#' @format A list containing the human cases (`cases'), source samples (`sources'), prevalences (`prev') and true values (`truevals').
#'
#' \strong{cases:} data frame with 364 rows and 4 variables:
#' \describe{
#'   \item{Human}{number of human cases of campylobacteriosis}
#'   \item{Time}{Time id for the samples}
#'   \item{Location}{Location id for the samples}
#'   \item{Type}{MLST type id for the samples}
#' }
#'
#' \strong{sources:} data frame with 1092 rows and 4 variables
#' \describe{
#'   \item{Count}{number of source samples positive for campylobacteriosis}
#'   \item{Time}{Time id for the samples}
#'   \item{Source}{Source id for the samples}
#'   \item{Type}{MLST type id for the samples}
#' }
#'
#' \strong{prev:} data frame with 12 rows and 3 variables
#' \describe{
#'   \item{Value}{Prevalence value (number of positive samples divided by total number of samples)}
#'   \item{Time}{Time id for the samples}
#'   \item{Source}{Source id for the samples}
#' }
#'
#' \strong{truevals:} list containing a long format data frame for each model parameter giving the true value of that parameter.
#' \describe{
#'   \item{alpha}{A dataframe with 24 rows and 4 variables: Value contains the true alpha values,
#'   Time, Location and Source contain the time, location and source id's respectively.}
#'   \item{q}{A dataframe with 91 rows and 2 variables: Value contains the true q values, and
#'   Type contains the type id's.}
#'   \item{lambda_i}{A dataframe with 364 rows and 4 variables: Value contains the true lambda_i values,
#'   Time, Location and Type contain the time, location and type id's respectively.}
#'   \item{xi}{A dataframe with 24 rows and 4 variables: Value contains the true xi values,
#'   Time, Location and Source contain the time, location and source id's respectively.}
#'   \item{r}{A dataframe with 2184 rows and 5 variables: Value contains the true r values,
#'   Time, Type, Location and Source contain the time, type, location and source id's respectively.}
#' }
"sim_SA"
