#' Cells data set
#'
#' @name cells
#' @docType data
#' @keywords datasets
#' @description Data set taken from \insertCite{2012_Crawley}{bellreg} and posteriorly analyzed by \insertCite{2020_Lemonte}{bellreg}. The data includes the count of infected blood cells per square millimetre on microscope slides prepared from n = 511 randomly selected individuals.
#' @format A data frame with 511 rows and 5 variables:
#' \itemize{
#'   \item cells: count of infected blood cells per square millimetre on microscope slides
#'   \item smoker: smoking status of the subject (0: smoker; 1: non smoker)
#'   \item gender: subject's gender (1: male; 0: female).
#'   \item age: subject's age categorized into three levels: young (\eqn{ \le 20}), mid (21 to 59), and old (\eqn{\ge 60}).
#'   \item weight: body mass score categorized into three levels: normal, overweight, obese.
#' }
#' @references
#'
#' \insertRef{2012_Crawley}{bellreg}
#'
#' \insertRef{2020_Lemonte}{bellreg}
#'
NULL


#' Faults data set
#'
#' @name faults
#' @docType data
#' @keywords datasets
#' @description Data set taken from \insertCite{1982_Hinde}{bellreg} and posteriorly analyzed by \insertCite{2018_Castellares}{bellreg}. The data contains the number of faults in rolls of fabric of different lengths.
#' @format A data frame with 32 rows and 2 variables:
#' \itemize{
#'   \item nf: number of faults in rolls of fabric of different lengths.
#'   \item lroll: length of the roll.
#' }
#'
#' @references
#' \insertAllCited{}
#'
NULL


