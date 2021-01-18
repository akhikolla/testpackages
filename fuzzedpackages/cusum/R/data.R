#' Non-Risk-adjusted Performance Data
#'
#' Generated performance data of indicator 17/1 54030: Preoperative stay over 24 hours for patients with proximal femur fracture.
#'
#' Patient outcomes were simulated based on average national failure rate. Two years are provided, so Phase I and Phase II can be defined.
#'
#' @format A data frame with 2000 rows and 3 variables:
#' \describe{
#'   \item{t}{Sequence of observations}
#'   \item{y}{Patient outcome}
#'   \item{year}{Year of treatment}
#' }
#' @source Data for simulation was provided by Bavarian Agency of Quality Assurance (BAQ), Munich Germany.
#'
#' Description of performance indicator (in German): \url{https://iqtig.org/downloads/auswertung/2016/17n1hftfrak/QSKH_17n1-HUEFTFRAK_2016_QIDB_V02_2017-04-26.pdf}
"cusum_example_data"

#' Risk-adjusted Performance Data
#'
#' Generated performance data of indicator: Ratio of observed to expected cases of severe stroke or death under open carotid stenosis surgery.
#'
#' Individual patient risk scores were drawn from actual hospital data and patient outcomes were simulated. Two years are provided, so Phase I and Phase II can be defined.

#' @format A data frame with 2000 rows and 4 variables:
#' \describe{
#'   \item{t}{Sequence of observations}
#'   \item{y}{Patient outcome}
#'   \item{score}{Patient risk score}
#'   \item{year}{Year of treatment}
#' }
#'
#' @source Data for simulation was provided by Bavarian Agency of Quality Assurance (BAQ), Munich Germany.
#'
#' Description of performance indicator (in German): \url{https://iqtig.org/downloads/auswertung/2016/10n2karot/QSKH_10n2-KAROT_2016_QIDB_V02_2017-04-26.pdf}
"racusum_example_data"

#' Group-sequential Non-Risk-adjusted Performance Data with Block Identifier
#'
#' Generated performance data of indicator 17/1 54030: Preoperative stay over 24 hours for patients with proximal femur fracture. 
#' 
#'
#' Patient outcomes were simulated based on average national failure rate. Two years are provided, so Phase I and Phase II can be defined.
#'
#' @format A data frame with 2000 rows and 4 variables:
#' \describe{
#'   \item{t}{Sequence of observations}
#'   \item{y}{Patient outcome}
#'   \item{year}{Year of treatment}
#'   \item{block_identifier}{Continuous block identifier}
#' }
#' @source Data for simulation was provided by Bavarian Agency of Quality Assurance (BAQ), Munich Germany.
#'
#' Description of performance indicator (in German): \url{https://iqtig.org/downloads/auswertung/2016/17n1hftfrak/QSKH_17n1-HUEFTFRAK_2016_QIDB_V02_2017-04-26.pdf}
"gscusum_example_data"


#' Group-sequential Risk-adjusted Performance Data with Block Identifier
#'
#' Generated performance data of indicator: Ratio of observed to expected cases of severe stroke or death under open carotid stenosis surgery.
#'
#' Individual patient risk scores were drawn from actual hospital data and patient outcomes were simulated. Two years are provided, so Phase I and Phase II can be defined.

#' @format A data frame with 2000 rows and 4 variables:
#' \describe{
#'   \item{t}{Sequence of observations}
#'   \item{y}{Patient outcome}
#'   \item{score}{Patient risk score}
#'   \item{year}{Year of treatment}
#'   \item{block_identifier}{Continuous block identifier}
#' }
#'
#' @source Data for simulation was provided by Bavarian Agency of Quality Assurance (BAQ), Munich Germany.
#'
#' Description of performance indicator (in German): \url{https://iqtig.org/downloads/auswertung/2016/10n2karot/QSKH_10n2-KAROT_2016_QIDB_V02_2017-04-26.pdf}
"ragscusum_example_data"