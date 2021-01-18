#' Travel Mode Choice
#'
#' The data set contains 210 observations on mode choice for travel between Sydney and Melbourne, Australia.
#'
#' @docType data
#' @usage data(TravelChoice)
#'
#' @format{A dataframe containing :
#' \describe{
#' \item{indv}{Id of the individual}
#' \item{mode}{available options: air, train, bus or car}
#' \item{choice}{a logical vector indicating as TRUE the transportation mode chosen by the traveler}
#' As category-specific variables:
#' \item{invt}{travel time in vehicle}
#' \item{gc}{generalized cost measure}
#' \item{ttme}{terminal waiting time for plane, train and bus; 0 for car}
#' \item{invc}{in vehicle cost}
#' As case-specific variables:
#' \item{hinc}{household income}
#' \item{psize}{traveling group size in mode chosen}
#' }
#' }
#'
#' @keywords datasets
#'
#' @references
#' Greene, W.H.  and  D.  Hensher (1997) \emph{Multinomial logit and discrete choice models} \emph{in}
#' Greene, W. H. (1997) \emph{LIMDEP version 7.0 user's manual revised}, Plainview, New York econometric software, Inc .
#' @source{
#' Download from on-line (18/09/2020) complements to Greene, W.H. (2011) Econometric Analysis, Prentice Hall, 7th Edition \url{http://people.stern.nyu.edu/wgreene/Text/Edition7/TableF18-2.csv}, Table F18-2.
#' }
#'
#' @examples
#' data(TravelChoice)
"TravelChoice"
