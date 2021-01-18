#' Spending on Food by Household Income
#'
#' Data frame on the proportion of food expenses per household income. 38 house rents were evaluated in a random sample from a large city in the United States.
#'@usage data("FoodExpenditure")
#'@source Taken from Griffiths et al. (1993, Table 15.4).
#'@details Originally, the \code{proportion} column did not exist, it was created by the bayesbr package.
#' @format A data frame containing 38 observations on 3 variables.
#'\describe{
#'  \item{food}{household expenditures for food.}
#'  \item{income}{household income.}
#'  \item{proportion}{proportion of household income spent on food.}
#'  \item{persons}{number of persons living in household.}
#'}
#'@references
#'\doi{10.18637/jss.v034.i02} Cribari-Neto, F., and Zeileis, A. (2010). Beta Regression in R.
#'\emph{Journal of Statistical Software}, \bold{34}(2), 1--24.
#'@references
#'\doi{10.1080/0266476042000214501} Ferrari, S.L.P., and Cribari-Neto, F. (2004).
#'Beta Regression for Modeling Rates and Proportions.
#'\emph{Journal of Applied Statistics}, \bold{31}(7), 799--815.
#'@references
#'doi{10.1002/jae.3950090208} Griffiths, W.E., Hill, R.C., and Judge, G.G. (1993).
#'\emph{Learning and Practicing Econometrics}
#'New York: John Wiley and Sons.
#' @examples
#'data("FoodExpenditure", package = "bayesbr")
#'
#'bbr <- bayesbr(proportion ~ income + persons, data = FoodExpenditure,
#'              iter=100)
#'residuals(bbr, type="quantile")
#'
#'\donttest{
#'pmse <-pmse(proportion ~ income + persons, test.set=0.4,
#'           data = FoodExpenditure, iter=100)$PMSE
#'}
"FoodExpenditure"

