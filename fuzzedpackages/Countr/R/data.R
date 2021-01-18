#' Fertility data
#'
#' Fertility data analysed by Winkelmann(1995). The data comes from the second
#' (1985) wave of German Socio-Economic Panel. The sample is formed by
#' 1,243 women aged 44 or older in 1985. The response variable is the number of
#' children per woman and explanatory variables are described in more details below.
#'
#' @format A data frame with 9 variables (5 factors, 4 integers) and 1243 observations:
#' \describe{
#' \item{\code{children}}{integer; response variable: number of children
#' per woman (integer).}
#' \item{\code{german}}{factor; is the mother German? (yes or no).}
#' \item{\code{years_school}}{integer; education measured as years of schooling.}
#' \item{\code{voc_train}}{factor; vocational training ? (yes or no)}
#' \item{\code{university}}{factor; university education ? (yes or no)}
#' \item{\code{religion}}{factor; mother's religion: Catholic, Protestant, Muslim
#' or Others (reference).}
#' \item{\code{rural}}{factor; rural (yes or no ?)}
#' \item{\code{year_birth}}{integer; year of birth (last 2 digits)}
#' \item{\code{age_marriage}}{integer; age at marriage}
#' }
#'
#' For further details, see Winlemann(1995).
#' @references
#' \insertRef{winkelmann1995duration}{Countr}
#'
"fertility"

#' Football data
#'
#' Final scores of all matches in the English Premier League from seasons
#' 2009/2010 to 2016/2017.
#'
#' The data were collected from \url{http://www.football-data.co.uk/englandm.php}
#' and slightly formatted and simplified.
#'
#' @format a data.frame with 6 columns and 1104 observations:
#' \describe{
#' \item{\code{seasonId}}{integer season identifier (year of the first month of competition).}
#' \item{\code{gameDate}}{POSIXct game date and time.}
#' \item{\code{homeTeam,awayTeam}}{character home and away team name}.
#' \item{\code{homeTeamGoals,awayTeamGoals}}{integer number of goals scored by the home
#' and the away team.}
#' }
#'
"football"
