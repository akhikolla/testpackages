#' Load Covid-19 Data
#'
#' This function pulls Covid-19 data up to a certain date, for a specified country (and state, if \code{country_name = "Brazil"}).
#' The output of this function is in the correct format to be used directly into the model adjustment function
#' \code{\link{pandemic_model}} included in this package.\cr
#' \cr
#' The user must online for this function to work.
#'
#' @param country_name Character string specifying the country of interest.
#' Check \code{country_list()} for the list of countries available in the database.
#' @param state_name Optional character string specifying the state of interest - only brazilians states currently
#' available in the database. \code{state_name} should be either \code{NULL} or a string of length 2.
#' Check \code{state_list()} for the state abbreviations that will be used and the corresponding state names.
#' @param last_date Optional date, character or factor argument specifying the last date in the data.
#' It should be in the YYYY-MM-DD or YYYY/MM/DD format.The default is the most recent date available in the database.
#'
#' @return An object of S3 class \code{pandemicData}. It is a list with 3 items:
#' \describe{
#'   \item{\code{data}:}{ data frame with the number of cumulative cases, new cases, cumulative deaths and new deaths associated
#'   with Covid-19 for each date, up to the \code{last_date} in the specified region.}
#'   \item{\code{name}:}{ character string with the country name (and state name, if available).}
#'   \item{\code{population}:}{numeric object that contains the population size of the given region.}
#' }
#'
#' @examples
#' \dontrun{
#' load_covid("Brazil","MG")
#' load_covid(country_name = "India", last_date = "2020-06-15")
#' load_covid("US")
#' load_covid(country_name = "italy")}
#'
#' @source \url{https://github.com/CSSEGISandData/COVID-19}\cr
#' \url{https://github.com/covid19br/covid19br.github.io}
#'
#' @references
#' CovidLP Team, 2020. CovidLP: Short and Long-term Prediction for COVID-19. Departamento de Estatistica. UFMG,
#' Brazil. URL: \url{http://est.ufmg.br/covidlp/home/en/}
#'
#' @seealso  \code{\link{country_list}}, \code{\link{state_list}},  \code{\link{pandemic_model}},
#' \code{\link{posterior_predict.pandemicEstimated}}, \code{\link{pandemic_stats}} and
#' \code{\link{plot.pandemicPredicted}}.
#'
#' @importFrom curl has_internet
#' @importFrom utils read.csv
#' @importFrom stats variable.names
#' @importFrom tidyr pivot_longer
#' @export

load_covid <- function(country_name, state_name = NULL, last_date){

  if(curl::has_internet() == F) stop("The user must be online to use this function")

  baseURLbr <-"https://raw.githubusercontent.com/covid19br/covid19br.github.io/master/dados"
  baseURL<-"https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series"
  country_list <- country_list()
  state_list<- state_list()[["state_abb"]]

  # Error checking
  if(nchar(state_name) == 2 && !is.null(state_name)) {state_name <- toupper(state_name)}
  if(class(country_name) == "character" ) country_name <- capitalize(tolower(country_name)) else stop("country_name must be a character string")
  if(!is.null(state_name) && class(state_name) != "character") stop("state_name must be a string of length 2")
  if(!(country_name %in% country_list)) stop("This country_name could not be found in the database. Use country_list() for available options")
  if(country_name != "Brazil" && !is.null(state_name)) warning("Selected country_name does not have state_name options available")
  if(nchar(state_name) != 2 && !is.null(state_name)) stop("state_name must be a string of length 2")
  if(length(country_name) != 1 || length(state_name) > 1) stop("country_name and state_name arguments cannot be vectors.")
  if(!(state_name %in% state_list) && !is.null(state_name)) stop("This state_name could not be found in the database. Use state_list() for available options")

  if (country_name != "Brazil") { # Loading Brazil country and states data. Different database as other countries
    covidworld <- try(utils::read.csv(file.path(baseURL,"time_series_covid19_confirmed_global.csv"), check.names=FALSE, stringsAsFactors=FALSE))
    if (class(covidworld) == "try-error") stop("Something went wrong retrieving the data from the repository. If the problem persists, please try again later or contact us at covidlp.team@gmail.com.")
    current_date <- as.Date(stats::variable.names(covidworld[dim(covidworld)[2]]), format = "%m/%d/%y")
  } else {
    covidbr <- try(utils::read.csv(file.path(baseURLbr,"EstadosCov19.csv"), check.names=FALSE, stringsAsFactors=FALSE))
    if (class(covidbr) == "try-error") stop("Something went wrong retrieving the data from the repository. If the problem persists, please try again later or contact us at covidlp.team@gmail.com.")
    current_date <- as.Date(max(covidbr$data))
  }
  if(missing(last_date)) last_date <- current_date
  if(last_date > current_date) warning(paste0("Invalid last_date. Database only contains data up to ", current_date))
  if(class(last_date) == "character" || class(last_date) == "factor") {
    last_date<- try(as.Date(last_date))
    if(class(last_date) == "try-error" || is.na(last_date))
      stop("last_date format must be YYYY-MM-DD or YYYY/MM/DD")
  }
  if(last_date < "2020-01-23") stop("last_date can't be earlier than 2020-01-23")

  if(country_name == "Brazil"){ # Data treatment

    if(is.null(state_name)){

      pop <- as.numeric(br_pop$pop[which(br_pop$uf == "BR")])

      Y <- utils::read.csv(file.path(baseURLbr,"BrasilCov19.csv"), check.names=FALSE, stringsAsFactors=FALSE)
      Y <- dplyr::rename(Y,date = "data",
               cases = "casos.acumulados",
               deaths = "obitos.acumulados",
               new_cases = "novos.casos",
               new_deaths = "obitos.novos")
      Y <- dplyr::mutate(Y,date = as.Date(Y$date))
      Y <- dplyr::select(Y,"date", "cases", "deaths", "new_cases", "new_deaths")
      Y <- dplyr::arrange(Y,Y$date)
      Y <- dplyr::filter(Y,Y$date >= '2020-01-23' & Y$date <= last_date)
    } else {
      pop <- as.numeric(br_pop$pop[which(br_pop$uf == state_name)])

      Y <- covidbr
      Y <- dplyr::rename(Y,name = "estado",
               date = "data",
               cases = "casos.acumulados",
               deaths = "obitos.acumulados",
               new_cases = "novos.casos",
               new_deaths = "obitos.novos")
      Y <- dplyr::mutate(Y,date = as.Date(Y$date))
      Y <- dplyr::arrange(Y,Y$name, Y$date)
      Y <- dplyr::filter(Y,Y$date >= '2020-01-23' & Y$date <= last_date & Y$name == state_name)
      Y <- dplyr::select(Y,"date", "cases", "deaths", "new_cases", "new_deaths")
    }
  } else{
    covid19_confirm <- covidworld
    covid19_confirm <- dplyr::select(covid19_confirm,-"Lat", -"Long")
    covid19_confirm <- tidyr::pivot_longer(covid19_confirm,-(1:2), names_to = "date", values_to = "confirmed")
    covid19_confirm <- dplyr::mutate(covid19_confirm,date = as.Date(covid19_confirm$date, format="%m/%d/%y"))
    covid19_confirm <- dplyr::rename(covid19_confirm,country = 'Country/Region', state = 'Province/State')

    covid19_deaths <- utils::read.csv(file.path(baseURL,"time_series_covid19_deaths_global.csv"), check.names=FALSE, stringsAsFactors=FALSE)
    covid19_deaths <- dplyr::select(covid19_deaths,-"Lat", -"Long")
    covid19_deaths <- tidyr::pivot_longer(covid19_deaths,-(1:2), names_to = "date", values_to ="deaths")
    covid19_deaths <- dplyr::mutate(covid19_deaths,date = as.Date(covid19_deaths$date, format="%m/%d/%y"))
    covid19_deaths <- dplyr::rename(covid19_deaths,country = 'Country/Region', state = 'Province/State')

    covid19 <- dplyr::left_join(covid19_confirm, covid19_deaths, by = c("state", "country", "date"))

    pop <- country_pop$pop[which(country_pop$country == country_name)]

    eval(parse(text="
    Y <- covid19
    Y <- dplyr::filter(Y,country == country_name)
    Y <- dplyr::mutate(Y,confirmed_new = confirmed - lag(confirmed, default=0),
             deaths_new = deaths - lag(deaths, default=0))
    Y <- dplyr::arrange(Y,date,state)
    Y <- dplyr::group_by(Y,date, country)
    Y <- dplyr::summarize(Y,cases = sum(confirmed, na.rm = T),
                deaths = sum(deaths, na.rm = T),
                new_cases = as.integer(sum(confirmed_new, na.rm = T)),
                new_deaths = as.integer(sum(deaths_new, na.rm = T)))
    Y <- dplyr::select(Y,date, cases, deaths, new_cases, new_deaths)
    Y <- dplyr::arrange(Y,date)
    Y <- dplyr::filter(Y,date >= '2020-01-23' & date <= last_date)
    "))
  }

  # Inconsistencies fail safe
  while(any(Y$new_cases < 0)){
    pos <- which(Y$new_cases < 0)
    for(j in pos){
      Y$new_cases[j-1] = Y$new_cases[j] + Y$new_cases[j-1]
      Y$new_cases[j] = 0
      Y$cases[j-1] = Y$cases[j]
    }
  }

  while(any(Y$new_deaths < 0)){
    pos <- which(Y$new_deaths < 0)
    for(j in pos){
      Y$new_deaths[j-1] = Y$new_deaths[j] + Y$new_deaths[j-1]
      Y$new_deaths[j] = 0
      Y$deaths[j-1] = Y$deaths[j]
    }
  }

  if(dim(Y)[1] == 0) warning("last_date assignment resulted in an empty data frame")

  list_out = list(data = as.data.frame(Y),
                  name = ifelse(is.null(state_name), paste0(country_name), paste0(country_name,"_",state_name)),
                  population = pop)
  class(list_out) = "pandemicData"
  return(list_out)

}

#' Capitalize country names
#'
#' This function capitalizes country_name inputs in order to match the capitalization in the database.
#'
#' @param x A character string
#'
#' @return Capitalized character string x to match the spelling and capitalization in the database.
#'
#' @noRd
#'

capitalize <- function(x) {
  s <- gsub("\\b(\\w)", "\\U\\1", x, perl = TRUE)
  s<- gsub("\\bAnd\\b", "and", s)
  if(length(s)== 1 && s == "Us") {s <-"US"}
  if(length(s)== 1 && s =="Cote D'ivoire") {s <-"Cote d'Ivoire"}
  s
}
