#' Plot pandemic data
#'
#' S3 method that plots the predicted data into an interactive graphic.
#' @param x Output of function \code{\link{load_covid}}.
#' @param y This parameter does nothing
#' @param cases A string which indicates whether new cases, cumulative cases or both plots should be generated.
#' The argument must be a string being either 'new', 'cumulative' or 'both'.
#' @param color A string which indicates whether the plot should be colorful or in gray scales.
#' The argument must be a string being either TRUE or FALSE.
#' @param ... Currently unused.
#' @return A list containing two objects, new and cumulative "new" shows the plot for the new cases/deaths
#' and "cumulative" shows the plot for the cumulative cases/deaths. If any of them did not get plotted due to lack of prediction or due to the \code{cases}
#' argument, its value will return NULL.
#'
#' By default both plot are generated.
#' This is changed with the \code{cases} argument.
#' \item{\code{new}}{
#'   The plotted new cases/deaths. The plots are for daily new confirmed cases and daily new deaths.
#'   }
#'   \item{\code{cumulative}}{
#'   The plotted cumulative cases/deaths. The plots are for daily cumulative cases or daily cumulative deaths.
#'   }
#' @seealso \code{\link{load_covid}}.
#' @references
#' CovidLP Team, 2020. CovidLP: Short and Long-term Prediction for COVID-19. Departamento de Estatistica. UFMG,
#' Brazil. URL: \url{http://est.ufmg.br/covidlp/home/en/}
#' @examples
#' \dontrun{
#' dataMG = load_covid(country_name="Brazil",state_name = "MG",last_date="2020-10-01")
#' plot(dataMG)
#' dataJapan = load_covid(country_name="Japan",last_date="2020-10-01")
#' plot(dataJapan)
#' }
#' @importFrom plotly plot_ly add_trace layout
#'
#' @method plot pandemicData
#' @export
plot.pandemicData <- function(x,y,cases = "new",color = TRUE,...){

  if(class(x) == "pandemicData" | class(x) == "list") {
    test <- 'ok'
  } else {
    stop("Please use the output of the load_covid() function or a list")
  }

  cases = tolower(cases)
  if(!(cases %in% c("both","new","cumulative"))) stop("Invalid \'cases\' argument. Please read \'help(plot.pandemicData)\' for available options.")

  cat(paste0("Plotting Data \n"))

  if(cases == "both") {
    terms <- c("new","cumulative")
  } else
    terms = cases

  plots = list()
  dados = x$data
  blu <- 'rgb(100, 140, 240)'
  dblu <- 'rgb(0, 0, 102)'
  red <- 'rgb(200, 30, 30)'
  dred <- 'rgb(100, 30, 30)'
  dt_format <- "%d/%b/%y"

  title <- ifelse(grepl("_", x$name, fixed=TRUE) == TRUE, paste(substr(x$name,1,5), "/", substr(x$name,8,9),sep = ""),
                  x$name)

  for(selTerm in terms){

    if(selTerm == "new"){
      data_plot_cases = dados$new_cases
      data_plot_deaths = dados$new_deaths
    } else {
      data_plot_cases = dados$cases
      data_plot_deaths = dados$deaths
    }

    if(selTerm == "new"){

    } else {

    }

    yaxis <- ifelse(selTerm == "new","New Cases per Day", "Cumulative Cases")

    fig2 <- plotly::add_trace(plotly::plot_ly(dados),x = dados$date, y = data_plot_cases,
                              type = 'scatter', mode = 'lines+markers', name = "Confirmed cases",
                              hoverinfo = "x+y",
                              marker = list(
                                color =  ifelse(color == TRUE,blu,'rgb(160,160,160)'),
                                line = list(color = ifelse(color == TRUE,dblu,'rgb(160,160,160)'), width = 1),
                                size = 5
                              ),
                              line = list(
                                color = ifelse(color == TRUE,blu,'rgb(160,160,160)'),
                                width = 1.5
                              )
    )

    fig2 <- plotly::add_trace(fig2,x = dados$date, y = data_plot_deaths,
                              type = 'scatter', mode = 'lines+markers', name = "Deaths",
                              hoverinfo = "x+y",
                              marker = list(
                                color = ifelse(color == TRUE,red,'rgb(96,96,96)'),
                                line = list(color = ifelse(color == TRUE,dred,'rgb(96,96,96)'), width = 1),
                                size = 5
                              ),
                              line = list(
                                color = ifelse(color == TRUE,red,'rgb(96,96,96)'),
                                width = 1.5
                              )
    )

    fig2 <- plotly::layout(fig2,title = title,
                           ## Add pred dates to the x axis
                           xaxis = list(title = "",
                                        tickangle = -90,
                                        ## Add pred dates to the x axis
                                        dtick = 14*86400000,
                                        tickformat = dt_format
                           ),
                           yaxis = list(
                             title = yaxis,
                             hoverformat = '.0f', hoverinfo = "x+y"
                           ),           # Show only date/value on hover
                           legend = list(
                             x = 0.03,
                             y = 0.97,                             # Legend position
                             bgcolor = 'rgba(240, 240, 240, 0.5)'
                           ),
                           font = list(family = "Arial", size = 10)
                           )

    plots[[selTerm]] = fig2
  }

  if(cases == "both") terms = "newcumulative" else terms = cases
  cat("The generated plot(s) can be stored in a variable.\n")
  if (grepl("new",terms)) cat("New Cases plot: variable$new\n")
  if (grepl("cumulative",terms)) cat("Cumulative Cases plot: variable$cumulative\n")

  class(plots) = "plottedPandemicData"
  return(plots)
}
