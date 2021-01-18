#' Plot pandemic predictions
#'
#' S3 method that plots the predicted data into an interactive graphic.
#' @param x Output of function \code{\link{posterior_predict.pandemicEstimated}}.
#' @param y This parameter does nothing
#' @param term A string which indicates whether long term, short term or both plots should be generated.
#' The argument must be a string being either 'long', 'short' or 'both'.
#' @param color A logical variable indicating whether the plot should be colorful or in gray scales.
#' If TRUE then the plot is colorful otherwise the plot will be gray.
#' @param summary A logical variable indicating whether the plot should contains the summary statistics about the pandemic or not.
#' If TRUE then the plot shows the statistics otherwise no.
#' @param ... Currently unused.
#' @return A list containing two objects, short and long. "long" shows the long-term predictions
#' and "short" the short-term predictions. If any of them did not get plotted due to lack of prediction or due to the \code{term}
#' argument, its value will return NULL.
#'
#' By default only the long-term plot is generated, which means that the short term plot returns NULL.
#' This is changed with the \code{term} argument.
#' \item{\code{long}}{
#'   The plotted long-term prediction. The predictions are made on daily new confirmed cases or daily new deaths.
#'   }
#'   \item{\code{short}}{
#'   The plotted short-term prediction. The predictions are made on daily cumulative cases or daily cumulative deaths.
#'   }
#' @seealso \code{\link{posterior_predict.pandemicEstimated}}, \code{\link{pandemic_stats}} and \code{\link{plottedPandemic-objects}}.
#' @references
#' CovidLP Team, 2020. CovidLP: Short and Long-term Prediction for COVID-19. Departamento de Estatistica. UFMG,
#' Brazil. URL: \url{http://est.ufmg.br/covidlp/home/en/}
#' @examples
#' \dontrun{
#' dataMG = load_covid("Brazil","MG")
#' estimMG = pandemic_model(dataMG)
#' predMG = posterior_predict(estimMG)
#' statsMG = pandemic_stats(predMG)
#' plot(predMG)}
#' @importFrom plotly plot_ly add_trace layout
#'
#' @method plot pandemicPredicted
#' @export
plot.pandemicPredicted <- function(x,y,term = "long",color = TRUE,summary = TRUE,...){
  if(class(x) != "pandemicPredicted") stop("Please use the output of the posterior_predict.pandemicEstimated method()")
  term = tolower(term)
  if(!(term %in% c("long","short","both"))) stop("Invalid \'term\' argument. Please read \'help(plot.pandemicPredicted)\' for available options.")

  cat(paste0("Plotting ",x$cases_type," cases\n"))

  if(term == "both") terms = "longshort" else terms = term

  plots = list()
  for (selTerm in c("short","long")){
    if (!grepl(selTerm,terms)) next
    test_long = selTerm == "long"
    test_short = selTerm == "short"
    if(ncol(x$predictive_Long) < 1 & test_long){
      message("There are no long term predictions.")
      plots$long = NULL
    }
    if(ncol(x$predictive_Short) < 1 & test_short){
      message("There are no short term predictions.")
      plots$short = NULL
      return(plots)
    }

    outputs <- pandemic_stats(x)

    test_conf_death = x$cases_type == "confirmed"
    dados <- outputs$data$data
    if (test_long) data = outputs$LT_predict else data = outputs$ST_predict
    pred_summary <- outputs$LT_summary
    metric_leg <- ifelse(test_conf_death,"cases","deaths")
    last_date_n <- min(data$date) - 1
    pred_time <- length(data$date)
    mu_plot <- outputs$mu
    if (test_long) past = dados[[paste0("new_",metric_leg)]] else past = dados[[metric_leg]]

    ## some graphs colors and fonts
    blu <- 'rgb(100, 140, 240)'
    dblu <- 'rgb(0, 0, 102)'
    red <- 'rgb(200, 30, 30)'
    dred <- 'rgb(100, 30, 30)'

    color_lines <- ifelse(color == TRUE,ifelse(test_conf_death,dblu,dred),'rgb(192,192,192)')
    color_markers <- ifelse(color == TRUE,ifelse(test_conf_death,blu,red),'rgb(160,160,160)')
    color_mu <- ifelse(color==TRUE,'rgb(230, 115, 0)','rgb(96,96,96)')
    dt_format <- "%d/%b/%y"
    f1 <- list(family = "Arial", size = 10, color = "rgb(30, 30, 30)")

    # y_t median for future times
    fig2 <- plotly::add_trace(plotly::plot_ly(data),x = ~date, y = ~med, type= "scatter", name = "Prediction", mode = "markers",
                              marker=list(color='rgb(0,0,0)', dash='solid', width=2.5))

    # past observations
    fig2 <- plotly::add_trace(fig2,x = dados$date, y = past, type = 'scatter', mode = 'lines+markers', name = "Observed Data",
                              hoverinfo = "x+y",
                              marker = list(
                                color = color_markers,
                                line = list(color = 'rgb(0, 0, 0)', width = 1)),
                              line = list(
                                color = color_lines,
                                width = 1.5)
    )

    # Textual summary statistics and layout
    if (test_long) {
      fig2 <- plotly::layout(fig2,title = paste0("Prediction of ",ifelse(test_conf_death,"New Cases","Deaths")," - ",x$location),
                             xaxis = list(title = "",
                                          tickangle = -90,
                                          ## Add pred dates to the x axis
                                          dtick = 14*86400000,
                                          tickformat = dt_format
                             ),
                             yaxis = list(title = paste0(ifelse(test_conf_death,"New Cases","New Deaths")," per Day")),
                             shapes = list(type = "line", opacity = 0.7, line = list(color = "black", width = 1),
                                           y0 = 0, y1 = 1, yref = "paper",
                                           x0 = last_date_n, x1 = last_date_n),
                             legend = list(x = 0.03,
                                           y = 0.97,                             # Legend position
                                           bgcolor = 'rgba(240, 240, 240, 0.5)',font = f1),
                             font = list(family = "Arial", size = 10)

      )

      if(summary == TRUE){
        fig2 <- plotly::layout(fig2,annotations = list(text = paste0(
                                "<b>Peak:</b> ", format(pred_summary$peak_date_med,"%d/%b/%y"),
                                "<br><b>Peak 95% CI:</b> ",
                                paste0("(", format(pred_summary$peak_date_LB,"%d/%b/%y"),", ",
                                       format(pred_summary$peak_date_UB,"%d/%b/%y"), ")"),

                                "<br><b>Total number of ", metric_leg, ":</b> ",
                                round(pred_summary$total_cases_med,0), "<br>",
                                "<b>95% CI:</b> ",
                                paste0("(", round(pred_summary$total_cases_LB, 0),
                                       ", ", round(pred_summary$total_cases_UB, 0), ")"),
                                "<br><b>End (99%) of ",metric_leg,":</b> ",
                                format(pred_summary$end_date_med,"%d/%b/%y"),
                                "<br><b>95% CI:</b> ",
                                paste0("(", format(pred_summary$end_date_LB,"%d/%b/%y"),", ",
                                       format(pred_summary$end_date_UB,"%d/%b/%y"),")" ),
                                "</span>"
                              ),x = 0.97, y = 0.97, xref = "paper", yref = "paper",
                              font = list(family = "Arial", size = 10), align = "right",
                              showarrow = FALSE))
      }
      } else
        fig2 <- plotly::layout(fig2,title = paste("Prediction of",ifelse(test_conf_death,"New Cases","Deaths"),"- ",x$location),
                               xaxis = list(title = "",
                                            tickangle = -90,
                                            ## Add pred dates to the x axis
                                            dtick = 14*86400000,
                                            tickformat = dt_format),
                               yaxis = list(title = paste0("Cumulative ",ifelse(test_conf_death,"New Cases","New Deaths"))),
                               shapes = list(type = "line", opacity = 0.7, line = list(color = "black", width = 1),
                                             y0 = 0, y1 = 1, yref = "paper",
                                             x0 = last_date_n, x1 = last_date_n),
                               legend = list(x = 0.03, y = 0.97, bgcolor = 'rgba(240, 240, 240, 0.5)',
                                             font = f1),
                               font = list(family = "Arial", size = 10)
        )

    # Predicted 2.5% quantile curve
    fig2 <- plotly::add_trace(fig2,
                              x = c(dados$date[which(dados$date == last_date_n)],
                                    data$date[1:pred_time]),
                              y = c(past[which(dados$date == last_date_n)],
                                    data$q2.5[1:pred_time]),
                              showlegend = F,
                              name = "95% CI",
                              type = 'scatter',
                              mode = 'lines', hoverinfo = "x+y",
                              #fill='tonexty',
                              #dash='solid',
                              line = list(color = 'rgba(0, 0, 0, 1)',width=0)
    )

    # Predicted 97.5% quantile curve
    fig2 <- plotly::add_trace(fig2,
                              x = c(dados$date[which(dados$date == last_date_n)],
                                    data$date[1:pred_time]),
                              y = c(past[which(dados$date == last_date_n)],
                                    data$q97.5[1:pred_time]),
                              type = 'scatter',
                              mode = 'lines', hoverinfo = "x+y",
                              fill = 'tonexty',
                              name = "95% CI",
                              # fillcolor = 'rgba(150, 150, 150, 0.5)',
                              fillcolor = 'rgba(100, 100, 100, 0.5)',
                              line = list(color = 'rgba(0, 0, 0, 1)', width = 0)
    )

    # mu(t) median for past times
    if (test_long)
      fig2 <- plotly::add_trace(fig2,
                                x = mu_plot$date, y = mu_plot$mu,
                                type = 'scatter', mode = 'lines', hoverinfo = "none",
                                name = ("Estimated Mean"),
                                line = list(color = color_mu, dash = 'solid',
                                            width=1.5) ## width = 2.5)
      )

    cat(ifelse(test_long,"Long","Short"),"term plot created successfully\n")
    plots[[selTerm]] <- fig2
  }

  cat("The generated plot(s) can be stored in a variable.\n")
  if (grepl("long",terms)) cat("Long term plot: variable$long\n")
  if (grepl("short",terms)) cat("Short term plot: variable$short\n")

  class(plots) = "plottedPandemic"
  return(plots)
}
