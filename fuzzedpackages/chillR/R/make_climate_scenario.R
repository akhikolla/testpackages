#' Make climate scenario
#' 
#' Function to make climate scenarios for plotting from a list of climate metric data, e.g.
#' produced by tempResponse_daily_list.
#' 
#' @param metric_summary character string specifying the folder holding the files, from which the 
#' scenario is to be built.
#' @param caption vector of up to three character strings indicating the caption to be displayed
#' in the respective plot panel; the elements of this vector are displayed on different lines.
#' If caption_above==TRUE in plot_climate_scenario, only the first element is displayed.
#' @param labels numeric vector containing labels for the scenarios. This defaults to the names
#' of elements in metric_summary.
#' @param time_series Boolean, indicating if the scenario contains a time series.
#' @param historic_data a data.frame containing a dataset of historic observations that is similar
#' in structure to metric_summary (should have column indicating the year and the metric to be
#' plotted, with identical names to metric_summary). Defaults to NULL, which means that no
#' historic data is included.
#' @param add_to list of climate scenarios that the newly created one is to be added to.
#' 
#' 
#' @return a list of climate scenario objects, which can be supplied to plot_climate_scenarios.
#' 
#' @author Eike Luedeling
#' @keywords utility
#' @examples
#' 
#' chill<-chilling(stack_hourly_temps(fix_weather(KA_weather[which(KA_weather$Year>1990),]),
#'    latitude=50.4))
#' multi_chills<-list('2001'=chill,'2005'=chill,'2009'=chill)
#' chills_to_plot<-make_climate_scenario(multi_chills,caption=c("Historic","data"),
#'    time_series=TRUE,historic_data=chill)
#' chills_to_plot<-make_climate_scenario(multi_chills,caption=c("Future1"),add_to=chills_to_plot)
#' chills_to_plot<-make_climate_scenario(multi_chills,caption=c("Future2"),add_to=chills_to_plot)
#' plot_climate_scenarios(chills_to_plot,metric="Chill_portions",metric_label="Chill Portions")
#' 
#'   
#' @export make_climate_scenario
make_climate_scenario<-function(metric_summary,
                                caption=NULL,labels=names(metric_summary),time_series=FALSE,
                                historic_data=NULL,add_to=NULL)
{
  
  if(!is.null(labels))
    {if(!length(labels)==length(metric_summary))
      stop("number of labels doesn't match number of elements",call. = FALSE)
    if(time_series) labels<-as.numeric(labels)
    }
  
  if(!is.null(historic_data))
    if(!is.data.frame(historic_data))
      if(is.data.frame(historic_data[[1]]))
        historic_data<-historic_data[[1]]
  
  if(is.null(add_to))
    output<-list(list(data=metric_summary,caption=caption,time_series=time_series,
            labels=labels,
            historic_data=historic_data)) else
            {
              if(!is.list(add_to)) stop("add_to needs to be a list of climate scenarios",call. = FALSE)
              add_to[[length(add_to)+1]]<-list(data=metric_summary,caption=caption,time_series=time_series,
                                               labels=labels,
                                               historic_data=historic_data) 
              output<-add_to
            }
  return(output)
}
  
