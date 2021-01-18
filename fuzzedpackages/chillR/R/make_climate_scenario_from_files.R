#' Make climate scenario from multiple saved csv files
#' 
#' Many climate scenarios we may want to plot consist of data stored across many files. These
#' files typically contain certain character strings that mark, e.g. the RCP scenario or the
#' point in time. This function facilitates accessing such files by allowing the specification
#' of search string (criteria_list), according to which files are selected. They are then
#' converted into climate_scenario files that can become part of a list passed to
#' plot_climate_scenarios for plotting.
#' 
#' @param metric_folder character string specifying the folder holding the files, from which the 
#' scenario is to be built.
#' @param criteria_list list of character vectors that specify parts of the file names that are
#' common to all files of a particular scenario. These can be single strings or vectors of string.
#' In the latter case, occurrence of either of the elements in a file name is sufficient. The
#' selection criteria are applied iteratively, i.e. first all files containing the first element of
#' 'criteria_list' are selected, then those containing the second element, and so forth.
#' @param caption vector of up to three character strings indicating the caption to be displayed
#' in the respective plot panel; the elements of this vector are displayed on different lines.
#' If caption_above==TRUE in plot_climate_scenario, only the first element is displayed.
#' @param time_series Boolean, indicating if the scenario contains a time series.
#' @param labels numeric vector containing labels for the time scenarios - only used for time series.
#' @param historic_data a data.frame containing at least two columns named the same as 'metric'
#' and 'year_name'.
#' 
#' 
#' @return a climate scenario object, which can be part of a list supplied to plot_climate_scenarios.
#' 
#' @author Eike Luedeling
#' @keywords utility
#' @examples
#' 
#' # historic_scenario<-make_climate_scenario(metric_folder=chillout_folder,
#' #                                          criteria_list=list(cult,c(1975,2000,2015)),
#' #                                          caption=c("Historic","data"),
#' #                                         time_series=TRUE,
#' #                                         labels=c(1975,2000,2015),
#' #                                         historic_data=historic_data)
#'   
#' @export make_climate_scenario_from_files
make_climate_scenario_from_files<-function(metric_folder,criteria_list,caption=NULL,time_series=FALSE,labels=NULL,historic_data=NULL)
{
  filelist<-list.files(metric_folder,full.names = TRUE)
  for(cr in criteria_list)
  {tokeep<-c()
  for(i in 1:length(cr)) tokeep<-c(tokeep,filelist[grep(cr[i],filelist)])
  filelist<-tokeep}
  
  dataset<-list()
  for (i in 1:length(filelist))
    dataset[[i]]<-read.csv(filelist[i])
  climate_scenario<-list(data=dataset,caption=caption,time_series=time_series,labels=labels,
                         historic_data=historic_data)
  return(climate_scenario)
}
