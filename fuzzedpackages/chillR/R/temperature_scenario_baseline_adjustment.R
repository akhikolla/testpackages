#' Make temperature scenario relative to a particular baseline
#' 
#' When interpreting future (or past) temperature scenarios that provide absolute temperatures,
#' it is important to consider the temperature baseline, i.e. a temperature scenario produced with
#' similar models and methods that corresponds to the current temperature regime. Such baselines are
#' normally available from the same source that provided the future scenarios. This function implements
#' this adjustment.
#' The function can be used for two situations:
#' 1) two absolute temperature scenarios: the output is the difference between the scenarios, i.e.
#' a relative temperature scenario describing the difference between monthly temperature extreme
#' means between the two scenarios.
#' 2) two relative temperature scenarios: the output is a relative temperature scenario that describes
#' the difference between the scenario year of the temperature_scenario and the baseline year of the
#' baseline_temperature_scenario. This only works if the scenario_year of the baseline_temperature_scenario
#' is the same as the reference_year of the temperature_scenario.
#' 
#' @param baseline_temperature_scenario baseline temperature scenario (e.g. produced with
#' 'extract_temperatures_from_grids'). This is a temperature scenario object, consisting of the
#' following elements: 'data' = data.frame with two columns Tmin and Tmax containing absolute
#' (normally monthly) mean minimum and maximum temperatures; 'reference_year' = the year the scenario
#' refers to (this is normally NA for absolute temperature scenarios, because they don't require
#' considering a reference scenario); 'scenario_type' = the scenario type, normally "absolute" (but can
#' also be "relative" or NA - then the type is automatically assigned); 'labels' = elements
#' attached to the input temperature_scenario. A subset of these elements can
#' also be specified, but 'data' must be present.
#' 
#' corresponds to the period, for
#' which observed weather records are available. The first step in this is to compute the change
#' generation. This is needed for making extracting information on prospective changes from sense of future climate scenarios, climate change analyses
#' @param temperature_scenario can be one of three options:
#' 1) a data.frame with two columns Tmin and Tmax and n_intervals (default: 12) rows containing
#' temperature changes for all time intervals, or absolute temperatures for these intervals.
#' 2) a temperature scenario object, consisting of the following elements: 'data' = a data frame with
#' n_intervals elements containing the absolute or relative temperature information (as in input option 1);
#' 'scenario_year' = the year the scenario is representative of; 'reference_year' = the year the
#' scenario is representative of; 'scenario_type' = the scenario type ('absolute' or 'relative' - if
#' NA, this is assigned automatically); 'labels' = and elements attached to the input
#' temperature_scenario as an element names 'labels'. A subset of these elements can also be specified,
#' but 'data' must be present.
#' 3) a list of elements of type 1 or 2. Then the adjustment is done for all elements.
#' @param  temperature_check_args list of arguments to be passed to the check_temperature_scenario function.
#' Check documentation of that function for details.
#' @param warn_me boolean variable specifying whether warnings should be shown. Defaults to TRUE.
#' @param required_variables character vectors containing names of variables that must be included in the
#' scenario.
#' 
#' @return temperature scenario object, consisting of the following elements: 'data' = a data frame with
#' n_intervals elements containing the absolute or relative temperature information. 'reference_year' =
#' the year the scenario is representative of. 'scenario_type' = the scenario type ('absolute' or 'relative');
#' 'labels' = and elements attached to the input temperature_scenario as an element names 'labels'.
#' 
#' The function also returns warnings, where elements are missing or the scenario_type appears to be
#' wrong, and it stops with an error, if the scenario isn't specified in a format that is usable by
#' chillR.
#' 
#' @author Eike Luedeling
#' @keywords utility
#' @examples
#' 
#' baseline_temperature_scenario<-list(data=data.frame(Tmin=c(1,1,1,1,1,1,1,1,1,1,1,1),
#'                                                     Tmax=c(1,1,1,1,1,1,1,1,1,1,1,1)),
#'                                                     scenario_year=1990,
#'                                                     reference_year=1975,
#'                                                     scenario_type="relative")
#'                                                     
#' temperature_scenario<-list(data=data.frame(Tmin=c(4,4,4,4,4,4,4,4,4,4,4,4),
#'                                            Tmax=c(4,4,4,4,4,4,4,4,4,4,4,4)),
#'                                            scenario_year=2000,
#'                                            reference_year=1990,
#'                                            scenario_type="relative")
#'                                            
#' relative_temperature_scenario<-temperature_scenario_baseline_adjustment(
#'                      baseline_temperature_scenario,temperature_scenario,
#'                      temperature_check_args=NULL)
#'                                                                        
#' baseline_temperature_scenario<-list(data=data.frame(Tmin=c(-5,-2,2,5,10,12,15,15,12,10,5,1),
#'                                                     Tmax=c( 1, 4,7,10,15,18,22,24,17,15,11,6)),
#'                                                     scenario_year=1980,
#'                                                     reference_year=NA,
#'                                                     scenario_type="absolute")
#'                                                     
#' temperature_scenario<-list(data=data.frame(Tmin=c(-3,0,4,7,12,14,17,17,14,12,7,3),
#'                                                     Tmax=c(3,6,9,12,17,20,24,26,19,17,13,8)),
#'                                            scenario_year=2000,
#'                                            reference_year=NA,
#'                                            scenario_type="absolute")
#'                                            
#' relative_temperature_scenario<-temperature_scenario_baseline_adjustment(
#'                      baseline_temperature_scenario,temperature_scenario,
#'                      temperature_check_args=NULL)
#'                                             
#'  
#' @export temperature_scenario_baseline_adjustment
temperature_scenario_baseline_adjustment<-function(baseline_temperature_scenario,temperature_scenario,
                                                   temperature_check_args=NULL,warn_me=TRUE,required_variables=c("Tmin","Tmax"))
{
  if(!is.data.frame(baseline_temperature_scenario[[1]]))
    baseline_temperature_scenario<-baseline_temperature_scenario[[1]]
  
  if(is.data.frame(temperature_scenario[[1]]))
    temperature_scenario<-list(temperature_scenario)
  
  if(required_variables[1]=="monthly_min_max_temps")
    required_variables<-c(paste("tasmin",1:12,sep=""),
                          paste("tasmax",1:12,sep=""))
  
  temperature_scenario_check_n_intervals<-12
  temperature_scenario_check_check_scenario_type<-TRUE
  temperature_scenario_check_scenario_check_thresholds<-c(-5,10)
  temperature_scenario_check_update_scenario_type<-TRUE
  temperature_scenario_check_warn_me<-TRUE
  
  if (!is.null(temperature_check_args))
  {if("n_intervals" %in% names(temperature_check_args))
    temperature_scenario_check_n_intervals<-temperature_check_args$n_intervals
  if("check_scenario_type" %in% names(temperature_check_args))
    temperature_scenario_check_check_scenario_type<-temperature_check_args$check_scenario_type
  if("scenario_check_thresholds" %in% names(temperature_check_args))
    temperature_scenario_check_scenario_check_thresholds<-temperature_check_args$scenario_check_thresholds
  if("update_scenario_type" %in% names(temperature_check_args))
    temperature_scenario_check_update_scenario_type<-temperature_check_args$update_scenario_type
  if("warn_me" %in% names(temperature_check_args))
    temperature_scenario_check_warn_me<-temperature_check_args$warn_me
  }
  
  if(is.null(temperature_scenario)) stop("No temperature scenario provided",call. = FALSE)   
  
  baseline_temperature_scenario<-check_temperature_scenario(baseline_temperature_scenario,
                                                            n_intervals=temperature_scenario_check_n_intervals,
                                                            check_scenario_type=temperature_scenario_check_check_scenario_type,
                                                            scenario_check_thresholds=temperature_scenario_check_scenario_check_thresholds,
                                                            update_scenario_type=temperature_scenario_check_update_scenario_type,
                                                            warn_me=temperature_scenario_check_warn_me,
                                                            required_variables=required_variables)
  for(i in 1:length(temperature_scenario))
    temperature_scenario[[i]]<-check_temperature_scenario(temperature_scenario[[i]],
                                                          n_intervals=temperature_scenario_check_n_intervals,
                                                          check_scenario_type=temperature_scenario_check_check_scenario_type,
                                                          scenario_check_thresholds=temperature_scenario_check_scenario_check_thresholds,
                                                          update_scenario_type=temperature_scenario_check_update_scenario_type,
                                                          warn_me=temperature_scenario_check_warn_me,
                                                          required_variables=required_variables)  
  
  out_scenario<-list()
  for(i in 1:length(temperature_scenario))
  {
    
    #subtracting two absolute temperature scenarios
    if(baseline_temperature_scenario$scenario_type=="absolute"&temperature_scenario[[i]]$scenario_type=="absolute")
    {out_scenario[[i]]<-temperature_scenario[[i]]
    if("GCM" %in% colnames(out_scenario[[i]]$data))
    {out_scenario[[i]]$data<-temperature_scenario[[i]]$data
    out_scenario[[i]]$data[,which(!colnames(out_scenario[[i]]$data)=="GCM")]<-temperature_scenario[[i]]$data[,which(!colnames(temperature_scenario[[i]]$data)=="GCM")]-
      baseline_temperature_scenario$data[,which(!colnames(temperature_scenario[[i]]$data)=="GCM")]} else
        out_scenario[[i]]$data<-temperature_scenario[[i]]$data-baseline_temperature_scenario$data
      
      out_scenario[[i]]$scenario_year<-temperature_scenario[[i]]$scenario_year
      out_scenario[[i]]$reference_year<-baseline_temperature_scenario$scenario_year
      out_scenario[[i]]$scenario_type<-"relative"
    }
    
    #adjusting the baseline of relative temperature scenarios (both are relative)
    #this adjusts the reference year
    if(baseline_temperature_scenario$scenario_type=="relative"&temperature_scenario[[i]]$scenario_type=="relative")
    {if(is.na(baseline_temperature_scenario$scenario_year)|is.na(temperature_scenario[[i]]$reference_year))
    {if (warn_me) warning("scenario year of baseline scenario and/or reference year of the temperature scenario",
                          " not specified - can't verify whether this is a valid transaction!",call. = FALSE)} else
                            if(!baseline_temperature_scenario$scenario_year==temperature_scenario[[i]]$reference_year)
                              stop("scenario year of the baseline scenario isn't equal to the reference year of the temperature scenario - ",
                                   "these scenarios aren't compatible.",call. = FALSE)
      out_scenario[[i]]<-temperature_scenario[[i]]
      if("GCM" %in% colnames(out_scenario[[i]]$data))
      {out_scenario[[i]]$data<-temperature_scenario[[i]]$data
      out_scenario[[i]]$data[,which(!colnames(out_scenario[[i]]$data)=="GCM")]<-temperature_scenario[[i]]$data[,which(!colnames(temperature_scenario[[i]]$data)=="GCM")]+
        baseline_temperature_scenario$data[,which(!colnames(temperature_scenario[[i]]$data)=="GCM")]} else
          out_scenario[[i]]$data<-temperature_scenario[[i]]$data+baseline_temperature_scenario$data
      out_scenario[[i]]$scenario_year<-temperature_scenario[[i]]$scenario_year
      out_scenario[[i]]$reference_year<-baseline_temperature_scenario$reference_year
      out_scenario[[i]]$scenario_type<-"relative"
    }
  }
  names(out_scenario)<-names(temperature_scenario)
  
  return(out_scenario)
}