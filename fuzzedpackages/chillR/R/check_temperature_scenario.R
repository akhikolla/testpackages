#' Check temperature scenario for consistency
#' 
#' chillR's temperature generation procedures require absolute or relative temperature scenarios.
#' This function checks these scenarios for consistency, regarding the data format, the reference
#' year, and whether they are relative or absolute scenarios (based on specified criteria).
#' 
#' Besides being able to validate classic temperature scenarios consisting of "Tmin" and "Tmax"
#' data, the function can also validate other datasets (e.g. outputs of the getClimateWizardData
#' function). To do this, the required variables should be provided as "required_variables" parameter.
#' If there is no column "GCM" in the data element of the scenario, then the check_scenario_type
#' parameter should be set to FALSE.
#' 
#' @param temperature_scenario can be one of two options:
#' 1) a data.frame with two columns Tmin and Tmax and n_intervals (default: 12) rows containing
#' temperature changes for all time intervals, or absolute temperatures for these intervals.
#' 2) a temperature scenario object, consisting of the following elements: 'data' = a data frame with
#' n_intervals elements containing the absolute or relative temperature information (as in input option 1);
#' 'scenario_year' = the year the scenario is representative of; 'reference_year' = the year the scenario is representative of; 'scenario_type' = the scenario type
#' ('absolute' or 'relative' - if NA, this is assigned automatically); 'labels' = and elements
#' attached to the input temperature_scenario as an element names 'labels'. A subset of these elements can
#' also be specified, but 'data' must be present.
#' @param n_intervals the number of time intervals specified in the temperature scenarios. This
#' is often the number of months in a year, so the default is 12. If the temperature scenario is
#' specified for a different number of time intervals, this should be adjusted.
#' @param check_scenario_type boolean variable indicating whether the specified (or unspecified)
#' scenario type should be verified, i.e. whether the scenario is a relative or absolute temperature
#' scenario.
#' @param scenario_check_thresholds vector with two numeric elements specifying the thresholds
#' for checking whether the scenario is an absolute or relative temperature scenario. These are the
#' minimum (first value) and maximum (second value) plausible changes in a relative temperature
#' scenario. The test only works in settings where either the lowest mean minimum temperature
#' across all time intervals is below the stated minimum threshold or the highest mean maximum
#' temperature across all time intervals is above the maximum threshold. With the default values
#' c(-5,10), this should be the case for most locations on Earth, but in extreme cases (either for
#' extreme change scenarios or where all monthly minimum and maximum temperatures are between -5
#' and 10 degrees), this may need adjustment. This is only used if check_scenario_type==TRUE.
#' @param update_scenario_type boolean variable stating whether, if scenario type is found to be
#' inconsistent with the numbers, the scenario_type should be updated. Defaults to TRUE and is only
#' used if check_scenario_type==TRUE.
#' @param warn_me boolean variable specifying whether warnings should be shown. Defaults to TRUE.
#' @param required_variables character vector containing the names of columns that are required.
#' This defaults to c("Tmin","Tmax").
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
#' temperature_scenario<-list(data=data.frame(Tmin=c(-5,-2,0, 4, 9,12,15,13,12, 9, 4,0),
#'                                        Tmax=c( 0, 4,8,12,15,18,21,19,17,14,11,5)),
#'                                        reference_year=1975,scenario_type="absolute",
#'                                        labels=list(GCM="none",RCM="none",Time="1950-2000"))
#' 
#' checked_temperature_scenario<-check_temperature_scenario(temperature_scenario,n_intervals=12,
#'    check_scenario_type=FALSE,scenario_check_thresholds=c(-5,10),update_scenario_type=FALSE)
#'                                             
#' checked_temperature_scenario<-check_temperature_scenario(temperature_scenario,n_intervals=12,
#'    check_scenario_type=TRUE,scenario_check_thresholds=c(-5,10),update_scenario_type=FALSE)
#'                                             
#'checked_temperature_scenario<-check_temperature_scenario(temperature_scenario,n_intervals=12,
#'   check_scenario_type=TRUE,scenario_check_thresholds=c(-5,10),update_scenario_type=TRUE)
#'                                             
#'                                             
#'  
#' @export check_temperature_scenario
check_temperature_scenario<-function(temperature_scenario,n_intervals=12,check_scenario_type=TRUE,
                                     scenario_check_thresholds=c(-5,10),update_scenario_type=TRUE,
                                     warn_me=TRUE,required_variables=c("Tmin","Tmax"))
{
  assertthat::assert_that(!is.null(temperature_scenario),msg="temperature_scenario is NULL")
  
  #if(is.null(temperature_scenario)) stop("temperature_scenario is NULL",call. = FALSE)
  
  if(length(which(required_variables %in% names(temperature_scenario)))==length(required_variables))
    {if (warn_me) warning(paste("scenario doesn't contain named elements - consider using the",
                  "following element names: 'data', 'reference_year','scenario_type','labels'"),call. = FALSE)
     data<-temperature_scenario
     scen_year<-NA
     ref_year<-NA
     scen_type<-NA
     labels<-NA} else
       {if("data" %in% names(temperature_scenario))
          data<-temperature_scenario$data else
            stop("scenario isn't an unnamed temperature scenario or provided as 'data' element in temperature_scenario",call. = FALSE)
          if("scenario_year" %in% names(temperature_scenario))
            scen_year<-temperature_scenario$scenario_year else
              {scen_year<-NA
               if (warn_me) warning("no 'scenario_year' element provided in temperature scenario",call. = FALSE)}
          if("reference_year" %in% names(temperature_scenario))
            ref_year<-temperature_scenario$reference_year else
              {ref_year<-NA
               if (warn_me) warning("no 'reference_year' element provided in temperature scenario",call. = FALSE)}
          if("scenario_type" %in% names(temperature_scenario))
            scen_type<-temperature_scenario$scenario_type else
              {scen_type<-NA
               if (warn_me) warning("no 'scenario_type' element provided in temperature scenario",call. = FALSE)}
          if("labels" %in% names(temperature_scenario))
            labels<-temperature_scenario$labels else
              {labels<-NA
               if (warn_me) warning("no 'labels' element provided in temperature scenario",call. = FALSE)}}

   if(!is.data.frame(data))
     stop("specified temperature scenario not provided in valid format - not a data frame",call. = FALSE)
   if(!(length(which(required_variables %in% colnames(data)))==length(required_variables)))
     stop("specified temperature scenario not provided in valid format - at least one of the following columns is missing: ",toString(required_variables),call. = FALSE)
   if(!(is.numeric(unlist(data[,required_variables]))))
     stop("specified temperature scenario not provided in valid format - one of the following columns is not numeric: ",toString(required_variables),call. = FALSE)
   if(!(sum(is.na(unlist(data[,required_variables])))==0))
     stop("specified temperature scenario contains NA values",call. = FALSE)
   if(!"GCM" %in% colnames(data))
     if(!(nrow(data)==n_intervals))
       stop(paste("wrong number of time intervals in temperature scenario - should be",n_intervals),call. = FALSE)
  
   if(!scen_type %in% c("relative","absolute",NA))
     stop("scenario_type must be either 'relative' or 'absolute'",call. = FALSE)
   
   if("GCM" %in% colnames(data)) check_scenario_type<-FALSE

   if(check_scenario_type)
     {if (max(data[,"Tmax"])<=scenario_check_thresholds[2]&
         min(data[,"Tmin"])>=scenario_check_thresholds[1])
            type_guess<-"relative" else type_guess<-"absolute"
      if(!is.na(scen_type))
        {if(!type_guess==scen_type)
           {if(scen_type=="relative")
             if (warn_me) warning(paste("scenario_type doesn't look right - is this really a relative scenario?"),call. = FALSE)
           if(scen_type=="absolute")
             if (warn_me) warning(paste("scenario_type doesn't look right - is this really an absolute scenario?"),call. = FALSE)
           if(update_scenario_type)
             {scen_type<-type_guess
             if(type_guess=="relative")
               {if (warn_me) warning(paste("updating scenario_type to 'relative'"),call. = FALSE)} else
                 if(warn_me) warning(paste("updating scenario_type to 'absolute'"),call. = FALSE)
             }
         }} else
           {scen_type<-type_guess
            if(type_guess=="relative")
             {if (warn_me) warning(paste("setting scenario_type to 'relative'"),call. = FALSE)} else
                if (warn_me) warning(paste("setting scenario_type to 'absolute'"),call. = FALSE)
           }
     }
       
    return(list(data=data,
                scenario_year=scen_year,
                reference_year=ref_year,
                scenario_type=scen_type,
                labels=labels))
}
