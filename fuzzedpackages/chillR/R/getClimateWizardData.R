#' Extract climate data from the ClimateWizard database
#' 
#' This function makes use of an API provided by the International
#' Center for Tropical Agriculture (CIAT) to access climate scenario
#' data for a location of interest. Climate model runs are queried
#' and data returned and summarized according to the specified parameters.
#' A number of metrics are available for several climate models, which are listed in
#' https://github.com/CIAT-DAPA/climate_wizard_api.
#' Refer to this document for details on what can be downloaded. This function
#' provides the additional option of automatically retrieving all data referring to changes
#' in daily temperature extremes (by month), by setting the ```metric``` parameter to
#' "monthly_min_max_temps". It also offers the option to automatically obtain data for
#' all climate models included in the database (as of January 2018).
#' #' 
#' @param coordinates position of the point of interest, specified by a vector
#' with two elements that are called longitude and latitude (e.g. c(longitude=10,
#' latitude=20)).
#' @param scenario representative concentration pathway scenario. Can only be
#' "historical", "rcp45" or "rcp85".
#' @param start_year start year of the interval, for which data is to be summarized.
#' @param end_year end year of the interval, for which data is to be summarized.
#' @param baseline numeric vector of length 2 indicating the time interval to be used
#' as baseline for the climate scenario. The function then returns projected values relative
#' to this baseline. Defaults to c(1950,2005) for the standard
#' baseline of the ClimateWizard dataset. This can also assume different values, but it must
#' span an interval of at least 20 years within the [1950; 2005] interval. Needs
#' to be set to NA for the function to return absolute values.
#' @param metric vector of metrics to output, from a list specified in the reference
#' provided above. This can also be "monthly_min_max_temps", which returns all
#' mean monthly minimum and maximum temperatures, or "precipitation" for precipitation
#' data for all months, or "monthly_tmean" for the mean monthly temperatures of all months.
#' @param GCMs vector of GCMs to be accessed, from a list specified in the above
#' reference. This can also be "all" for all available GCMs (as of January 2018).
#' @param temperature_generation_scenarios parameter to indicate whether the scenarios to be
#' generated should be formatted in such a way that they are direclty usable by
#' chillR's temperature_generation function. This is only applicable, when metric==
#' 'monthly_min_max_temps'.
#' @return data.frame containing the requested information.
#' 
#' @references Girvetz E, Ramirez-Villegas J, Navarro C, Rodriguez C, Tarapues J, undated.
#' ClimateWizard REST API for querying climate change data.
#' https://github.com/CIAT-DAPA/climate_wizard_api
#' 
#' @author Eike Luedeling
#' @keywords utility
#' @examples
#' 
#' # the example is #d out, since the download request sometimes times out, and that
#' # causes problems with CRAN approval of the package
#' 
#' # getClimateWizardData(coordinates=c(longitude=10.613975,latitude=34.933439),
#' #   scenario="rcp45", start_year=2020, end_year=2050,
#' #   metric=c("CD18","R02"), GCMs=c("bcc-csm1-1","BNU-ESM"))
#' 
#' 
#' @export getClimateWizardData
getClimateWizardData<-function(coordinates,scenario,start_year,end_year,
                               baseline=c(1950,2005),metric="monthly_min_max_temps",GCMs="all",
                               temperature_generation_scenarios=FALSE)
{
  assertthat::assert_that(length(coordinates)==2,msg="coordinates not of length 2")
  assertthat::assert_that(all(is.numeric(coordinates)),msg="not all coordinates are numeric")
  assertthat::assert_that(all(scenario %in% c("historical","rcp45","rcp85")),
                          msg="scenario must be 'historical', 'rcp45' or 'rcp85'") 
  assertthat::assert_that(is.numeric(start_year),msg="Start year not numeric")
  assertthat::assert_that(is.numeric(end_year),msg="End year not numeric")
  assertthat::assert_that(length(start_year)==1 & length(end_year)==1,
                          msg="Start or end year object has too many elements")
  assertthat::assert_that(is.character(metric), msg="metric not a character vector")
  assertthat::assert_that(is.character(GCMs), msg="metric not a character vector")
  assertthat::assert_that(baseline[2]<=2005, msg="baseline period cannot end after 2005")
  assertthat::assert_that(baseline[1]>=1950, msg="baseline period cannot begin before 1950")
  assertthat::assert_that((baseline[2]-baseline[1])>=20, msg="baseline period must span at least 20 years")
  
  
  
  coordinates_usable<-FALSE
  if(is.null(names(coordinates)))
    if(length(coordinates)==2) {longitude<-coordinates[1]
    latitude<-coordinates[2]
    warning("Coordinates not named. Interpreting the first ",
            "element (",longitude,") as longitude and the second (",
            latitude,") as latitude.")
    coordinates_usable<-TRUE}
  if(!is.null(names(coordinates)))
  {if(!length(coordinates)==2)
  {warning("Unable to interpret coordinates parameter. ",
           "Needs to be a vector with two numeric elements, ",
           "ideally called 'x' and 'y' or 'longitude' and 'latitude'")
    return()}
    if(length(which(c("latitude","longitude") %in% names(coordinates)))==2)
    {longitude<-as.numeric(coordinates["longitude"])
    latitude<-as.numeric(coordinates["latitude"])
    coordinates_usable<-TRUE}
    if(length(which(c("x","y") %in% names(coordinates)))==2)
    {longitude<-as.numeric(coordinates["x"])
    latitude<-as.numeric(coordinates["y"])
    coordinates_usable<-TRUE}
    if(length(which(c("long","lat") %in% names(coordinates)))==2)
    {longitude<-as.numeric(coordinates["long"])
    latitude<-as.numeric(coordinates["lat"])
    coordinates_usable<-TRUE}
  }
  if(coordinates_usable)
  {if(is.na(longitude)) coordinates_usable<-FALSE
  if(is.na(latitude)) coordinates_usable<-FALSE}
  if(!coordinates_usable)
  {warning("Unable to interpret coordinates parameter. ",
           "Needs to be a vector with two numeric elements, ",
           "ideally called 'x' and 'y' or 'longitude' and 'latitude'")
    return()}
  
  if(!scenario %in% c("historical","rcp45","rcp85"))
  {warning('Scenario not available. Must be "historical", "rcp45" or "rcp85".')
    return()}
  
  all_available_gcms<-c("bcc-csm1-1","BNU-ESM","CanESM2","CESM1-BGC","MIROC-ESM","CNRM-CM5","ACCESS1-0","CSIRO-Mk3-6-0",
                        "GFDL-CM3","GFDL-ESM2G","GFDL-ESM2M","inmcm4",
                        "IPSL-CM5A-LR","IPSL-CM5A-MR",
                        "CCSM4")
  if(GCMs[1]=="all") gcms<-all_available_gcms else
  {if(length(which(!GCMs %in% all_available_gcms))>0)
  {warning(GCMs[which(!GCMs %in% all_available_gcms)]," not contained in the database (as far as the getClimateWizardData knows). Models that are available are ",
           toString(all_available_gcms))
    return()}
    
    gcms<-GCMs}
  
  metric_GCMs<-data.frame(GCM=gcms,metric_available=NA)
  
  for(gcm in gcms)
  {
    if(metric[1]=="monthly_min_max_temps")
    {available<-TRUE
    for(Textreme in c("tasmin","tasmax"))
      for(i in 1:12)
      {if(available)
      {tmetric<-paste(Textreme,i,sep="")
      if(!is.na(baseline[1]))
        checkgcm<-jsonlite::fromJSON(paste("http://climatewizard.ccafs-climate.org/service?lat=",
                                           latitude,"&lon=",longitude,
                                           "&index=",tmetric,"&scenario=",scenario,"&gcm=",gcm,"&range=",
                                           start_year,"-",end_year,"&baseline=",baseline[1],"-",baseline[2],"&avg=true",sep="")
        ) else
          checkgcm<-jsonlite::fromJSON(paste("http://climatewizard.ccafs-climate.org/service?lat=",
                                             latitude,"&lon=",longitude,
                                             "&index=",tmetric,"&scenario=",scenario,"&gcm=",gcm,"&range=",
                                             start_year,"-",end_year,"&avg=true",sep=""))
        if(length(checkgcm)==1) {if(!is.null(checkgcm$error)) available<-FALSE}
        if (available)
          if(checkgcm$values[[2]][1]=="out of range")
          {warning("Time interval not (fully) contained in the dataset - no data retrieved")}
        if(available) metric_GCMs[which(metric_GCMs$GCM==gcm),tmetric]<-as.numeric(checkgcm$values[[2]][1])}
      }
    } 
    
    if(metric[1]=="precipitation")
    {available<-TRUE
    for(i in 1:12)
    {if(available)
    {tmetric<-paste("pr",i,sep="")
    if(!is.na(baseline[1]))
      checkgcm<-jsonlite::fromJSON(paste("http://climatewizard.ccafs-climate.org/service?lat=",
                                         latitude,"&lon=",longitude,
                                         "&index=",tmetric,"&scenario=",scenario,"&gcm=",gcm,"&range=",
                                         start_year,"-",end_year,"&baseline=",baseline[1],"-",baseline[2],"&avg=true",sep="")
      ) else
        checkgcm<-jsonlite::fromJSON(paste("http://climatewizard.ccafs-climate.org/service?lat=",
                                           latitude,"&lon=",longitude,
                                           "&index=",tmetric,"&scenario=",scenario,"&gcm=",gcm,"&range=",
                                           start_year,"-",end_year,"&avg=true",sep=""))
      if(length(checkgcm)==1) {if(!is.null(checkgcm$error)) available<-FALSE}
      if (available)
        if(checkgcm$values[[2]][1]=="out of range")
        {warning("Time interval not (fully) contained in the dataset - no data retrieved")}
      if(available) metric_GCMs[which(metric_GCMs$GCM==gcm),tmetric]<-as.numeric(checkgcm$values[[2]][1])}
    }
    }
    
    if(metric[1]=="monthly_tmean")
    {available<-TRUE
    for(i in 1:12)
    {if(available)
    {tmetric<-paste("tas",i,sep="")
    if(!is.na(baseline[1]))
      checkgcm<-jsonlite::fromJSON(paste("http://climatewizard.ccafs-climate.org/service?lat=",
                                         latitude,"&lon=",longitude,
                                         "&index=",tmetric,"&scenario=",scenario,"&gcm=",gcm,"&range=",
                                         start_year,"-",end_year,"&baseline=",baseline[1],"-",baseline[2],"&avg=true",sep="")
      ) else
        checkgcm<-jsonlite::fromJSON(paste("http://climatewizard.ccafs-climate.org/service?lat=",
                                           latitude,"&lon=",longitude,
                                           "&index=",tmetric,"&scenario=",scenario,"&gcm=",gcm,"&range=",
                                           start_year,"-",end_year,"&avg=true",sep=""))
      if(length(checkgcm)==1) {if(!is.null(checkgcm$error)) available<-FALSE}
      if (available)
        if(checkgcm$values[[2]][1]=="out of range")
        {warning("Time interval not (fully) contained in the dataset - no data retrieved")}
      if(available) metric_GCMs[which(metric_GCMs$GCM==gcm),tmetric]<-as.numeric(checkgcm$values[[2]][1])}
    }
    }
    
    
    
    if(!metric[1] %in% c("precipitation","monthly_min_max_temps","monthly_tmean"))
    {for (met in metric)
    {available=TRUE
    if(!is.na(baseline[1]))
      checkgcm<-jsonlite::fromJSON(paste("http://climatewizard.ccafs-climate.org/service?lat=",
                                         latitude,"&lon=",longitude,
                                         "&index=",met,"&scenario=",scenario,"&gcm=",gcm,"&range=",
                                         start_year,"-",end_year,"&baseline=",baseline[1],"-",baseline[2],"&avg=true",sep="")) else
                                           checkgcm<-jsonlite::fromJSON(paste("http://climatewizard.ccafs-climate.org/service?lat=",
                                                                              latitude,"&lon=",longitude,
                                                                              "&index=",met,"&scenario=",scenario,"&gcm=",gcm,"&range=",
                                                                              start_year,"-",end_year,"&avg=true",sep=""))
                                         
                                         if(length(checkgcm)==1) {if(!is.null(checkgcm$error)) available<-FALSE}
                                         if (available) if(checkgcm$values[[2]][1]=="out of range") {warning("Time interval not (fully) contained in the dataset - no data retrieved")}
                                         if(available) metric_GCMs[which(metric_GCMs$GCM==gcm),met]<-as.numeric(checkgcm$values[[2]][1])
    }
    }
    if(length(metric)==1) metric_GCMs[which(metric_GCMs$GCM==gcm),"metric_available"]<-available
  }
  #if(is.na(metric_GCMs$metric_available[1]))
  metric_GCMs<-metric_GCMs[,-2]
  
  checkmetrics<-metric
  if(checkmetrics[1]=="monthly_min_max_temps")
    checkmetrics<-c(paste("tasmin",1:12,sep=""),paste("tasmax",1:12,sep=""))
  if(checkmetrics[1]=="precipitation")
    checkmetrics<-c(paste("pr",1:12,sep=""))
  if(checkmetrics[1]=="monthly_tmean")
    checkmetrics<-c(paste("tas",1:12,sep=""))
  
  if(length(which(!checkmetrics %in% colnames(metric_GCMs))>0))
    warning("Data for ",toString(checkmetrics[which(!checkmetrics %in% colnames(metric_GCMs))])," not available")
  
  
  output<-list(data=metric_GCMs,scenario=scenario,start_year=start_year,end_year=end_year,scenario_year=median(c(start_year,end_year)),
               reference_year=NA,scenario_type="absolute",labels=list())
  if(!is.na(baseline[1])) {output$reference_year<-median(baseline[1]:baseline[2])
  output$scenario_type<-"relative"}
  
  if(metric[1]=="monthly_min_max_temps" & temperature_generation_scenarios)
  {outputs<-list()
  for (i in 1:nrow(output$data))
  {outputs[[i]]<-output
  outputs[[i]]$data<-data.frame(Tmin=as.numeric(output$data[i,2:13]),Tmax=as.numeric(output$data[i,14:25]))
  outputs[[i]]$labels<-as.character(output$data[i,1])
  }
  names(outputs)<-output$data$GCM
  output<-outputs
  }
  
  return(output)
}
