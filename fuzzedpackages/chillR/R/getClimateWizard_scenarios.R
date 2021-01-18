#' Extract mutltiple scenarios from the ClimateWizard database
#' 
#' This function is a wrapper for the getClimateWizardData function to access climate
#' scenario data for a location of interest. Climate model runs are queried
#' and data returned and summarized according to the specified parameters.
#' A number of metrics are available for several climate models, which are listed in
#' https://github.com/CIAT-DAPA/climate_wizard_api.
#' This function can download data for multiple climate scenarios, saving users the
#' effort to retrieve them separately.
#' 
#' Note that this function lacks quality checks. If something goes wrong, you may
#' consider checking individual scenarios with the getClimateWizardData function.
#' 
#' @param coordinates position of the point of interest, specified by a vector
#' with two elements that are called longitude and latitude (e.g. c(longitude=10,
#' latitude=20)).
#' @param scenarios vector of representative concentration pathway scenarios. Can only be
#' "historical", "rcp45" or "rcp85".
#' @param start_years vector of start year of the intervals, for which data is to be summarized.
#' Must be of same length as scenarios.
#' @param end_years vector of end years of the intervals, for which data is to be summarized.
#' Must be of same length as scenarios.
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
#' #example is #d out, because of runtime issues.
#' #getC<-getClimateWizard_scenarios(coordinates=c(longitude=6.99,latitude=50.62),
#' #                                scenarios=c("rcp85","rcp45"),
#' #                                start_years=c(2070,2035),
#' #                                end_years=c(2100,2065),
#' #                                metric=c("monthly_tmean"),
#' #                                GCMs=c("all"))
#' 
#' 
#' @export getClimateWizard_scenarios
getClimateWizard_scenarios<-function(coordinates,scenarios,start_years,end_years,
                               baseline=c(1950,2005),metric="monthly_min_max_temps",GCMs="all")
{
  assertthat::assert_that(length(coordinates)==2,msg="coordinates not of length 2")
  assertthat::assert_that(all(is.numeric(coordinates)),msg="not all coordinates are numeric")
  assertthat::assert_that(all(scenarios %in% c("historical","rcp45","rcp85")),
                          msg="not all scenarios as 'historical', 'rcp45' or 'rcp85'") 
  assertthat::assert_that(is.numeric(start_years),msg="Start years not numeric")
  assertthat::assert_that(is.numeric(end_years),msg="End years not numeric")
  assertthat::assert_that(length(start_years)==length(end_years),
                          msg="Start and end year vectors have different lengths")
  assertthat::assert_that(is.character(metric), msg="metric not a character vector")
  assertthat::assert_that(is.character(GCMs), msg="metric not a character vector")
  assertthat::assert_that(baseline[2]<=2005, msg="baseline period cannot end after 2005")
  assertthat::assert_that(baseline[2]<=1950, msg="baseline period cannot begin before 1950")
  assertthat::assert_that((baseline[2]-baseline[1])>=20, msg="baseline period must span at least 20 years")
  
  
  results<-list()
  
  for(i in 1:length(scenarios))
    results[[i]]<-getClimateWizardData(coordinates=coordinates,
                                           scenario=scenarios[i],
                                           start_year=start_years[i],
                                           end_year=end_years[i],
                                           baseline=baseline,
                                           metric=metric,
                                           GCMs=GCMs)
  return(results)
}
  

  
#' Plot mutltiple ClimateWizard scenarios obtained with getClimateWizard_scenarios
#' 
#' This function plots multiple scenarios obtained with the getClimateWizard_scenarios
#' function.
#' 
#' @param getscenarios_element outputs from the getClimateWizard_scenarios function
#' @param low_filter numeric value specifying the lowest plausible value for the
#' variable of interest. This is sometimes necessary to exclude erroneous values
#' in the ClimateWizard database. 
#' @param high_filter numeric value specifying the highest plausible value for the
#' variable of interest. This is sometimes necessary to exclude erroneous values
#' in the ClimateWizard database. 
#' @param color color to be used for the plots.
#' @import ggplot2
#' 
#' @return returns nothing, but a plot is produced as a side effect.
#' 
#' @references Girvetz E, Ramirez-Villegas J, Navarro C, Rodriguez C, Tarapues J, undated.
#' ClimateWizard REST API for querying climate change data.
#' https://github.com/CIAT-DAPA/climate_wizard_api
#' 
#' 
#' @author Eike Luedeling
#' @keywords utility
#' @examples
#' 
#' #example is #d out, because of runtime issues.
#' #getC<-getClimateWizard_scenarios(coordinates=c(longitude=6.99,latitude=50.62),
#' #                                scenarios=c("rcp45","rcp45","rcp85","rcp85"),
#' #                                start_years=c(2035,2070,2035,2070),
#' #                                end_years=c(2065,2100,2065,2100),
#' #                                metric=c("monthly_tmean"),
#' #                                GCMs=c("all"))
#' #plot_climateWizard_scenarios(getC,low_filter=-6,high_filter=6,color="red")
#' 
plot_climateWizard_scenarios<-function(getscenarios_element,low_filter=-1000,high_filter=1000,color="cadetblue")
{
  allfuture<-reshape2::melt(getscenarios_element[[1]]$data)
  
  allfuture<-allfuture[which(allfuture$value>low_filter),]
  allfuture<-allfuture[which(allfuture$value<high_filter),]
  allfuture[,"Year"]<-median(c(getscenarios_element[[1]]$start_year,getscenarios_element[[1]]$end_year))
  allfuture[,"RCP"]<-getscenarios_element[[1]]$scenario
  
  if (length(getscenarios_element)>1)
      for (i in 2:length(getscenarios_element))
      { 
        data<-reshape2::melt(getscenarios_element[[i]]$data)
        data<-data[which(data$value>low_filter),]
        data<-data[which(data$value<high_filter),]
        data[,"Year"]<-median(c(getscenarios_element[[i]]$start_year,getscenarios_element[[i]]$end_year))
        data[,"RCP"]<-getscenarios_element[[i]]$scenario
        allfuture<-rbind(allfuture,data)
      }
  
  #the following two lines are necessary to make the ggplot call pass the CRAN checks
  Year<-9999
  RCP<-"funny"
  
  varname<-identify_common_string(as.character(allfuture$variable))
  allfuture$variable<-sapply(allfuture$variable,function(x) {as.numeric(strsplit(as.character(x),varname)[[1]][2])} )
  
  ggplot(data=allfuture, aes_(factor(~variable),~value)) + 
    geom_violin(mapping = NULL, data = NULL, stat = "ydensity",
                draw_quantiles = NULL, trim = TRUE,
                scale = "width", na.rm = FALSE, show.legend = NA,
                inherit.aes = TRUE,fill=color) +
    geom_jitter(height = 0, width = 0.1) +facet_grid(cols=vars(Year),rows=vars(RCP))+
    theme_grey(base_size = 22)+
    labs(x="Month",y="Projected change")
 
}

