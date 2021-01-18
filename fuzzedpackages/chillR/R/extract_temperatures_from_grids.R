#' Extract temperature information from gridded dataset
#' 
#' Temperature data is often available in gridded format, and records for particular points must
#' be extracted for work on site-specific issues (such as chill calculation). This function
#' implements this, for certain types of gridded data.
#' 
#' @param coordinates numeric vector specifying coordinates for the point location of interest.
#' These coordinates have to use the same coordinate system as the grids, from which data are to
#' be extracted. The elements can be named as 'longitude' and 'latitude', or provided as unnamed
#' elements. In the latter case, the first element is interpreted as the x-coordinate (e.g. longitude
#' or Easting) and the second element as the y-coordinate (e.g. latitude or Northing).
#' @param grid_format character string specifying the type of raster data. See details below.
#' @param grid_specifications list of specifications that instruct the function on where to find
#' the temperature grids. See grid_format descriptions for what is required here.
#' @param scenario_year year the temperature scenario is representative of, e.g. 2050, 2080. If the
#' scenario period is an interval, this should be the median of all years in this interval.
#' @param reference_year year of reference for the gridded climate data. This is only important
#' for relative temperature scenarios. If the reference period is an interval, this should be the
#' median of all years in this interval.
#' @param scenario_type character string specifying whether the climate data contains a relative
#' or absolute temperature sceanario. Accordingly, this should be 'relative' or 'absolute'. Can also
#' be NA, which is the default, in which case the function makes a guess on which type applies. This
#' guess is directed by the temperature_check_args.
#' @param labels list of labels to be passed to the labels argument of the resulting temperature_scenario
#' @param  temperature_check_args list of arguments to be passed to the check_temperature_scenario function.
#' Check documentation of that function for details.
#' 
#' @details The following climate data formats are supported: "AFRICLIM" - data downloaded from
#' https://www.york.ac.uk/environment/research/kite/resources/; "CCAFS" - data downloaded from
#' http://ccafs-climate.org/data_spatial_downscaling/; "WorldClim" - data downloaded from
#' http://www.worldclim.org/. All these databases provide separate zipped files for monthly minimum
#' and monthly maximum temperatures, but they differ slightly in format and structure. If you want to
#' see additional formats included, please send me a message.
#' 
#' 
#' @return temperature scenario object extracted from the grids, consisting of the following elements:
#' 'data' = a data frame with n_intervals elements containing the absolute or relative temperature
#' information. 'reference_year' = the year the scenario is representative of. 'scenario_type' =
#' the scenario type ('absolute' or 'relative'); 'labels' = and elements attached to the input
#' temperature_scenario as an element names 'labels'.
#' 
#' The function generates errors, when problems arise.
#' 
#' @importFrom raster raster extract
#' 
#' @author Eike Luedeling
#' @keywords utility
#' @examples
#' 
#'   coordinates<-c(10.6082,34.9411)
#'  # grid_specifications<-list(base_folder="D:/DATA/AFRICLIM/GeoTIFF_30s/future_scenarios/",
#'  #                           minfile="tasmin_rcp45_2055_CCCma-CanESM2_CCCma-CanRCM4_wc30s.zip",
#'  #                           maxfile="tasmax_rcp45_2055_CCCma-CanESM2_CCCma-CanRCM4_wc30s.zip")
#'                             
#'  # extract_temperatures_from_grids(coordinates,grid_format="AFRICLIM",grid_specifications,
#'  #    scenario_type="relative",scenario_year=2055)
#'                  
#'  # grid_specifications<-list(base_folder="D:/DATA/CCAFS_climate/",
#'  #                           minfile="bcc_csm1_1_rcp2_6_2030s_tmin_30s_r1i1p1_b4_asc.zip",
#'  #                           maxfile="bcc_csm1_1_rcp2_6_2030s_tmax_30s_r1i1p1_b4_asc.zip")
#'  #temps<-extract_temperatures_from_grids(coordinates,grid_format="CCAFS",grid_specifications,
#'  #                                       scenario_type="relative",scenario_year=2035)
#'  
#' @export extract_temperatures_from_grids
extract_temperatures_from_grids<-function(coordinates,grid_format,grid_specifications,scenario_year=NA,
                                          reference_year=NA,scenario_type=NA,labels=NA,
                                          temperature_check_args=NULL)
{
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
  
  if(is.null(coordinates)) stop("no coordinates provided",call. = FALSE)
  if(length(coordinates)<2) stop("no useable coordinates provided",call. = FALSE)

  if("longitude" %in% names(coordinates))
    longitude<-as.numeric(coordinates["longitude"]) else
      longitude<-as.numeric(coordinates[1])
  if("latitude" %in% names(coordinates))
    latitude<-as.numeric(coordinates["latitude"]) else
      as.numeric(latitude<-coordinates[2])
  
  if(!is.numeric(c(longitude,latitude))) stop("coordinates not numeric",call. = FALSE)
  
  if(grid_format %in% c("AFRICLIM","CCAFS","WorldClim"))
    {wd<-getwd()
    if(!dir.exists(grid_specifications$base_folder)) stop("base folder ",grid_specifications$base_folder," not found",call. = FALSE)
    setwd(grid_specifications$base_folder)
    if(!file.exists(grid_specifications$minfile)) stop("minimum temperature file",grid_specifications$minfile," not found",call. = FALSE)
    if(!file.exists(grid_specifications$maxfile)) stop("maximum temperature file",grid_specifications$maxfile," not found",call. = FALSE)
  
    if(grid_format=="AFRICLIM")
       {Tmax_files<-unzip(grid_specifications$maxfile)
        Tmax<-ordered_climate_list(Tmax_files,"tif")
        Tmin_files<-unzip(grid_specifications$minfile)
        Tmin<-ordered_climate_list(Tmin_files,"tif")}
    if(grid_format=="CCAFS")
       {Tmax_files<-unzip(grid_specifications$maxfile)
        Tmax<-ordered_climate_list(Tmax_files,"asc")
        dir_name_max<-strsplit(Tmax[1],"/tmax_1.asc")[[1]]
        Tmin_files<-unzip(grid_specifications$minfile)
        Tmin<-ordered_climate_list(Tmin_files,"asc")
        dir_name_min<-strsplit(Tmin[1],"/tmin_1.asc")[[1]]}
    if(grid_format=="WorldClim")
       {Tmax_files<-unzip(grid_specifications$maxfile)
        Tmax<-ordered_climate_list(Tmax_files,"tif")
        Tmin_files<-unzip(grid_specifications$minfile)
        Tmin<-ordered_climate_list(Tmin_files,"tif")}
    
    temperatures<-data.frame(Tmin=rep(NA,12),Tmax=rep(NA,12))
    for (i in 1:12)
    {temperatures[i,"Tmin"]<-extract(raster(Tmin[i]),as.matrix(data.frame(x=longitude,y=latitude)))/10
    temperatures[i,"Tmax"]<-extract(raster(Tmax[i]),as.matrix(data.frame(x=longitude,y=latitude)))/10}
    file.remove(Tmin_files)
    file.remove(Tmax_files)
    if(grid_format=="CCAFS") {unlink(dir_name_max,recursive=TRUE);unlink(dir_name_min,recursive=TRUE)}
    
    setwd(wd)} else
               stop("undefined grid format",call. = FALSE)
  
  temperature_scenario<-list(data=temperatures,scenario_year=scenario_year,reference_year=reference_year,scenario_type=scenario_type,
                             labels=labels)
  
  temperature_scenario<-check_temperature_scenario(temperature_scenario,
                                                   n_intervals=temperature_scenario_check_n_intervals,
                                                   check_scenario_type=temperature_scenario_check_check_scenario_type,
                                                   scenario_check_thresholds=temperature_scenario_check_scenario_check_thresholds,
                                                   update_scenario_type=temperature_scenario_check_update_scenario_type,
                                                   warn_me=temperature_scenario_check_warn_me)
  
  return(temperature_scenario)}
