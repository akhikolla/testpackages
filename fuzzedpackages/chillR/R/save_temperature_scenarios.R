#' Save temperature scenarios generated with temperature_generation
#' 
#' The temperature_generation can produce synthetic temperature scenarios, but it can take a while
#' to run, especially for large ensembles of climate scenarios. The save_temperature_scenarios
#' function can then save these scenarios to disk as a series of .csv files, so that they can
#' later be used again, without re-running the generation function.
#' Conversely, the load_temperature_scenarios function allows reading the data back into R.
#' This function also works with any other list of data.frames.
#' 
#' @param generated_temperatures list of temperature scenarios produced with the
#' temperature_generation function.
#' @param path character string indicating the file path where the files are to be written.
#' @param prefix character string specifying the prefix for all files.
#' @return no values are returned, but files are written as a side_effect.
#' 
#' @author Eike Luedeling
#' @keywords utility
#' @examples
#' 
#' temps<-list(Element1=data.frame(a=1,b=2),Element2=data.frame(a=c(2,3),b=c(8,4)))
#' # save_temperature_scenarios(temps,path=getwd(),prefix="temperatures")
#' # temps_reloaded<-load_temperature_scenarios(path=getwd(),prefix="temperatures")
#' 
#'  
#' @export save_temperature_scenarios
save_temperature_scenarios<-function(generated_temperatures,path,prefix)
{
  if(!dir.exists(path)) dir.create(path)
  if(is.list(generated_temperatures))
    for (i in 1:length(generated_temperatures))
      write.csv(generated_temperatures[[i]],
                file.path(path,
                          paste(prefix,"_",i,"_",names(generated_temperatures)[[i]],".csv",sep="")),
                row.names=FALSE)
}
  
#' Load temperature scenarios
#' 
#' The temperature_generation can produce synthetic temperature scenarios, but it can take a while
#' to run, especially for large ensembles of climate scenarios. The save_temperature_scenarios
#' function can then save these scenarios to disk as a series of .csv files, so that they can
#' later be used again, without re-running the generation function.
#' Conversely, the load_temperature_scenarios function allows reading the data back into R.
#' This function also works with any other list of data.frames.
#' 
#' @param path character string indicating the file path where the files are to be written.
#' @param prefix character string specifying the prefix for all files.
#' @return a list of temperature scenarios.
#' 
#' @author Eike Luedeling
#' @keywords utility
#' @examples
#' 
#' temps<-list(Element1=data.frame(a=1,b=2),Element2=data.frame(a=c(2,3),b=c(8,4)))
#' # save_temperature_scenarios(temps,path=getwd(),prefix="temperatures")
#' # temps_reloaded<-load_temperature_scenarios(path=getwd(),prefix="temperatures")
#' 
#'  
#' @export load_temperature_scenarios
load_temperature_scenarios<-function(path,prefix)
{
  files<-list.files(path)
  file_prefixes<-lapply(files,function(x) substr(x,1,nchar(prefix)))
  scenario_files<-files[which(file_prefixes==prefix)]
  scenario_num_names<-lapply(scenario_files,function(x) substr(x,nchar(prefix)+2,nchar(x)))
  nums<-unlist(lapply(scenario_num_names,
             function(x) as.numeric(strsplit(x,"_")[[1]][1])))
  ordered_scenarios<-scenario_files[order(nums)]
  ordered_scenario_names<-scenario_num_names[order(nums)]
  scennames<-lapply(ordered_scenario_names,
                function(x)
                  {num<-strsplit(x,"_")[[1]][1]
                   substr(x,nchar(num)+2,nchar(x)-4)})
  output<-list()
  for(i in 1:length(ordered_scenarios))
    output[[i]]<-read.csv(file.path(path,ordered_scenarios[[i]]))
  names(output)<-scennames
  return(output)
}

#' Load climate wizard scenarios
#' 
#' This is a slightly modified version of the load_temperature_scenarios function that can
#' load climate scenarios downloaded with the getClimateWizardData and saved with the
#' save_temperature_scenarios function. This separate function is necessary, because the
#' climate scenarios are expressed as lists, with one element being a data.frame.
#' 
#' @param path character string indicating the file path where the files are to be written.
#' @param prefix character string specifying the prefix for all files.
#' @return a list of temperature scenarios.
#' 
#' @author Eike Luedeling
#' @keywords utility
#' @examples
#' 
#' temps<-list(Element1=data.frame(a=1,b=2),Element2=data.frame(a=c(2,3),b=c(8,4)))
#' # save_temperature_scenarios(temps,path=getwd(),prefix="temperatures")
#' # temps_reloaded<-load_temperature_scenarios(path=getwd(),prefix="temperatures")
#' 
#'  
#' @export load_ClimateWizard_scenarios
load_ClimateWizard_scenarios<-function(path,prefix)
{
  files<-list.files(path)
  file_prefixes<-lapply(files,function(x) substr(x,1,nchar(prefix)))
  scenario_files<-files[which(file_prefixes==prefix)]
  scenario_num_names<-lapply(scenario_files,function(x) substr(x,nchar(prefix)+2,nchar(x)))
  nums<-unlist(lapply(scenario_num_names,
                      function(x) as.numeric(strsplit(x,"_")[[1]][1])))
  ordered_scenarios<-scenario_files[order(nums)]
  ordered_scenario_names<-scenario_num_names[order(nums)]
  scennames<-lapply(ordered_scenario_names,
                    function(x)
                    {num<-strsplit(x,"_")[[1]][1]
                    substr(x,nchar(num)+2,nchar(x)-4)})
  output<-list()
  for(i in 1:length(ordered_scenarios))
    output[[i]]<-read.csv(file.path(path,ordered_scenarios[[i]]))
  names(output)<-scennames
  
  for(o in 1:length(output))
    output[[o]]<-list(data=data.frame(Tmin=output[[o]]$data.Tmin,
                       Tmax=output[[o]]$data.Tmax),
            scenario=as.character(output[[o]]$scenario[1]),
            start_year=output[[o]]$start_year[1],
            end_year=output[[o]]$end_year[1],
            scenario_year=output[[o]]$scenario_year[1],
            reference_year=output[[o]]$reference_year[1],
            scenario_type=as.character(output[[o]]$scenario_type[1]),
            labels=as.character(output[[o]]$labels[1]))
            
            
    
  return(output)
}

