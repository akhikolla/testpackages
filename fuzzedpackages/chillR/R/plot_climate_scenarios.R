#' Plot multiple chilling scenario groups (or for other metrics)
#' 
#' For quantifying climate risks, it is useful to generate many version of plausible weather for
#' particular climate scenarios. This can, for example, be done with the temperature_generation function.
#' This function facilitates illustration of these results by providing various options to show them as
#' boxplots. The function can plot either a single panel of climate scenarios or multiple panels
#' side by side.
#' 
#' @param climate_scenario_list list of lists containing information about the chill scenarios.
#' These lists must have an element named 'data' which should contain a data.frame with a column
#' named the same as the 'metric' argument, which contains (numeric) information to be plotted.
#' Additional optional elements are 'time_series' (Boolean, indicating if a time series is to be plotted),
#' 'labels' (vector of length 'length(data)' containing labels for the scenarios - if this is a time
#' series scenario, these must be numeric; if the data are not a time series, the labels aren't shown in the
#' plot, because there wouldn't normally be enough space - only numbers are shown there, and the legend
#' is provided in the value returned by this function), 'caption' (up to three character strings indicating
#' the caption to be displayed in the respective plot panel; the elements of this vector are displayed on
#' different lines. If caption_above==TRUE, only the first element is displayed) and 'historic_data'
#' (a data.frame containing at least two columns named the same as 'metric' and 'year_name').
#' documentation of 'make_chill_scenario_plot' for details on these.
#' @param metric character string corresponding to the name of the column that contains the data of
#' interest in the climate_scenario_list data.frames (and if applicable the historic_data data.frame).
#' @param metric_label character string specifying the y-axis label.
#' @param year_name character string indicating the name of the time column in the historic_data data.frame.
#' @param label_sides indicates what sides of the plot y-axis labels are to be drawn. Can be "left', "right"
#' or "both". If label_sides assumes any other value, no labels are plotted.
#' @param ylim numeric vector of length 2, specifying the lower and upper limits of the y-axis. If either of
#' these two values is NA, it is automatically selected based on the data range. 
#' @param reference_line numeric vector of length 1, 2 or 3, specifying a horizontal reference bar to be drawn
#' across the plot (e.g. to indicate exceedance of a threshold). A reference_line argument of length 1 is
#' interpreted by drawing a line across the plot at the specified value. If length(reference_line)==2, the
#' values are interpreted as lower and upper limit of a rectangular threshold area. If
#' length(reference_line)==3, the lowest and highest values are used to draw a rectangle and the median value
#' to draw a line (e.g. to show a best estimate and a confidence interval around it).
#' @param col_rect color code or name for the color of the reference_line rectangle.
#' @param col_line color code or name for the color of the reference_line line.
#' @param hist_col color code or name for the color of the historic data points.
#' @param texcex numeric variable indicating character size (cex for all text elements in the plot).
#' @param caption_above Boolean variable indicating whether the caption should be drawn above (TRUE) or
#' inside the figure.
#' @param family character string specifying the font family ('serif', 'sans' or 'mono').
#' @param no_scenario_numbers Boolean variable indicating whether climate scenarios should be numbered
#' in the plot (this can clutter the figure).
#' 
#' @importFrom grDevices col2rgb
#' @importFrom graphics boxplot strheight strwidth
#' 
#' @return List of legends for the different panels of the plot. This list reads 'time series labels'
#' for time series plot, 'no adequate labels provided' for unlabeled collections of boxplot, and a data.frame
#' explaining the number codes used as the legend in labeled collections of boxplots.
#' As a side effect, a plot of the climate scenarios is drawn.
#' 
#' The function generates errors, when problems arise.
#' 
#' @importFrom raster raster extract
#' 
#' @author Eike Luedeling
#' @keywords utility
#' @examples
#'
#'   #making 3 identical objects as scenarios; let's assume these represent the
#'   #years 2000, 2005 and 2010. 
#'   
#'   models<-list(Chilling_Hours=Chilling_Hours,Utah_Chill_Units=Utah_Model,Chill_Portions=
#'      Dynamic_Model,GDH=GDH)
#'   
#'   chill<-tempResponse(stack_hourly_temps(
#'      fix_weather(KA_weather[which(KA_weather$Year>2003),]),latitude=50.4),
#'      Start_JDay = 305,End_JDay = 60,models)
#'   scenario_results<-list(chill,chill,chill)
#'   
#'   climate_scenario_list<-list(list(data=scenario_results,
#'                                    caption=c("Historic","data"),
#'                                    time_series=TRUE,
#'                                    labels=c(2000,2005,2010),
#'                                    historic_data=chill),
#'                               list(data=scenario_results,
#'                                    caption=c("Scenario","1"),
#'                                    labels=c("Climate model 1",
#'                                             "Climate model 2",
#'                                             "Climate model 3")),
#'                               list(data=scenario_results,
#'                                    caption=c("Scenario","2")),
#'                               list(data=scenario_results,
#'                                    caption=c("Scenario","3")))
#'                                    
#'                                    
#'  plot_climate_scenarios(climate_scenario_list,metric="Chill_Portions",
#'                                  metric_label="Chill Portions",
#'                                  year_name="End_year",label_sides="both",
#'                                  reference_line=c(40,45,50),col_rect=NA,col_line=NA,
#'                                  texcex=2,caption_above=FALSE)
#'                                  
#'  plot_climate_scenarios(climate_scenario_list,"Chill_Portions","Chill Portions",
#'     texcex=1)
#'
#'   
#'  
#' @export plot_climate_scenarios
plot_climate_scenarios<-function(climate_scenario_list,metric,metric_label,
                                          year_name="End_year",label_sides="both",ylim=c(0,NA),
                                          reference_line=NULL,col_rect=NA,col_line=NA,hist_col=NA,
                                          texcex=2,caption_above=FALSE,family='serif',no_scenario_numbers=FALSE)
{
  layout(matrix(c(1:length(climate_scenario_list)),ncol=length(climate_scenario_list),byrow=FALSE),
         widths=rep(6,length(climate_scenario_list)))
  
  par(family=family)
  ylimplot<-ylim
  if(is.na(ylimplot[2]))
    ylimplot[2]<-max(sapply(climate_scenario_list,function(x) max(sapply(x[["data"]],function(x) max(x[,metric])))))
  if(is.na(ylimplot[1]))
    ylimplot[1]<-min(sapply(climate_scenario_list,function(x) min(sapply(x[["data"]],function(x) min(x[,metric])))))
  
  text_line_height<-par('cin')[2]*par('cex')*par('lheight')
  
  scenario_caption_max_length<-max(sapply(climate_scenario_list,function(x) if("caption" %in% names(x)) length(x$caption) else 0),na.rm=TRUE)
  
  if(!caption_above)
    {caption_space_needs<-(scenario_caption_max_length+0.5)*texcex*text_line_height
     plot_height<-par("pin")[2]
     ylimplot[2]<-ylimplot[1]+(ylimplot[2]-ylimplot[1])/(plot_height-caption_space_needs)*plot_height}
  
  if(label_sides %in% c("left","right","both"))
    axis_space_need<-strheight(metric_label,cex=texcex,font=2,units="inches")+
    max(strwidth(pretty(ylimplot),units="inches",cex=texcex,font=2))+
    2.5*text_line_height # 2 lines of empty space
  
  paromi_old<-par("omi")
  if(label_sides %in% c("left","both")) left_space<-axis_space_need else left_space<-0
  if(label_sides %in% c("right","both")) right_space<-axis_space_need else right_space<-0
  
  par(omi=c(strheight(metric_label,cex=texcex,font=2,units="inches")+text_line_height*2,left_space,
           strheight(metric_label,cex=texcex,font=2,units="inches")+text_line_height*2,right_space))
  par(mar=c(0,0,0,0))
  left_label<-FALSE; right_label<-FALSE
  if(label_sides=="both") {left_label<-TRUE; right_label<-TRUE}
  if(label_sides=="left") left_label<-TRUE
  if(label_sides=="right") right_label<-TRUE
  
  legend<-list()
  
  for(scen in 1:length(climate_scenario_list))
  {
    if("data" %in% names(climate_scenario_list[[scen]]))
      plotdata<-climate_scenario_list[[scen]]$data else stop("Climate scenario ",scen," doesn't have a data element",call. = FALSE)
    if("caption" %in% names(climate_scenario_list[[scen]]))
      caption<-climate_scenario_list[[scen]]$caption else caption<-NULL
    if("time_series" %in% names(climate_scenario_list[[scen]]))
      time_series<-climate_scenario_list[[scen]]$time_series else time_series<-FALSE
    if("labels" %in% names(climate_scenario_list[[scen]]))
      plotlabels<-climate_scenario_list[[scen]]$labels else plotlabels<-NA
    if("historic_data" %in% names(climate_scenario_list[[scen]]))
      historic_data<-climate_scenario_list[[scen]]$historic_data else historic_data<-NULL
  
    for(i in 1:length(climate_scenario_list))
    {if(is.null(climate_scenario_list[[i]])) stop("no data in climate_scenario_list element ",i,call. = FALSE)
      for(j in 1:length(climate_scenario_list[[i]]$data))
        {if(!is.data.frame(climate_scenario_list[[i]]$data[[j]])) stop("scenario ",j," in climate_scenario_list element ",i," is not a data.frame",call. = FALSE)
         if(!metric %in% colnames(climate_scenario_list[[i]]$data[[j]])) stop("scenario ",j," in climate_scenario_list element ",i," doesn't have a column ",metric,call. = FALSE)
         if(!is.numeric(climate_scenario_list[[i]]$data[[j]][,metric])) stop(metric," data in scenario ",j," in climate_scenario_list element ",i," isn't numeric",call. = FALSE)
    }}
    
    
    
    if(is.na(metric)) stop("no metric for plotting specified",call. = FALSE)
    if(!metric %in% colnames(climate_scenario_list[[1]]$data[[1]])) stop("data.frames for plotting don't have a column ",metric,call. = FALSE)
    
    
    if(time_series) xlimplot<-c(min(plotlabels)-1,max(plotlabels)+1)
    if(!time_series) xlimplot<-c(1-0.5,length(plotdata)+0.5)
    
    if(time_series==TRUE)
      if(!is.null(historic_data))
      {if(!is.data.frame(historic_data)) stop("historic_data argument is not a data.frame",call. = FALSE)
        if(!metric %in% colnames(historic_data)) stop(metric," column missing in historic_data data.frame",call. = FALSE)
        if(!year_name %in% colnames(historic_data)) stop(year_name," column missing in historic_data data.frame",call. = FALSE)
        if(!is.numeric(unlist(historic_data[metric]))) stop(metric," data not numeric",call. = FALSE)
        if(!is.numeric(unlist(historic_data[year_name]))) stop(year_name," data not numeric",call. = FALSE)
        xlimplot[1]<-min(c(xlimplot[1],min(historic_data[year_name])),na.rm=TRUE)
        xlimplot[2]<-max(c(xlimplot[2],max(historic_data[year_name])),na.rm=TRUE)
      }
    
    
    if(!time_series)
      boxplot(lapply(plotdata,function(x) x[[metric]]),
              ylab=NA,ylim=ylimplot,xlim=xlimplot,axes=FALSE,yaxs="i",col="LIGHT YELLOW",lwd=texcex*3/4)
    if(time_series)
    {if(!length(plotdata)==length(plotlabels)) stop("number of points in time doesn't correspond to number of scenarios",call. = FALSE)
      if(!is.numeric(plotlabels)) stop("labels for time series aren't numeric")
      boxplot(lapply(plotdata,function(x) x[[metric]]),
              ylab=NA,ylim=ylimplot,xlim=xlimplot,yaxs="i",las=1,at=as.numeric(plotlabels),
              axes=FALSE,col="LIGHT YELLOW",lwd=texcex*3/4,boxwex=2)}
    
    
    if(!is.null(reference_line))
    {if(!is.numeric(reference_line)) stop("information for reference_line isn't numeric")
      if(is.na(col_rect)) col_rect<-rgb(66/255,72/255,244/255,0.20) else
      {colelem<-as.numeric(col2rgb(col_rect)/255)
      col_rect<-rgb(colelem[1],colelem[2],colelem[3],0.2)}
      if(is.na(col_line)) col_line<-rgb(66/255,72/255,244/255)    
      if(length(reference_line) %in% c(2,3))
        rect(par("usr")[1],sort(reference_line)[1],
             par("usr")[2],sort(reference_line)[3],
             col= col_rect,border=NA)
      if(length(reference_line) %in% c(1,3))
        arrows(x0=par("usr")[1],y0=sort(reference_line)[2],
               x1=par("usr")[2],y1=sort(reference_line)[2],
               lwd=texcex,col=col_line,length=0)}
    
    
    if(time_series)
      if(!is.null(historic_data))
        {if(is.na(hist_col)) hist_col<-"RED"
        points(historic_data[,metric]~historic_data[,year_name],pch=16,cex=0.75,col=hist_col)}
    if(time_series) axis(1,cex.axis=texcex,lwd=texcex,padj=0.5)
    if(!time_series)
      if(length(plotlabels)==length(plotdata))
        {labelpos<-unique(round(pretty(c(0.5,1:length(plotdata)))))
        if(no_scenario_numbers) axis(1,labels=FALSE,at=labelpos[which(labelpos>0)],cex.axis=texcex,lwd=texcex,padj=0.5,tck=-0.05) else
          axis(1,at=labelpos[which(labelpos>0)],cex.axis=texcex,lwd=texcex,padj=0.5,tck=-0.05)
        axis(1,labels=FALSE,at=c(1:length(plotdata)),cex.axis=texcex,lwd=texcex,padj=0.5,tck=-0.05,tck=-0.025)
        legend[[scen]]<-data.frame(code=1:length(plotlabels),Label=plotlabels)} else
                    legend[[scen]]<-"no adequate labels provided" else legend[[scen]]<-"time series labels"
    
    
    
    if(scen==1)
     if(left_label)
       {axis(2,tck=0.02,lwd=texcex,las=1,cex.axis=texcex)
        axis(2,tck=-0.02,lwd=texcex,las=1,cex.axis=texcex,labels=FALSE)
        if(!is.na(metric_label)) mtext(metric_label,2,
                                   line=(par('omi')[2]/text_line_height)-1-
                                     strheight(metric_label,cex=texcex,font=2,units="inches")/text_line_height,
                                   cex=texcex,font=2,outer=TRUE)}
    if(scen==length(climate_scenario_list))
      if(right_label)
      {axis(4,tck=0.02,lwd=texcex,las=1,cex.axis=texcex)
        axis(4,tck=-0.02,lwd=texcex,las=1,cex.axis=texcex,labels=FALSE)
        if(!is.na(metric_label)) mtext(metric_label,4,line=(par('omi')[4]/text_line_height)-1.5,cex=texcex,font=2,
                                       outer=TRUE)}
    axis(2,labels=FALSE,tck=0.02,lwd=texcex)
    axis(4,labels=FALSE,tck=0.02,lwd=texcex)
    box(lwd=texcex)  
    
    if(!is.null(caption))
    {if(caption_above)
    {if(length(caption)>1) warning("too many captions provided - 1 is the max for captions above the plot")
      mtext(caption[1],3,1,font=2,cex=texcex)} else
      {if(length(caption)>0) mtext(caption[1],3,-1.1*texcex,font=2,cex=texcex)
        if(length(caption)>1) mtext(caption[2],3,-2.2*texcex,font=2,cex=texcex)
        if(length(caption)>2) mtext(caption[3],3,-3.3*texcex,font=2,cex=texcex)
        if(length(caption)>3) warning("too many captions provided - 3 is the max")}}
    
    
  }
  return(legend=legend)
}
    
    
    
