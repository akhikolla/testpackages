#' Plot daily climate metric accumulation throughout the year
#' 
#' This function generates a plot of the accumulation of a climate metric
#' throughout the year. Its standard output are the mean daily accumulation and
#' the standard deviation. It is also possible to add one or several so-called
#' focusyears to add the daily accumulation during these years to the plots.
#' Plots can be produced in R or directly exported as .png files.
#' 
#' Plots daily accumulation of climatic metrics, such as winter chill, as daily
#' accumulation rates or as cumulative accumulation. A legend is only added,
#' when focusyears are also shown. Otherwise the plot is reasonably
#' self-explanatory.
#' 
#' @param daily_chill a daily chill object generated with the daily_chill
#' function, which can calculate several standard chilling metrics or be
#' supplied with user-written temperature models. Since the format for the
#' input file must meet certain requirements, I recommend that you follow the
#' steps shown in the example below to prepare it.
#' @param metrics list of the metrics to be evaluated. This defaults to NA, in
#' which case the function makes a guess on what metrics you want to
#' calculated. This is done by choosing all column headers that are not
#' required for a daily_chill object.
#' @param startdate the first day of the season for which the metrics are to be
#' summarized (as a Julian date = day of the year)
#' @param enddate the last day of the season for which the metrics are to be
#' summarized (as a Julian date = day of the year)
#' @param useyears if only certain years are to be used, these can be provided
#' here as a numeric vector. Defaults to NA, which means all years in the
#' daily_chill object are used.
#' @param metriclabels Character vector with labels for each metric to be
#' analyzed. Defaults to NA, which means that the strings passed as metrics
#' will be used.
#' @param focusyears Numeric vector containing the years that are to be
#' highlighted in the plot. Years for which no data are available are
#' automatically removed.
#' @param cumulative Boolean argument (TRUE or FALSE) indicating whether the
#' climate metric should be shown as daily accumulation rates or as cumulative
#' accumulation.
#' @param image_type Character string indicating the file format that should be
#' output. Image files are only produced for the moment, if this is "png". All
#' other values, as well as the default NA lead to output as an R plot only.
#' @param outpath Path to the folder where the images should be saved. Should
#' include a trailing "/". The folder must already exists.
#' @param filename Suffix of the filenames for output graph files. These will
#' be amended by the name of the metric and by the file extension.
#' @param fonttype The type of font to be used for the figures. Can be 'serif'
#' (default) for a Times New Roman like font, 'sans' for an Arial type font or
#' 'mono' for a typewriter type font.
#' @param title title of the plot (if unhappy with the default).
#' @param plotylim numeric vector of length 2 indicating the extent of the y axis.
#' Defaults to NA, which means that y limits are determined automatically.
#' @return The main purpose of the function is a side effect - plots of daily
#' climate metric accumulation. However, all the data used for making the plots
#' is returned as a list containing an element for each metric, which consists
#' of a data.table with the daily means, standard deviation and daily values
#' for all focusyears.
#' @author Eike Luedeling
#' @keywords utility
#' @examples
#' 
#' day_chill<-make_daily_chill_plot(daily_chill(stack_hourly_temps(fix_weather(
#'   KA_weather[which(KA_weather$Year>2005),])),
#'   running_mean=11),focusyears=c(2001,2005),cumulative=TRUE,startdate=300,enddate=30)
#' 
#'  
#' @export make_daily_chill_plot
make_daily_chill_plot <-function(daily_chill,metrics=NA,startdate=1,enddate=366,useyears=NA,metriclabels=NA,
                                 focusyears="none",cumulative=FALSE,image_type=NA,outpath=NA,filename=NA,fonttype='serif',
                                 title=NA,plotylim=NA)
{

 make_plot<-function(plottable,name,ylabel,focusyears,imageout,fonttype="default")
   {if (imageout) {cex.main<-4;linelwd=3} else {cex.main<-1;linelwd=2}

    JDay<-plottable$JDay
    means=plottable$Mean
    sds=plottable$Sd
    months<-c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
    leg<-strptime(strptime(paste("01/",1,"/2001",sep=""),format="%d/%m/%Y")-86400+JDay*86400,"%Y-%m-%d")
    tick_marks<-leg[sapply(strsplit(as.character(leg),"-"),"[",3)=="01"]
    tick_label_pos<-leg[sapply(strsplit(as.character(leg),"-"),"[",3)=="15"]
    tick_labels<-as.Date(tick_label_pos,origin="1999-1-2")
    tick_labels<-as.POSIXlt(tick_labels)$mon
    tick_labels<-months[tick_labels+1]
    tick_marks<-as.numeric(format(strptime(tick_marks,"%Y-%m-%d"),format="%j"))
    if(enddate<=startdate) tick_marks[which(tick_marks>=startdate)]<-tick_marks[which(tick_marks>=startdate)]-366
    ticks_labels<-as.numeric(format(strptime(tick_label_pos,"%Y-%m-%d"),format="%j"))
    if(enddate<=startdate) ticks_labels[which(ticks_labels>=startdate)]<-ticks_labels[which(ticks_labels>=startdate)]-366

    if(!(focusyears[1]=="none")) layout(matrix(c(1,2),1,2,byrow=TRUE),widths=c(3,1))   else par(mfcol=c(1,1))

    if(imageout) par(mar=c(5.1,7.1,4.1,0.5))  else par(mar=c(0,0,1.5,0))
    suppressWarnings(if(!(focusdata[1]=="none")) {ymin<-min(c(means-sds,min(plottable[,as.character(focusyears)],na.rm=TRUE)),na.rm=TRUE)
       ymax<-max(c(means+sds,max(plottable[,as.character(focusyears)],na.rm=TRUE)),na.rm=TRUE)} else
      {ymin<-min(means-sds,na.rm=TRUE)
       ymax<-max(means+sds,na.rm=TRUE)})
    if(!is.na(plotylim[1]))
      {ymin<-plotylim[1]
       ymax<-plotylim[2]}
    ymin<-ymin-(ymax-ymin)/20
    ymax<-ymax+(ymax-ymin)/20
    if(is.na(title)) main<-paste(name,"accumulation") else main<-title
    if(length(which(is.na(means)))==0)
      plot(JDay,means,main=main,cex.main=cex.main,ylab=NA,xlab=NA,xaxs="i",yaxs="i",xaxt="n",yaxt="n",type="l",col="BLACK",ylim=c(ymin,ymax)) else
    plot(JDay,rep(NA,length(JDay)),main=paste(name,"accumulation"),cex.main=cex.main,ylab=NA,xlab=NA,xaxs="i",yaxs="i",xaxt="n",yaxt="n",type="l",col="BLACK",ylim=c(0,10))
    if(length(which(is.na(means)))==0)  arrows(JDay,means+sds,JDay,means-sds, angle=90, code=3,lwd=linelwd,length=0,col="GRAY")
    if(length(which(is.na(means)))==0)  lines(JDay,means,lwd=linelwd)
    if(length(which(is.na(means)))==0)  suppressWarnings(if(!(focusyears[1]=="none"))
    for (focusyear in focusyears)  lines(JDay,plottable[,as.character(focusyear)],lwd=linelwd,col=rainbow(length(focusyears))[which(focusyear==focusyears)]))

    if (!imageout) {
      axis(1,labels=FALSE,at=tick_marks,padj=1,lwd.ticks=3,lwd=2)
      axis(2,lwd=2,lwd.ticks=2)
      axis(1,lwd.ticks=0,at=ticks_labels,labels=tick_labels,lwd=2)
      mtext(side=2,text=ylabel,line=3,font=2)
      box(which="plot",lwd=2)}

    if(imageout) {
      axis(1,labels=FALSE,at=tick_marks,padj=1,lwd.ticks=3,lwd=3)
      axis(2,lwd=3,lwd.ticks=3,cex.axis=3,padj=0)
      axis(1,lwd.ticks=0,at=ticks_labels,labels=tick_labels,lwd=3,cex.axis=3,padj=1)
      mtext(side=2,text=ylabel,line=4,cex=4,font=2)
      box(which="plot",lwd=3)}



    if(length(which(is.na(means)))>0)  text(x=mean(JDay),y=5,"no data to show")
    if(!(focusyears[1]=="none"))
    {
      if(!imageout)
        {par(mar=c(3.1,0,2.1,0))
        yend<-min(-(length(focusyears)+3),-15)
        plot(-100,ylim=c(yend,0),xlim=c(0,5),axes=FALSE,xlab="")
        for (focusyear in focusyears)  lines(c(0.5,2),y=rep(-3-which(focusyear==focusyears),2),
            lwd=2,col=rainbow(length(focusyears))[which(focusyear==focusyears)])
        text(y=-3-1:length(focusyears),x=2.5,labels=focusyears,adj=0)
        lines(c(0.5,2),y=c(-1,-1),lwd=2,col="BLACK")
        text(y=-1,x=2.5,labels="Mean",adj=0)
        lines(c(0.5,2),y=c(-2,-2),lwd=2,col="GRAY")
         text(y=-2,x=2.5,labels="Std. dev.",adj=0)}

      if(imageout)
      {par(mar=c(5.1,0,4.1,0))
        yend<-min(-(length(focusyears)+3),-15)
        plot(-100,ylim=c(yend,0),xlim=c(0,5),axes=FALSE,xlab="")
        for (focusyear in focusyears)  lines(c(0.5,2),y=rep(-3-which(focusyear==focusyears),2),
                                             lwd=4,col=rainbow(length(focusyears))[which(focusyear==focusyears)])
        text(y=-3-1:length(focusyears),x=2.5,labels=focusyears,adj=0,cex=3)
        lines(c(0.5,2),y=c(-1,-1),lwd=4,col="BLACK")
        text(y=-1,x=2.5,labels="Mean",adj=0,cex=3)
        lines(c(0.5,2),y=c(-2,-2),lwd=4,col="GRAY")
        if(fonttype=="mono")  text(y=-2,x=2.5,labels="St.dev.",adj=0,cex=3) else
          text(y=-2,x=2.5,labels="Std. dev.",adj=0,cex=3)}

    }
 }





dc<-daily_chill$daily_chill
dc[,"JDay"]<-strptime(paste(dc$Month,"/",dc$Day,"/",dc$Year,sep=""),"%m/%d/%Y")$yday+1
if("no_Tmin" %in% colnames(dc) & "no_Tmax" %in% colnames(dc))
dc<-dc[1:max(which(!(dc$no_Tmin|dc$no_Tmax))),]
if(enddate>startdate) relevant_days<-c(startdate:enddate) else relevant_days<-c(startdate:366,1:enddate)
relevant_days<-relevant_days[which(relevant_days %in% unique(dc$JDay))]

if(enddate>startdate) dc[which(dc$JDay %in% relevant_days),"End_year"]<-dc[which(dc$JDay %in% relevant_days),"Year"] else
{dc[which(dc$JDay %in% c(1:enddate)),"End_year"]<-dc[which(dc$JDay %in% c(1:enddate)),"Year"]
dc[which(dc$JDay %in% c(startdate:366)),"End_year"]<-dc[which(dc$JDay %in% c(startdate:366)),"Year"]+1}
if(is.na(useyears[1])) useyears<-unique(dc$End_year)[which(!is.na(unique(dc$End_year)))]
dc<-dc[which(dc$End_year %in% useyears),]
if(!focusyears[1]=="none") focusyears<-focusyears[which(focusyears %in% dc$End_year)]
if(length(focusyears)==0) focusyears<-"none"

if(is.na(metrics[1])) metrics<-colnames(dc) [which(!colnames(dc) %in% c("YYMMDD","Year","Month","Day","Tmean","JDay","End_year","no_Tmin","no_Tmax"))]
if(is.na(metriclabels[1])) metriclabels<-metrics

df<-list()
focusdata=list()

if(startdate<enddate) Jdays<-relevant_days else
  Jdays<-c((c(startdate:366)-366),c(1:enddate))

for (met in metrics)
{

  if(cumulative==TRUE)
        for(e in unique(dc$End_year))
         dc[which(dc$End_year==e&(dc$JDay %in% relevant_days)),met]<-cumsum(dc[which(dc$End_year==e&(dc$JDay %in% relevant_days)),met])
  dc[which(!dc$JDay %in% relevant_days)]<-0

  df[[met]]<-data.frame(JDay=NA,Mean=NA,Sd=NA)

  if(cumulative==TRUE)
  {
  complete<-data.frame(End_year=unique(dc$End_year),ndays=NA)
  for(e in unique(dc$End_year))
   complete[which(complete$End_year==e),"ndays"]<-length(which(dc$End_year==e&dc$JDay %in% relevant_days))

  exclude_from_mean<-complete[which(complete$ndays<max(complete$ndays)*0.9),"End_year"]
  } else exclude_from_mean=NULL

  for (i in 1:length(relevant_days)) {df[[met]][i, ] <- c(Jdays[i],mean(dc[which((!dc$End_year %in% exclude_from_mean)&(dc$JDay == relevant_days[i])),met],na.rm=TRUE),
                                                          sd(dc[which((!dc$End_year %in% exclude_from_mean)&(dc$JDay == relevant_days[i])),met],na.rm=TRUE))
  df[[met]][which(is.na(df[[met]][,"Sd"]))  ,"Sd" ]<-0
  }


  if(!is.null(focusyears)) if(!(focusyears[1]=="none")) {
  for(focusyear in focusyears)
  {
  focusdata<-dc[which(dc$End_year==focusyear),]
  if(startdate<enddate) focusdata[,"Jdays"]<-focusdata[,"JDay"] else
  {focusdata[which(focusdata$JDay<startdate),"Jdays"]<-focusdata$JDay[which(focusdata$JDay<startdate)]
  focusdata$Jdays[which(focusdata$JDay>=startdate)]<-focusdata$JDay[which(focusdata$JDay>=startdate)]-366
  }
  df[[met]][which(df[[met]]$JDay %in% focusdata$Jdays),as.character(focusyear)]<-focusdata[which(focusdata$Jdays %in% df[[met]]$JDay),met]
  }} else focusdata[[met]]<-"none" else focusdata[[met]]<-"none"
}


for(met in metrics)
{metlab<-metriclabels[which(met==metrics)]
if(!is.na(image_type)) if(image_type=="png") {suppressWarnings(dir.create(outpath))
  png(file.path(outpath,paste(filename,"_",met,".png",sep="")),width=1250,height=1000)
  imageout<-TRUE} else imageout<-FALSE else imageout<-FALSE
par(family=fonttype)
if (!cumulative)  make_plot(plottable=df[[met]],name=metlab,ylabel=paste(metlab,"per day"),focusyears,imageout,fonttype)
if (cumulative)   make_plot(plottable=df[[met]],name=metlab,ylabel=paste("Cumulative", metlab),focusyears,imageout,fonttype)
if(!imageout) if(length(metrics)>1&!which(metrics==met)==length(metrics)) readline()
if(!is.na(image_type)) if(image_type=="png") dev.off()
}

  #else "no data to show"
return(df)
}





#' Plot daily climate metric accumulation throughout the year (2)
#' 
#' This function generates a plot of the accumulation of a climate metric
#' throughout the year. Its standard output are the mean daily accumulation and
#' the standard deviation. It is also possible to add one or several so-called
#' focusyears to add the daily accumulation during these years to the plots.
#' Plots can be produced in R or directly exported as .png files.
#' 
#' Plots daily accumulation of climatic metrics, such as winter chill, as daily
#' accumulation rates or as cumulative accumulation. A legend is only added,
#' when focusyears are also shown. Otherwise the plot is reasonably
#' self-explanatory.
#' 
#' @param daily an object generated with the daily_chill
#' function, which can calculate several standard chilling metrics or be
#' supplied with user-written temperature models. Since the format for the
#' input file must meet certain requirements, I recommend that you follow the
#' steps shown in the example below to prepare it.
#' @param metrics list of the metrics to be evaluated. This defaults to NA, in
#' which case the function makes a guess on what metrics you want to
#' calculated. This is done by choosing all column headers that are not
#' required for a daily_chill object.
#' @param startdate the first day of the season for which the metrics are to be
#' summarized (as a Julian date = day of the year)
#' @param enddate the last day of the season for which the metrics are to be
#' summarized (as a Julian date = day of the year)
#' @param useyears if only certain years are to be used, these can be provided
#' here as a numeric vector. Defaults to NA, which means all years in the
#' daily_chill object are used.
#' @param metriclabels Character vector with labels for each metric to be
#' analyzed. Defaults to NA, which means that the strings passed as metrics
#' will be used.
#' @param focusyears Numeric vector containing the years that are to be
#' highlighted in the plot. Years for which no data are available are
#' automatically removed.
#' @param cumulative Boolean argument (TRUE or FALSE) indicating whether the
#' climate metric should be shown as daily accumulation rates or as cumulative
#' accumulation.
#' @param fix_leap boolean parameter indicating whether the anomaly that can
#' originate when leaf years are present in the data should be smoothed by
#' interpolating between Dec 30 and Jan 1 in leap years.
#' @return The main purpose of the function is a side effect - plots of daily
#' climate metric accumulation. However, all the data used for making the plots
#' is returned as a list containing an element for each metric, which consists
#' of a data.table with the daily means, standard deviation and daily values
#' for all focusyears.
#' @import ggplot2
#' @author Eike Luedeling
#' @keywords utility
#' @examples
#' 
#' 
#' daily<-daily_chill(stack_hourly_temps(fix_weather(
#'   KA_weather[which(KA_weather$Year>2005),])),running_mean=11)
#' 
#' make_daily_chill_plot2(daily,metrics=c("Chill_Portions","GDH"),cumulative=TRUE,
#'    startdate=300,enddate=30,focusyears=c(2009,2008))
#' 
#' 
#'  
#' @export make_daily_chill_plot2
make_daily_chill_plot2 <-function(daily,metrics=NA,startdate=1,enddate=366,useyears=NA,metriclabels=NA,
                                 focusyears="none",cumulative=FALSE,fix_leap=TRUE)
{
  if(length(metriclabels)==length(metrics))
    metriclabels[which(is.na(metriclabels))]<-metrics[which(is.na(metriclabels))] else
      metriclabels<-metrics
  
  dc<-daily$daily_chill
  dc[,"JDay"]<-strptime(paste(dc$Month,"/",dc$Day,"/",dc$Year,sep=""),"%m/%d/%Y")$yday+1
  if("no_Tmin" %in% colnames(dc) & "no_Tmax" %in% colnames(dc))
    dc<-dc[1:max(which(!(dc$no_Tmin|dc$no_Tmax))),]
  if(enddate>startdate) relevant_days<-c(startdate:enddate) else relevant_days<-c(startdate:366,1:enddate)
  relevant_days<-relevant_days[which(relevant_days %in% unique(dc$JDay))]
  if(is.na(metrics[1])) metrics<-colnames(dc) [which(!colnames(dc) %in% c("YYMMDD","Year","Month","Day","Tmean","JDay","End_year","no_Tmin","no_Tmax"))]
  if(is.na(metriclabels[1])) metriclabels<-metrics
  

  if(startdate<enddate) Jdays<-relevant_days else
    Jdays<-c((c(startdate:366)-366),c(1:enddate))
  
  dc<-make_JDay(dc)
  dc[,"seasonday"]<-unlist(sapply(dc$JDay,function(x) if(x %in% relevant_days) Jdays[which(relevant_days==x)] else NA))
  dc[,"End_year"]<-dc$Year+(dc$seasonday<=0)
  
  if(cumulative)
   for (met in metrics)
    {for(e in unique(dc$End_year))
      dc[which(dc$End_year==e&(dc$JDay %in% relevant_days)),met]<-cumsum(dc[which(dc$End_year==e&(dc$JDay %in% relevant_days)),met])
    dc[which(is.na(dc$End_year)),met]<-NA}
  
  
  endyears<-as.numeric(names(table(dc$End_year))[table(dc$End_year) %in% c(max(table(dc$End_year)),max(table(dc$End_year)-1))])
  dc<-dc[which(dc$End_year %in% endyears),]
  
  output<-list()
  
  for (met in metrics)
    {dc_summ <- plyr::ddply(dc[which(!is.na(dc[,met])&(dc$End_year %in% endyears)),], "seasonday", .drop=TRUE,
                   .fun = function(xx, col) {
                            c(N    = sum(!is.na(xx[[col]])),
                              mean = mean   (xx[[col]], na.rm=TRUE),
                              sd   = sd     (xx[[col]], na.rm=TRUE))},
                   met)
     dc_summ$sd[which(is.na(dc_summ$sd))]<-0
     dc_summ[,"se"]<-dc_summ$sd/sqrt(dc_summ$N)
     if(fix_leap) dc_summ$mean[which(dc_summ$seasonday==0)]<-mean(c(dc_summ$mean[which(dc_summ$seasonday %in% c(-1,1))]))
     
     focusyears<-focusyears[which(focusyears %in% endyears)]
     
     if(!length(focusyears)==0)
       {
         for(foc in focusyears)
           {dc_summ<-merge(dc_summ,dc[which(dc$End_year==foc),c("End_year","seasonday",met)])
            dc_summ<-dc_summ[,c(1:(ncol(dc_summ)-2),ncol(dc_summ))]
            colnames(dc_summ)[ncol(dc_summ)]<-foc}
       foci<-reshape2::melt(dc_summ[,c("seasonday",focusyears)],id.vars="seasonday")
       colnames(foci)<-c("seasonday","Year","value")
       foci[,"Date"]<-as.Date(foci$seasonday,origin=as.Date("2011-12-31"))
       }
       
     dc_summ[,"Date"]<-as.Date(dc_summ$seasonday,origin=as.Date("2011-12-31"))

     
     all_dates<-as.Date(dc_summ$seasonday,origin=as.Date("2011-12-31"))
     months<-as.numeric(format(all_dates, "%m"))
     days<-as.numeric(format(all_dates, "%d"))
     griddatelist<-all_dates[which((months<3&days==1)|(months>2&days<3))]
     days<-as.numeric(format(griddatelist, "%d"))
     gridlabels<-format(griddatelist, "%b")
     gridlabels[which(days==2)]<-""
     
     if(cumulative) ylabe<-paste("Cumulative",metriclabels[which(metrics==met)]) else ylabe<-metriclabels[which(metrics==met)]

     
     
     if(!length(focusyears)==0)
       print(ggplot2::ggplot(dc_summ,aes_(x=~Date,y=~mean)) +
         geom_errorbar(aes_(ymin=~mean-sd, ymax=~mean+sd), width=.1,lwd=2,colour="grey") +
         geom_line(colour="black",lwd=1) +
         geom_line(data=foci,aes_(x=~Date,y=~value,colour=~Year),lwd=2) +
         ylab(ylabe) + xlab("Date") +
         theme_gray(base_size=22) +
         scale_x_date(breaks=griddatelist,labels=gridlabels,expand=expand_scale(add=0),minor_breaks=NULL)
       )
     
     if(length(focusyears)==0)
       print(ggplot2::ggplot(dc_summ,aes_(x=~Date,y=~mean)) +
               geom_errorbar(aes_(ymin=~mean-sd, ymax=~mean+sd), width=.1,lwd=2,colour="grey") +
               geom_line(colour="black",lwd=1) +
               ylab(ylabe) + xlab("Date") +
               theme_gray(base_size=22) +
               scale_x_date(breaks=griddatelist,labels=gridlabels,expand=expand_scale(add=0),minor_breaks=NULL) 
       )

          output[[met]]<-dc_summ[,1:ncol(dc_summ)-1]
     
     if(!met==metrics[length(metrics)]) readline(prompt="Press [enter] to continue")
     }
     
  return(output)
}


