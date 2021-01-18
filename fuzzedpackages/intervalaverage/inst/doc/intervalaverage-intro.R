## ---- include = FALSE----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup,echo=FALSE,message=FALSE--------------------------------------
library(intervalaverage)
set.seed(1)

## ------------------------------------------------------------------------
exposure_dataset <- data.table(location_id=1, start=seq(as.IDate("2000-01-01"),by=7,length=4),
                end=seq(as.IDate("2000-01-07"),by=7,length=4),pm25=rnorm(4,mean=15), no2=rnorm(4,mean=25))
exposure_dataset

## ------------------------------------------------------------------------
exposure_dataset[start %in% as.IDate(c("2000-01-01","2000-01-08")),mean(pm25)]

## ------------------------------------------------------------------------
exposure_dataset[start %in% as.IDate(c("2000-01-01","2000-01-08")),weighted.mean(pm25,w=c(7/10,3/10))]

## ------------------------------------------------------------------------
averaging_periods <- data.table(start=seq(as.IDate("2000-01-01"),by=10,length=3),
                end=seq(as.IDate("2000-01-10"),by=10,length=3))
averaging_periods

## ------------------------------------------------------------------------
averaged_exposures <- intervalaverage(x=exposure_dataset,y=averaging_periods,
                                      interval_vars=c("start","end"),value_vars=c("pm25","no2")
                                      )
averaged_exposures[, list(start,end,pm25,no2)]

## ------------------------------------------------------------------------
averaged_exposures

## ------------------------------------------------------------------------
intervalaverage(x=exposure_dataset,
                y=averaging_periods,
                interval_vars=c("start","end"),
                value_vars=c("pm25","no2"),
                required_percentage = 75)

## ------------------------------------------------------------------------
exposure_dataset2 <- rbindlist(lapply(1:3, function(z){
  data.table(location_id=z, 
             start=seq(as.IDate("2000-01-01"),by=7,length=4),
             end=seq(as.IDate("2000-01-07"),by=7,length=4),pm25=rnorm(4,mean=15),
             no2=rnorm(4,mean=25))} ))
exposure_dataset2

## ------------------------------------------------------------------------
#unexpanded:
averaging_periods
#expanded to every location_id:
rbindlist(lapply(1:3, function(z)copy(averaging_periods)[,location_id:=z][])) 

## ------------------------------------------------------------------------
exposure_dataset2_unique_locs <- data.table(location_id=unique(exposure_dataset2$location_id))
averaging_periods2 <- CJ.dt(averaging_periods, exposure_dataset2_unique_locs)
averaging_periods2

#or, more concisely:
averaging_periods2 <- CJ.dt(averaging_periods, unique(exposure_dataset2[,list(location_id)]))

## ------------------------------------------------------------------------
intervalaverage(x=exposure_dataset2,
                y=averaging_periods2,
                interval_vars=c("start","end"),
                value_vars=c("pm25","no2"),
                group_vars="location_id",
                required_percentage = 75)[, list(location_id, start,end, pm25,no2)]

## ------------------------------------------------------------------------
exposure_dataset_overlapping <- rbindlist(lapply(1:3, function(z){
  data.table(location_id=z, 
             start=seq(as.IDate("2000-01-01"),by=7,length=4),
             end=seq(as.IDate("2000-01-08"),by=7,length=4),
             pm25=rnorm(4,mean=15),
             no2=rnorm(4,mean=25) )
} ))
exposure_dataset_overlapping

## ----error=TRUE----------------------------------------------------------
intervalaverage(exposure_dataset_overlapping,averaging_periods2,
                interval_vars=c("start","end"),
                value_vars=c("pm25","no2"),
                group_vars="location_id",
                required_percentage = 75)

## ------------------------------------------------------------------------
is.overlapping(exposure_dataset_overlapping, interval_vars=c('start','end'),group_vars="location_id")

## ------------------------------------------------------------------------
exposure_dataset_isolated <- isolateoverlaps(exposure_dataset_overlapping, 
                                             interval_vars=c("start","end"), 
                                             group_vars="location_id",
                                             interval_vars_out=c("start2","end2"))
exposure_dataset_isolated[1:15] #only show the first 15 rows

## ------------------------------------------------------------------------
exposure_dataset_overlaps_averaged <-
  exposure_dataset_isolated[, list(pm25=mean(pm25),no2=mean(no2)),by=c("location_id","start2","end2")]

setnames(exposure_dataset_overlaps_averaged, c("start2","end2"),c("start","end"))
exposure_dataset_overlaps_averaged

## ------------------------------------------------------------------------
intervalaverage(exposure_dataset_overlaps_averaged,
                averaging_periods2,
                interval_vars=c("start","end"),
                value_vars=c("pm25","no2"),
                group_vars="location_id",
                required_percentage = 75)[,list(location_id, start,end,pm25,no2)]

## ------------------------------------------------------------------------

overlapping_averaging_periods <- data.table(start=as.IDate(c("2000-01-01","2000-01-01")),
                                            end=as.IDate(c("2000-01-10","2000-01-20"))
                                            )
overlapping_averaging_periods
overlapping_averaging_periods_expanded <-
  CJ.dt(overlapping_averaging_periods,unique(exposure_dataset2[,list(location_id)]))

overlapping_averaging_periods_expanded


intervalaverage(exposure_dataset_overlaps_averaged,
                overlapping_averaging_periods_expanded,
                interval_vars=c("start","end"),
                value_vars=c("pm25","no2"),
                group_vars="location_id",
                required_percentage = 75)[,list(location_id, start,end,pm25,no2)]

## ------------------------------------------------------------------------
n_locs <- 2000
n_weeks <- 1000
exposure_dataset3 <- rbindlist(
  lapply(1:n_locs, function(id) {
    data.table(
      location_id = id,
      start = seq(as.IDate("2000-01-01"), by = 7, length = n_weeks),
      end = seq(as.IDate("2000-01-07"), by = 7, length = n_weeks),
      pm25 = rnorm(n_weeks, mean = 15),
      no2 = rnorm(n_weeks, mean = 25)
    )
  })
)

exposure_dataset3

## ------------------------------------------------------------------------
averaging_periods3 <- data.table(location_id=1:n_locs,
                                 start=sample(
                                   x=seq(as.IDate("2000-01-01"),as.IDate("2019-12-31"),by=1),
                                   size=n_locs
                                 )
)
averaging_periods3[,end:=start+round(3*365.25)]
averaging_periods3

## ------------------------------------------------------------------------
intervalaverage(exposure_dataset3,
                averaging_periods3,
                interval_vars=c("start","end"),
                value_vars=c("pm25","no2"),
                group_vars="location_id")[,list(location_id,start,end,pm25,no2)]

## ------------------------------------------------------------------------
averaging_periods3[, avg3yr:=end-round(3*365.25)]
averaging_periods3[, avg2yr:=end-round(2*365.25)]
averaging_periods3[, avg1yr:=end-round(1*365.25)]
#reshape the data.table:
averaging_periods4 <- melt(averaging_periods3,id.vars=c("location_id","end"),
                           measure.vars = c("avg3yr","avg2yr","avg1yr")) 
setnames(averaging_periods4, "value","start")
setnames(averaging_periods4, "variable","averaging_period")
averaging_periods4

## ------------------------------------------------------------------------

intervalaverage(exposure_dataset3,averaging_periods4,interval_vars=c("start","end"),
                value_vars=c("pm25","no2"),
                group_vars=c("location_id"),
                required_percentage = 75)[,list(location_id,start,end,pm25,no2)]

## ----echo=FALSE----------------------------------------------------------
address_history0 <- data.table(addr_id=c(1,2,2,3,5),
                 ppt_id=c(1,1,1,2,2),
                 addr_start=c(1L,10L,12L,1L,13L),
                 addr_end=c(9L,11L,14L,12L,15L)) 
exposure_dataset5 <- data.table(addr_id=rep(1:4,each=3),
  exp_start=rep(c(1L,8L,15L),times=4),
  exp_end=rep(c(7L,14L,21L),times=4),
  exp_value=c(rnorm(12))
 )

## ------------------------------------------------------------------------
exposure_dataset5

## ------------------------------------------------------------------------
address_history0

## ------------------------------------------------------------------------

exposure_addresss_table <- intervalintersect(exposure_dataset5,
                                             address_history0,
                                             interval_vars=c(exp_start="addr_start",
                                                             exp_end="addr_end"),
                                             "addr_id")
exposure_addresss_table

## ------------------------------------------------------------------------
setdiff(address_history0$addr_id,exposure_addresss_table$addr_id)
setdiff(exposure_dataset5$addr_id,exposure_addresss_table$addr_id)

## ------------------------------------------------------------------------
n_ppt <- 300
addr_history <- data.table(ppt_id=paste0("ppt",sprintf("%03d", 1:n_ppt)))
addr_history[, n_addr := rbinom(.N,size=length(unique(exposure_dataset3$location_id)),prob=.001)]
addr_history[n_addr <1L, n_addr := 1L] 

#repl=TRUE because it's possible for an address to be lived at multiple different time intervals:
addr_history <- addr_history[, 
                             list(location_id=sample(exposure_dataset3$location_id,n_addr,replace=TRUE)),
                             by="ppt_id"
                             ]   
addr_history

#note that not all of these 2000 locations in exposure_dataset3 were "lived at" in this cohort:
length(unique(addr_history$location_id)) 

#also note that it's possible for different participants to live at the same address.
addr_history[,list(loc_with_more_than_one_ppt=length(unique(ppt_id))>1),
             by=location_id][,
                             sum(loc_with_more_than_one_ppt)]
#Because of the way I generated this data, it's way more common than you'd expect in a real cohort 
#but it does happen especially in cohorts with familial recruitment or people living in nursing home complexes.

## ----echo=FALSE----------------------------------------------------------
#generate a vector from which dates will be sampled

sample_dates <- function(n){
  stopifnot(n%%2==0)
  dateseq <- seq(as.IDate("1960-01-01"),as.IDate("2015-01-01"),by=1)
  dates <- sort(sample(dateseq,n))
  #90% of the time, make the last date "9999-01-01" which represents that the currently
   #lives at that location and we're carrying that assumption forward
  if(runif(1)>.1){
    dates[length(dates)] <- as.IDate("9999-01-01") 
  }
  dates
}


addr_history_dates <- addr_history[,list(date=sample_dates(.N*2)) ,by="ppt_id"] #for every address, ppt needs two dates: start and end
addr_history_dates_wide <- addr_history_dates[, list(start=date[(1:.N)%%2==1],end=date[(1:.N)%%2==0]),by="ppt_id"]
addr_history <- cbind(addr_history,addr_history_dates_wide[, list(start,end)])
setnames(addr_history, c("start","end"),c("addr_start","addr_end"))
#addr_history[,any(end=="9999-01-01"),by="ppt_id"][, sum(V1)]
setkey(addr_history, ppt_id, addr_start)
#here i'm using a trick to map distinct values of location_id to integers (within ppt) by coercing to factor then back to numeric.
addr_history[, addr_id:=paste0(ppt_id, "_", as.numeric(as.factor(location_id))),by=ppt_id]

## ------------------------------------------------------------------------
addr_history

## ------------------------------------------------------------------------
is.overlapping(addr_history,interval_vars=c("addr_start","addr_end"),
               group_vars="ppt_id") #FALSE is good--it means there's no overlap of dates within ppt.

## ------------------------------------------------------------------------
exposure_dataset3_addr <- unique(addr_history[, list(location_id,addr_id)])[exposure_dataset3, 
                                                                            on=c("location_id"),
                                                                            allow.cartesian=TRUE,
                                                                            nomatch=NULL]
exposure_dataset3
exposure_dataset3_addr

## ------------------------------------------------------------------------
exposure_dataset3[, sum(duplicated(location_id)),by=c("start")][,max(V1)] #no duplicate locations at any date
exposure_dataset3_addr[, sum(duplicated(location_id)),by=c("start")][,max(V1)] 

## ------------------------------------------------------------------------
exposure_dataset3_addr[, location_id:=NULL]

## ------------------------------------------------------------------------
z <- intervalintersect(x=exposure_dataset3,
               y=addr_history,
               interval_vars=c(
                 start="addr_start",
                 end="addr_end"
                 ),
               group_vars=c("location_id"),
               interval_vars_out=c("start2","end2")
) 

z_addr <- intervalintersect(x=exposure_dataset3_addr,
               y=addr_history,
               interval_vars=c(
                 start="addr_start",
                 end="addr_end"
                 ),
               group_vars=c("addr_id"),
               interval_vars_out=c("start2","end2")
) 

setkey(z,ppt_id,start2,end2)
setkey(z_addr,ppt_id,start2,end2)
all.equal(z,z_addr)

z

## ------------------------------------------------------------------------

final_averaging_periods <- data.table(ppt_id=sort(unique(addr_history$ppt_id)))
final_averaging_periods[, end2:=sample(seq(as.IDate("2003-01-01"),as.IDate("2015-01-01"),by=1),.N)]
final_averaging_periods[,start2:=as.IDate(floor(as.numeric(end2-3*365.25)))]
final_averaging_periods

intervalaverage(z,final_averaging_periods, interval_vars=c("start2","end2"),
                value_vars=c("pm25","no2"),group_vars="ppt_id",required_percentage = 95
                )

