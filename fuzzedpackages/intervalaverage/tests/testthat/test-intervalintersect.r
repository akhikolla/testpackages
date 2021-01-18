
#full validation of intervalintersect via comparison with manual expansion of every time point on a random sample
test_that("intervalintersect",{

  set.seed(100)


##leverage the structure of the example from the vignette for some tests
  n_weeks <- 500
  n_loc <- 500
exposure_dataset3 <- rbindlist(lapply(1:n_loc, function(z){
  data.table(location_id=z, start=seq(as.IDate("2000-01-01"),by=7L,length=n_weeks),
             end=seq(as.IDate("2000-01-07"),by=7L,length=n_weeks),
             pm25=rnorm(n_weeks,mean=15),
             no2=rnorm(n_weeks,mean=25) )
} ))

n_ppt <- 50
addr_history <- data.table(ppt_id=paste0("ppt",1:n_ppt))
addr_history[, n_addr := rbinom(.N,size=length(unique(exposure_dataset3$location_id)),prob=.005)]

#repl=TRUE because it's possible for an address to be lived at multiple different time intervals:
addr_history <- addr_history[,
                             list(location_id=sample(exposure_dataset3$location_id,n_addr,replace=TRUE)),
                             by="ppt_id"
                             ]
addr_history

#note that not all of these n_loc locations in exposure_dataset3 were "lived at" in this cohort:
length(unique(addr_history$location_id))

#also note that it's possible for different participants to live at the same address.
addr_history[,list(loc_with_more_than_one_ppt=length(unique(ppt_id))>1),by=location_id][,sum(loc_with_more_than_one_ppt)]
#Because of the way I generated this data,
 #it's way more common than you'd expect in a real cohort
#but it does happen especially in cohorts
 #with familial recruitment or people living in nursing home complexes.

sample_dates <- function(n){
  stopifnot(n%%2==0)
  dateseq <- seq(as.IDate("1995-01-01"),as.IDate("2008-01-01"),by=1L)
  dates <- sort(sample(dateseq,n))
  #half of the time, make the last date "9999-01-01" which represents that the currently
  #lives at that location and we're carrying that assumption forward
  if(runif(1)>.5){
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

#FALSE is good--it means there's no overlap of dates within ppt.
test_that("example address history is non-overlapping",{
  expect_false(is.overlapping(addr_history,interval_vars=c("addr_start","addr_end"),
                              group_vars=c("ppt_id","location_id"))
  )
})


#examples are non-overlapping

  expect_false(is.overlapping(exposure_dataset3,interval_vars=c("start","end"),
                              group_vars="location_id")
  )

  #scramble the order and set different keys to check if state is returned after
  setkey(exposure_dataset3,no2)
  addr_history <- addr_history[sample(1:.N)]
  setkey(addr_history,addr_end)

  addr_history_original <- copy(addr_history)
  exposure_dataset3_original <- copy(exposure_dataset3)

z <- intervalintersect(x=exposure_dataset3,
                       y=addr_history,
                       interval_vars=c(
                         start="addr_start",
                         end="addr_end"
                       ),
                       group_vars=c("location_id"),
                       interval_vars_out=c("start2","end2")
)



##make sure state wasn't altered
expect_equal(exposure_dataset3_original, exposure_dataset3)
expect_equal(addr_history_original,addr_history)


##make sure that the order of x and y doesn't matter

z2 <- intervalintersect(x=addr_history,
                        y=exposure_dataset3,
                        interval_vars=c(
                          addr_start="start",
                          addr_end="end"
                        ),
                        group_vars=c("location_id"),
                        interval_vars_out=c("start2","end2")
)

setcolorder(z2, names(z))
setkey(z2,NULL)
setkeyv(z2,key(z))
expect_equal(z,z2)





###what happens if you have interval_vars with the same name?

z3 <- intervalintersect(x=
                          setnames(
                            copy(addr_history),
                            c("addr_start","addr_end"),c("start","end")
                            ),
                        y=exposure_dataset3,
                        interval_vars=c("start","end"),
                        group_vars=c("location_id"),
                        interval_vars_out=c("start2","end2")
)

expect_equal(z,z2)



###what happens if you have group_vars with different names?

z3 <- intervalintersect(x=
                          setnames(
                            copy(addr_history),
                            c("location_id"),c("addr_location_id")
                          ),
                        y=exposure_dataset3,
                        interval_vars=c(
                          addr_start="start",
                          addr_end="end"
                        ),
                        group_vars=c(addr_location_id="location_id"),
                        interval_vars_out=c("start2","end2")
)



#if the address history and the exposure datasets are both non-overlapping
#then the resulting intersect must also be non-overlapping

  expect_false(is.overlapping(z[,,],interval_vars=c("start2","end2"),
                              group_vars=c("location_id","ppt_id"))
  )






  #actually manually do the intersect by expansion:
  exposure_dataset3[,i:= 1:.N]
  exposure_dataset3_expanded <- exposure_dataset3[, list(date=seq(start,end,by=1),
                                                         location_id=rep(location_id,.N),
                                                         pm25=rep(pm25,.N),
                                                         no2=rep(no2,.N)),
                                                  by="i"]
  exposure_dataset3_expanded[,i:=NULL]

  addr_history[,i:=1:.N]


  addr_history_expanded <- addr_history[,list(date=seq(addr_start,addr_end,by=1),
                                              location_id=rep(location_id,.N),
                                              addr_id=rep(addr_id,.N),
                                              ppt_id=rep(ppt_id,.N)),
                                        by="i"]


  addr_history_expanded[,i:=NULL]


  setkey(addr_history_expanded,date,location_id)
  setkey(exposure_dataset3_expanded,date,location_id)

  p_expanded <- exposure_dataset3_expanded[addr_history_expanded,nomatch=NULL]
  z_expanded  <- z[,list(date=seq(start2,end2,by=1),
                         pm25=rep(pm25,.N),
                         no2=rep(no2,.N),
                         location_id=rep(location_id,.N)
  ) ,by=c("addr_id","ppt_id","start2","end2")]
  z_expanded[, start2:=NULL]
  z_expanded[, end2:=NULL]
  setcolorder(p_expanded, names(z_expanded))
  setkey(p_expanded,date, ppt_id, addr_id, location_id)
  setkey(z_expanded,date, ppt_id, addr_id, location_id)


  expect_equal(p_expanded,z_expanded)


})
