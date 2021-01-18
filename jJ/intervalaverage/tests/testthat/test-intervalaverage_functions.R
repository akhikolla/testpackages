



#####test isolateoverlaps ####
test_that("isolateoverlaps", {

  #isolateoverlaps simple example
  x <- data.table(
    start0 = c(1L, 5L, 5L),
    end0 = c(5L, 5L, 10L),
    id1 = "1",
    id2 = "1"
  )
  expect_equal(
    isolateoverlaps(
      x,
      interval_vars = c("start0", "end0"),
      group_vars = c("id1", "id2")
    ),
    structure(
      list(
        id1 = c("1", "1", "1", "1", "1"),
        id2 = c("1",
                "1", "1", "1", "1"),
        start = c(1, 5, 5, 5, 6),
        end = c(4, 5,
                5, 5, 10),
        start0 = c(1, 1, 5, 5, 5),
        end0 = c(5, 5, 5, 10, 10)
      ),
      row.names = c(NA,-5L),
      class = c("data.table", "data.frame")
    )
  )

  set.seed(90)
  n <- 1000
  x <- matrix(as.integer(round(runif(n=n*2L, 0L,1000L))),ncol=2)
  x <- as.data.table(t(apply(x,1,sort)))
  setnames(x,names(x),c("start0","end0"))
  x[, id1:=rbinom(n,3,prob=.3)]
  x[, id2:=rbinom(n,7,prob=.5)]
  x[, rnow:=1:.N]
  xd <- isolateoverlaps(x,interval_vars=c("start0","end0"),group_vars=c("id1","id2"))

  ##check that the union of the new intervals
   #span exactly the original intervals:
  expect_true(xd[,{
    x <- .SD[,list(v1=seq(start,end)),by=1:nrow(.SD)]$v1
    x <- unique(sort(x))
    list(V1= min(x)==start0 &
    max(x)==end0 &
    all(diff(x)==1) &
    sum(duplicated(.SD))==0) #within original intervals and grouping variables, no duplicated new intervals
  },by=c("id1","id2", "start0","end0"),.SDcols=c("start","end")][,all(V1)])

  expect_true(xd[,list(V1=all(start0[1]==start0&end0[1]==end0)),by=rnow][,all(V1)])
  expect_true(xd[,list(V1=identical(seq(min(start),max(end)),seq(start0[1],end0[1]))),by=rnow][,all(V1)])

  #grouping by rows in x, take every row in x and expand the start0 to end0, combine all those
  #expanded values into a single vector, sort that,
   #and compare that to the original expanded inteval
  system.time(expect_true(xd[,list(V1=identical(
                         sort(unlist(mapply(SIMPLIFY=FALSE,start,end,FUN=function(s,e){
                             seq(s,e)
                         }))),
                         seq(start0[1],end0[1]))
  ),
  by=rnow][,all(V1)]))


  #make sure that *every* original start or end point in x is either a start or end point
   #for the new start/end variables in xd,
  #if those original start/end dates overlap with
    #the new start dates at all (within groups)
  #if value is between start and end but not exactly start or end,
   #this means new start/end is missing a cut
  x_l <- melt(x,id.vars=c("id1","id2"),measure.vars=c("start0","end0"))
  x_l[,value2:=value]
  setkey(x_l, id1,id2,value,value2)
  setkey(xd, id1,id2,start,end)
  vr <- foverlaps(x_l,xd)
  expect_true(vr[, list(V1=value %in% c(start,end)),by=1:nrow(vr)][,all(V1)])



})

#####test CJ.dt for data.tables ####
test_that("CJ.dt", {

  ##simple example
  X <- data.table(x1=1:2,x2=2:3)
  Y <- data.table(y1=4:6,y2=5:7)
  expect_equal(CJ.dt(X,Y),
                      data.table(x1=rep(1:2,times=3),
                                 x2=rep(2:3,times=3),
                                 y1=rep(4:6,each=2),
                                 y2=rep(5:7,each=2)
                      ))
  ##check grouping
  X <- data.table(
    x1 = 1:8,
    x2 = 2:9,
    id1 = rep(1:2, each = 2),
    id2 = rep(1:2, times = 2)
  )
  Y <- data.table(
    y1 = 4:11,
    y2 = 5:12,
    id1 = rep(1:2, each = 2),
    id2 = rep(1:2, times = 2)
  )


  ##CJ.dt with groups should be like iterating through each combination
  #of groups and doing a cartesian merge/grid expand within that:
  templ <- list()
  counter <- 1
  for(i in unique(X$id1)){
    for(j in unique(X$id2)){
      templ[[counter]] <- CJ.dt(X[id1==i&id2==j],Y[id1==i&id2==j,list(y1,y2)])
      counter <- counter +1
    }
  }
  expect_equal(
    CJ.dt(X,Y,groups=c("id1","id2")),
    rbindlist(templ)
  )




  X <- data.table(
    x1 = 1,
    x2 = 2,
    id1 = 1,
    id2 = 1
  )
  Y <- data.table(
    y1 = 4:11,
    y2 = 5:12,
    id1 = rep(1:2, each = 2),
    id2 = rep(1:2, times = 2)
  )


  templ <- list()
  counter <- 1
  for(i in unique(X$id1)){
    for(j in unique(X$id2)){
      templ[[counter]] <- CJ.dt(X[id1==i&id2==j],Y[id1==i&id2==j,list(y1,y2)])
      counter <- counter +1
    }
  }
  expect_equal(
    CJ.dt(X,Y,groups=c("id1","id2")),
    rbindlist(templ)
  )

})


#### intervalaverage state restoration ####
test_that("intervalaveraging restores state", {
  #averaging intervals longer than observed period
  set.seed(72)
  a_start_date <- seq(structure(10590L, class = c("IDate", "Date")),
                      structure(17163L, class = c("IDate", "Date")),by=7L)

  a0 <- CJ(id1=1:3,id2=1:100, start_date=a_start_date)
  a0[, end_date:=start_date+6L]
  a0[, value1:=rnorm(.N)]
  a0[, value2:=rnorm(.N)]

  b0_temp <- data.table(start_date=as.IDate(paste0(1999:2017,"-01-01")),
                        end_date=as.IDate(paste0(1999:2017,"-12-31")))
  b0 <- CJ.dt(b0_temp, unique(a0[,list(id1,id2)]) )

  #two groups in x,two values
  a0[,neworder:=sample(1:.N)]
  setkey(a0,"neworder")
  a0_original <- copy(a0)

  b0[,neworder:=sample(1:.N)]
  setkey(b0,"neworder")
  b0_original <- copy(b0)

  q0_1 <- intervalaverage(x=a0,
                          y=b0,
                          interval_vars=c("start_date","end_date"),
                          value_vars=c("value1","value2"),
                          group_vars=c("id1","id2"))

  expect_equal(a0,a0_original)
  expect_equal(b0,b0_original)


})



test_that("intervalaveraging restores state even when it throws an error", {
  #averaging intervals longer than observed period
  set.seed(72)
  a_start_date <- seq(structure(10590L, class = c("IDate", "Date")),
                      structure(17163L, class = c("IDate", "Date")),by=7L)

  a0 <- CJ(id1=1:3,id2=1:100, start_date=a_start_date)
  a0[, end_date:=start_date+6L]
  a0[, value1:=rnorm(.N)]
  a0[, value2:=rnorm(.N)]
  #make it overlapping so it throws an error
  a0[id1==1&id2==1,start_date:=min(a0$start_date)]
  a0[id1==1&id2==1,end_date:=max(a0$end_date)]

  expect_true(is.overlapping(a0,interval_vars=c("start_date","end_date"),
                              group_vars=c("id1","id2")
                              ),"overlap")

  b0_temp <- data.table(start_date=as.IDate(paste0(1999:2017,"-01-01")),
                        end_date=as.IDate(paste0(1999:2017,"-12-31")))
  b0 <- CJ.dt(b0_temp, unique(a0[,list(id1,id2)]) )

  #two groups in x,two values
  a0[,neworder:=sample(1:.N)]
  setkey(a0,"neworder")
  a0_original <- copy(a0)

  b0[,neworder:=sample(1:.N)]
  setkey(b0,"neworder")
  b0_original <- copy(b0)

  expect_error(q0_1 <- intervalaverage(x=a0,
                          y=b0,
                          interval_vars=c("start_date","end_date"),
                          value_vars=c("value1","value2"),
                          group_vars=c("id1","id2")),"replicate/duplicate intervals")

  expect_equal(a0,a0_original)
  expect_equal(b0,b0_original)


})

##### intervalaverage and missingness #####

test_that("intervalaverage only partially observed schedule and missingness:",{

  x <- data.table(start=1L,end=9L,value=3)
  y <- data.table(start=1L,end=10L)
expect_equal(
  intervalaverage(x,y,interval_vars=c("start","end"),value_vars="value")$value,
  as.numeric(NA))

expect_equal(
  intervalaverage(x,y,interval_vars=c("start","end"),value_vars="value",required_percentage = 90)$value,
  3)
expect_equal(
  intervalaverage(x,y,interval_vars=c("start","end"),value_vars="value",required_percentage = 89)$value,
  3)
expect_equal(
  intervalaverage(x,y,interval_vars=c("start","end"),value_vars="value",required_percentage = 0)$value,
  3)

})


test_that("intervalaverage, entirely missing interval:",{

  x <- data.table(start=1L,end=10L,value=3)
  y <- data.table(start=20L,end=30L)
  expect_equal(
    intervalaverage(x,y,interval_vars=c("start","end"),value_vars="value")$value,
    as.numeric(NA)
  )

  expect_false(
    is.nan(intervalaverage(x,y,interval_vars=c("start","end"),value_vars="value")$value)
  )


})




test_that("intervalaveraging interval_vars can't have names",{

  #these should throw errors!



  #averaging intervals longer than observed period
  set.seed(72)
  a_start_date <- seq(structure(10590L, class = c("IDate", "Date")),
                      structure(17163L, class = c("IDate", "Date")),by=7L)

  a0 <- CJ(id1=1:3,id2=1:100, start_date=a_start_date)
  a0[, end_date:=start_date+6L]
  a0[, value1:=rnorm(.N)]
  a0[, value2:=rnorm(.N)]

  b0_temp <- data.table(start_date=as.IDate(paste0(1999:2017,"-01-01")),
                        end_date=as.IDate(paste0(1999:2017,"-12-31")))
  b0 <- CJ.dt(b0_temp, unique(a0[,list(id1,id2)]) )


  expect_error(intervalaverage(x=a0,
                               y=b0,
                               interval_vars=c(test="start_date",test2="end_date"),
                               value_vars=c("value1","value2"),
                               group_vars=c("id1","id2")))

  expect_error(intervalaverage(x=a0,
                               y=b0,
                               interval_vars=c(test="start_date","end_date"),
                               value_vars=c("value1","value2"),
                               group_vars=c("id1","id2")))


  expect_error(intervalaverage(x=a0,
                               y=b0,
                               interval_vars=c("start_date",test="end_date"),
                               value_vars=c("value1","value2"),
                               group_vars=c("id1","id2")))



  expect_error(intervalaverage(x=a0,
                               y=b0,
                               interval_vars=c(test="start_date","end_date"),
                               value_vars=c("value1","value2"),
                               group_vars=c(test="id1","id2")))



  expect_error(intervalaverage(x=a0,
                          y=b0,
                          interval_vars=c("start_date","end_date"),
                          value_vars=c("value1","value2"),
                          group_vars=c(test="id1","id2")))



  expect_error(intervalaverage(x=a0,
                               y=b0,
                               interval_vars=c("start_date","end_date"),
                               value_vars=c("value1","value2"),
                               group_vars=c("id1",test="id2")))




})



####test intervalaveraging function #########
test_that("intervalaveraging", {

  ##averaging intervals where groups=NULL by default
  set.seed(32)
  a0.0 <- data.table(start=seq(1L,100L,by=10L),value1=rnorm(10))
  a0.0[,end:=start+9L]
  b0.0 <- data.table(start=1L,end=25L)
  q0.01 <- intervalaverage(x=a0.0,
                           y=b0.0,
                           interval_vars=c("start","end"),
                           value_vars=c("value1"))


  q0.02 <- intervalaverage:::interval_weighted_avg_slow_f(x=a0.0,
                                        y=b0.0,
                                        interval_vars=c("start","end"),
                                        value_vars=c("value1"))

  expect_equal(q0.01,q0.02)
  expect_equal(nrow(q0.01),nrow(b0.0))





  #averaging intervals longer than observed period
  set.seed(72)
  a_start_date <- seq(structure(10590L, class = c("IDate", "Date")),
                      structure(17163L, class = c("IDate", "Date")),by=7L)

  a0 <- CJ(id1=1:3,id2=1:100, start_date=a_start_date)
  a0[, end_date:=start_date+6L]
  a0[, value1:=rnorm(.N)]
  a0[, value2:=rnorm(.N)]

  b0_temp <- data.table(start_date=as.IDate(paste0(1999:2017,"-01-01")),
                        end_date=as.IDate(paste0(1999:2017,"-12-31")))
  b0 <- CJ.dt(b0_temp, unique(a0[,list(id1,id2)]) )

  ###two groups in x,two values####

  q0_1 <- intervalaverage(x=a0,
                          y=b0,
                          interval_vars=c("start_date","end_date"),
                          value_vars=c("value1","value2"),
                          group_vars=c("id1","id2"))


  q0_2 <- intervalaverage:::interval_weighted_avg_slow_f(x=a0,
                                       y=b0,
                                       interval_vars=c("start_date","end_date"),
                                       value_vars=c("value1","value2"),
                                       group_vars=c("id1","id2"))

  expect_equal(q0_1,q0_2)
  expect_equal(nrow(q0_1),nrow(b0))


  ##make sure isolateoverlaps throws an error when "i." variables are provided
  a_overlap1 <- CJ(id1=1:3,id2=1:100, start_date=a_start_date)
  a_overlap1[, end_date:=start_date+10L]
  a_overlap1[, i.start_date:=TRUE]
  expect_error(
    isolateoverlaps(
      x=a_overlap1,
      interval_vars=c("start_date","end_date"),
      group_vars=c("id1","id2")
    ))





  #simple example:
  #in this example, b is a regular schedule over which we'd like to compute averages of values from a
  #a is also on a regular schedule that does not fit neatly into b's schedule
  #(some intervals in a partially overlap with two intervals in b )
  #additionally, the first interval in a doesn't overlap with
  #intervals in b at all and is entirely non-adjoining. since we're looking for averages over intervals in b,
  #we don't care about this interval and it should be left out of the results

  #b contains an interval that is observed in a for almost all but not the entire period
  #(partial overlap, we probably want an average for this period)
  #b also contains an interval that is barely observed except for a small overlap
  #(partial overlap, we probably want NA for the average of this period)



  #b finishes with an interval which contains no overlap with any intervals in a
  #we want averges for all intervals in b--
  #so the result should have a row for intervals in b whether there's a complete, partial, or no overlap with a:
  #intervals in a that don't overlap at all with intervals in b should not be returned--
  #because we have no interest in intervals not specified in b

  set.seed(12380)
  a <- CJ(id=1:2,id2=1:2,start=c(-13L,seq(1L,36L,by=7L)))
  a[, end:=start+6L]
  a[, value:=rbinom(.N,5,prob=.5)]
  a[, value2:=rbinom(.N,5,prob=.5)]
  #randomly insert some missing values in one column:
  a[as.logical(rbinom(.N,1,prob=.3)), value:=NA]

  b_temp <- data.table(start=seq(0L,56L,by=14L))
  b_temp[,end:=start+13L]

  b <- CJ.dt(b_temp, unique(a[,list(id,id2)]) )
  q1 <- intervalaverage(x=a,
                        y=b,
                        interval_vars=c("start","end"),
                        group_vars=c("id","id2"),
                        value_vars=c("value","value2")
  )
  q2 <- intervalaverage:::interval_weighted_avg_slow_f(x=a,
                                                       y=b,
                                                       interval_vars=c("start","end"),
                                                       group_vars=c("id","id2"),
                                                       value_vars=c("value","value2"))
  expect_equal(q1,q2)
  expect_equal(nrow(q1),nrow(b))





  ####group_vars=NULL
  q1n <- intervalaverage(x=a[id==1&id2==1],
                         y=b[id==1&id2==1],
                         interval_vars=c("start","end"),
                         group_vars=NULL,
                         value_vars=c("value","value2")
  )


  q2n <- intervalaverage:::interval_weighted_avg_slow_f(x=a[id==1&id2==1],
                                                        y=b[id==1&id2==1],
                                                        interval_vars=c("start","end"),
                                                        group_vars=NULL,
                                                        value_vars=c("value","value2"))

  expect_equal(q1n,q2n)
  setkey(q1,NULL)
  setkey(q1n,NULL)
  expect_equal(q1n,q1[id==1&id2==1,.SD,.SDcols=-c("id","id2")])
  expect_equal(nrow(q1),nrow(b))





  ##test missingness when required_percentage is not 100


  qm1 <- intervalaverage(x=a,
                         y=b,
                         interval_vars=c("start","end"),
                         group_vars=c("id","id2"),
                         value_vars=c("value","value2"),
                         required_percentage = 50
  )


  qm2 <- intervalaverage:::interval_weighted_avg_slow_f(x=a,
                                                        y=b,
                                                        interval_vars=c("start","end"),
                                                        group_vars=c("id","id2"),
                                                        value_vars=c("value","value2"),
                                                        required_percentage=50)
  expect_equal(qm1,qm2)
  expect_equal(nrow(qm1),nrow(b))






  ##averaging intervals shorter than observed period & includes dates not in observed period at all ####


  b2_temp <- data.table(start=seq(0L,56L,by=3L))
  b2_temp[,end:=start+2L]
  b2 <- CJ.dt(b2_temp, unique(a[,list(id,id2)]) )



  q2_1 <- intervalaverage(x=a,
                          y=b2,
                          interval_vars=c("start","end"),
                          group_vars=c("id","id2"),
                          value_vars=c("value","value2"),
                          required_percentage = 50
  )


  q2_2 <- intervalaverage:::interval_weighted_avg_slow_f(x=a,
                                                         y=b2,
                                                         interval_vars=c("start","end"),
                                                         group_vars=c("id","id2"),
                                                         value_vars=c("value","value2"),
                                                         required_percentage = 50
  )



  expect_equal(q2_1,q2_2)
  expect_equal(nrow(q2_1),nrow(b2))








  #1 groups in both x AND y with different desired averagin periods for each group
  #(e.g., different desired averaging period for each group))


  b3_11 <- data.table(id=1L,id2=1L,start=seq(0L,56L,by=14L))
  b3_11[,end:=start+13L]
  b3_12 <- data.table(id=1L,id2=2L,start=seq(0L,54L,by=3L))
  b3_12[,end:=start+2L]
  b3_21 <- data.table(id=2L,id2=1L,start=seq(3L,43L,by=20L))
  b3_21[,end:=start+19L]
  b3_22 <- data.table(id=c(2L),id2=c(2L),start=c(5L,100L),
                      end=c(12L,101L))

  b3 <- rbind(b3_11,b3_12,b3_21,b3_22)
  #note that for id=1,id2=1 and id=1,id2=2 these should be the same as
  #subsets of qm1 and q2_1 respectively



  q3_1 <- intervalaverage(x=a,
                          y=b3,
                          interval_vars=c("start","end"),
                          group_vars=c("id","id2"),
                          value_vars=c("value","value2"),
                          required_percentage = 50
  )


  q3_2 <- intervalaverage:::interval_weighted_avg_slow_f(x=a,
                                                         y=b3,
                                                         interval_vars=c("start","end"),
                                                         group_vars=c("id","id2"),
                                                         value_vars=c("value","value2"),
                                                         required_percentage = 50
  )
  expect_equal(q3_1[id==1&id2==1],qm1[id==1&id2==1])
  expect_equal(q3_1[id==1&id2==2],q2_1[id==1&id2==2])

  expect_equal(q3_1,q3_2,check.attributes=FALSE)
  expect_equal(nrow(q3_2),nrow(b3))






  #periods that have overlaps in y--this should work
  b_overlap_temp <- rbind(b2_temp,data.table(start=3L,end=18L))

  b_overlap <- CJ.dt(b_overlap_temp, unique(a[,list(id,id2)]) )

  q_overlap_y1 <- intervalaverage(x=a,
                                  y=b_overlap,
                                  interval_vars=c("start","end"),
                                  group_vars=c("id","id2"),
                                  value_vars=c("value","value2"),
                                  required_percentage = 50
  )


  q_overlap_y2 <- intervalaverage:::interval_weighted_avg_slow_f(x=a,
                                                                 y=b_overlap,
                                                                 interval_vars=c("start","end"),
                                                                 group_vars=c("id","id2"),
                                                                 value_vars=c("value","value2"),
                                                                 required_percentage = 50
  )


  expect_equal(q_overlap_y1,q_overlap_y2)
  expect_equal(nrow(q_overlap_y1),nrow(b_overlap))






  #periods that have partial overlaps in x--this should return an error
  set.seed(2340)
  a_overlap <- CJ(id=1:2,id2=1:2,start=c(-13L,seq(1L,36L,by=7L),2L))
  a_overlap[, end:=start+6L]
  a_overlap[, value:=rbinom(.N,5,prob=.5)]

  expect_error(
    intervalaverage(a_overlap,b,interval_vars=c("start","end"),
                    value_vars=c("value"),
                    group_vars=c("id","id2"))
  )




  ##periods that have exact overlaps but no partial overlaps--this return error##
  a_exact_overlap <- data.table(start_date=c(1L,1L,4L),end_date=c(3L,3L,10L),value1=c(0,5,10))
  b_exact_overlap <- data.table(start_date=2L:3L,end_date=8L:9L)

  expect_error(intervalaverage(a_exact_overlap,b_exact_overlap,interval_vars=c("start_date","end_date"),
                             value_vars=c("value1")))


  expect_error(intervalaverage:::interval_weighted_avg_slow_f(a_exact_overlap,b_exact_overlap,interval_vars=c("start_date","end_date"),
                                                              value_vars=c("value1")))




  #misspecified order of interval_vars should return an error:
  expect_error(
    intervalaverage(a,b,interval_vars=c("end_date","start_date"),
                    value_vars=c("value1","value2"),
                    group_vars=c("id1","id2")),
    )








  ####realistic example with overlaping values: deoverlap them then average to a period:
  set.seed(93450)
  a_overlap1 <- CJ(id1=1:3,id2=1:100, start=a_start_date)
  a_overlap1[, end:=start+10L]
  a_overlap1[, value1:=rnorm(.N)]
  a_overlap1[, value2:=rnorm(.N)]


  expect_true(
    is.overlapping(x=a_overlap1,
                   interval_vars = c("start","end"),
                   group_vars=c("id1","id2"))
  )

  a_overlap1_removed <- isolateoverlaps(
    x=a_overlap1,
    interval_vars=c("start","end"),
    group_vars=c("id1","id2"),
    interval_vars_out = c("start_date","end_date")
  )

  #average duplicate intervals together
  a_overlap1_removed <-
    a_overlap1_removed[,list(value1=mean(value1),value2=mean(value2)),
                       by=list(id1,id2,start_date,end_date)]

  expect_false(
    is.overlapping(a_overlap1_removed,
                   interval_vars = c("start_date","end_date"),
                   group_vars=c("id1","id2"))
  )

  #insert missingness
  a_overlap1_removed[sample(1:nrow(a_overlap1_removed),size=floor(.2*nrow(a_overlap1_removed))),
                     `:=`(value1=NA,value2=NA)]


  b_overlap1 <- CJ(id1=1:3,id2=1:100)[,
                                      rep(TRUE,sample(1:10,size=1)),
                                      by=list(id1,id2)]
  b_overlap1[, start_date:=sample(a_start_date,nrow(b_overlap1),replace=TRUE)+sample(-5L:5L,nrow(b_overlap1),replace=TRUE)]
  b_overlap1[, end_date:=start_date+sample(0:10,nrow(b_overlap1),replace=TRUE)]




  system.time(realistic1 <-   intervalaverage(x=a_overlap1_removed,
                                              b_overlap1,
                                              interval_vars=c("start_date","end_date"),
                                              value_vars=c("value1","value2"),
                                              group_vars=c("id1","id2"),
                                              skip_overlap_check=FALSE))

  system.time(realistic2 <-
                intervalaverage:::interval_weighted_avg_slow_f(a_overlap1_removed,
                                                               b_overlap1,
                                                               interval_vars=c("start_date","end_date"),
                                                               value_vars=c("value1","value2"),
                                                               group_vars=c("id1","id2"),
                                                               skip_overlap_check=FALSE))

  expect_equal(realistic1,realistic2)


})


test_that("intervalaveraging big group 1", {
  skip_on_cran()
  ##more realism
  set.seed(12323)
  n <- 1e5
  x <- matrix(as.integer(round(runif(n=n*2, 0L,1000L))),ncol=2)
  x <- as.data.table(t(apply(x,1,sort)))
  setnames(x,names(x),c("start0","end0"))
  x[, id1:=rbinom(n,3,prob=.3)]
  x[, id2:=rbinom(n,7,prob=.5)]
  x[,value1:=rnorm(n)]
  x[,value2:=rnorm(n)]
  setkey(x, id1,id2,start0,end0)


  a <- isolateoverlaps(x,interval_vars=c("start0","end0"),group_vars=c("id1","id2"),
                       c("start","end")
  )

  #average duplicate intervals together
  a <-
    a[,list(value1=mean(value1),value2=mean(value2)),by=list(id1,id2,start,end)]

  #insert missingness
  a[sample(1:nrow(a),size=floor(.2*nrow(a))),`:=`(value1=NA,value2=NA)]
  b <- matrix(as.integer(round(runif(n=n*2, 0L,1000L))),ncol=2)
  b <- as.data.table(t(apply(b,1,sort)))
  setnames(b,names(b),c("start","end"))
  b[, id1:=rbinom(n,3,prob=.3)]
  b[, id2:=rbinom(n,7,prob=.5)]



  ##for some reason, expect_warning make all the arguments copied rather than passed by ref
  ##here I want to test to make sure that when there are overlapping intervals in y, it still returns y to its original state
  #but in order to test this, I can't wrap in expect_warning

  options(warn=-1)
  bleh <- intervalaverage(x=a,y=b,interval_vars=c("start","end"),
                          value_vars=c("value1","value2"),
                          group_vars=c("id1","id2"),
                          skip_overlap_check=FALSE)


  expect_false("rowindex" %in% names(b))
  options(warn=0)


  expect_warning(zzz1 <- intervalaverage(x=a,y=b,interval_vars=c("start","end"),
                          value_vars=c("value1","value2"),
                          group_vars=c("id1","id2"),
                          skip_overlap_check=FALSE),
                 "removing these duplicated rows automatically")

  expect_warning(zzz2 <- intervalaverage:::interval_weighted_avg_slow_f(a,b,
                                                         interval_vars=c("start","end"),
                                                         value_vars=c("value1","value2"),
                                                         group_vars=c("id1","id2"),
                                                         skip_overlap_check=FALSE),
                 "removing these duplicated rows automatically")

  expect_equal(zzz1,zzz2)

})



test_that("intervalaveraging big group 2", {
  skip_on_cran()

  ####large dataset that's non-overlapping
  set.seed(18)
  az_start_date <- seq(structure(0L, class = c("IDate","Date")),
                       structure(1e5L, class = c("IDate","Date")),by=14)

  az <- CJ(id1=1:100, start_date=az_start_date)
  az[, end_date:=start_date+13L]
  az[, value1:=rnorm(.N)]
  az[, value2:=rnorm(.N)]

  bz_start_date <- seq(structure(2L, class = c("IDate","Date")),
                       structure(1e5L, class = c("IDate","Date")),by=7)

  bz <- CJ(id1=1:100, start_date=bz_start_date)
  bz[, end_date:=start_date+6L]

  zzbig1 <- intervalaverage(x=az,y=bz,interval_vars=c("start_date","end_date"),
                            value_vars=c("value1","value2"),
                            group_vars="id1",
                            skip_overlap_check=TRUE)

  zzbig2 <- intervalaverage:::interval_weighted_avg_slow_f(az,bz,interval_vars=c("start_date","end_date"),
                                                           value_vars=c("value1","value2"),
                                                           group_vars="id1",
                                                           skip_overlap_check=TRUE)

  expect_equal(zzbig1,zzbig2)



  ####test loop
  zzbig1_list <- list()
  uid <- unique(az$id1)
  for(i in 1:length(uid)){
    zzbig1_list[[i]] <- intervalaverage(az[id1==uid[i]],bz[id1==uid[i]],
                                        interval_vars=c("start_date","end_date"),
                                        value_vars=c("value1","value2"),
                                        group_vars=NULL,
                                        skip_overlap_check=TRUE)

    zzbig1_list[[i]][,id1:=uid[i]]
  }
  zzbig1l <- rbindlist(zzbig1_list)
  setcolorder(zzbig1l, "id1")
  setkey(zzbig1l,NULL)
  setkey(zzbig1,NULL)
  expect_equal(zzbig1, zzbig1l)

  ###
  test_x <- data.table(id1=c(1,1,1:4),region=c(1,1,1,1,2,NA),value=c(1:6),start_date=c(1L, 6L,11L,6L,11L,1L),end_date=c(5L, 10L,15L, 10L,15L,5L))
  test_y_temp <- data.table(id1=c(1,1,200),start_date=c(1L,50L,1L),end_date=c(7L,60L,7L))

  test_y <- CJ.dt(test_y_temp, unique(test_x[,list(region)]) )


  ff1 <- intervalaverage(test_x, test_y,
                         interval_vars=c("start_date","end_date"),
                         group_vars = c("id1","region"), value_vars=c("value"),required_percentage = 0)


  ff2 <- intervalaverage:::interval_weighted_avg_slow_f(test_x, test_y,
                                                        interval_vars=c("start_date","end_date"),
                                                        group_vars = c("id1","region"), value_vars=c("value"),required_percentage = 0)


  expect_equal(ff1,ff2)


  expect_warning(fg1 <- intervalaverage(x=test_x,
                         y=test_y,
                         interval_vars=c("start_date","end_date"),
                         group_vars = c("id1"),
                         value_vars=c("value"),
                         required_percentage = 0
  ),"removing these duplicated rows automatically")


  expect_warning(fg2 <-
    intervalaverage:::interval_weighted_avg_slow_f(test_x,
                                                   test_y,
                                                   interval_vars=c("start_date","end_date"),
                                                   group_vars = c("id1"),
                                                   value_vars=c("value"),
                                                   required_percentage = 0
    ),"removing these duplicated rows automatically")



  expect_equal(fg1,fg2)

})
