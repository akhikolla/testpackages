#this used to compare performance of various methods but since optimizing this makes
 #less sense to be a test
#keeping this in since it's a nice example of the different approaches to intervalaveraging
 #using data.table
#and also is a simple reproducible example of how the averaging works,
#separated from all the abstraction in the function
#plus it acts as unit test for the data.table package around the specific operations I'm doing.

test_that("performance", {
  skip_on_cran()

  uid <- 1L:1000L
  x <- CJ(id=uid,start=seq(1L,2000L,by=7L))
  x[, end:=start+6L]
  x[, value:=rnorm(.N)]

  ##when averaging intervals are smaller than observation intervals
    #the foverlaps approach creates a large allocation
  y <- CJ(id=uid, start=seq(1L,1991L,by=3L))
  y[,end:=start+9L]

  setkey(x,id,start,end)
  setkey(y,id,start,end)



  ##
  time_package_approach <- system.time({
    out_package_approach <-
      intervalaverage(x,
                      y,
                      interval_vars=c("start","end"),
                      group_vars=c("id"),
                      value_vars = "value"
      )
  })




  ##for comparison--naive expansion of intervals approach--large allocation by definition:
  time0 <- system.time({

    x_expanded <- x[,list(id=id, time=start:end,value=value),by=1:nrow(x)]
    x_expanded[, time2:=time]
    setkey(x_expanded, id,time,time2)
    z_0 <- foverlaps(y,x_expanded) #~6 million rows

    out0 <- z_0[,list(avg_value=mean(value)),by=c("id","start","end")]

  })
  ##foverlaps approach, using GFORCE but has large intermediate allocation (z)

  time1 <- system.time({
    z <- foverlaps(y,x) # ~1.5 million rows
    z[, new_start:=pmax(start,i.start)]
    z[, new_end:=pmin(end,i.end)]

    #calculate the weighted average in parts to make use of GFORCE
    z[, length_new_interval:=new_end-new_start+1L]
    z[, product:=value*length_new_interval]
    out1 <- z[, list(sum_lengths=sum(length_new_interval),sum_products=sum(product)),
      by=c("id","i.start","i.end")]


    out1[, avg_value:=sum_products/sum_lengths]
    out1[, sum_lengths:=NULL]
    out1[, sum_products:=NULL]
    setnames(out1, c("i.start","i.end"),c("start","end"))
  })

  #non-equijoin approach, using .EACHI to avoid large intermediate allocation


  #because I can never understand what a non-equi join does with the columns:
  setnames(x,c("start","end"),c("xstart","xend"))
  setnames(y,c("start","end"),c("ystart","yend"))
  x[, xstart2:=xstart]
  x[, xend2:=xend]
  y[, ystart2:=ystart]
  y[, yend2:=yend]

  time2a <- system.time({
  #x[i=y,by=.EACHI, on=c("id","xend>=ystart","xstart<=yend"),nomatch=NULL,allow.cartesian=TRUE]
    out2a <- x[i=y,by=.EACHI, on=c("id","xend>=ystart","xstart<=yend"),nomatch=NULL,
      j=list(avg_value=weighted.mean(x=value,
                           w=pmin(xend2,yend2)-pmax(xstart2,ystart2)+1L)
             )
      ]
  })


  #maybe the above is slow because weighted.mean is poorly optimized?
  #calculate the above weighted mean manually:

  time2b <- system.time({
    #z2 <- x[i=y, on=c("id","xend>=ystart","xstart<=yend"),nomatch=NULL,allow.cartesian=TRUE]
    #z2[, new_start:=pmax(xstart2,ystart2)]
    #z2[, new_end:=pmin(xend2,yend2)]
    #z2[, length_new_interval:=new_end-new_start+1L]

    out2b <- x[i=y,by=.EACHI, on=c("id","xend>=ystart","xstart<=yend"),nomatch=NULL,
              j=list(avg_value={
                length_new_interval <- pmin(xend2,yend2)-pmax(xstart2,ystart2)+1L
                sum(value*length_new_interval)/sum(length_new_interval)
              }
              )
              ]
  })


  #it turns out that basically none of that time is actually spent calclulating the mean
   #pmin()/pmax() is much slower using .EACHI rather than directly on the large intermediate allocation:
  time2_onlyintervallength <- system.time({
    #x[i=y,by=.EACHI, on=c("id","xend>=ystart","xstart<=yend"),nomatch=NULL,allow.cartesian=TRUE]
    out2_intervallength <- x[i=y,by=.EACHI, on=c("id","xend>=ystart","xstart<=yend"),nomatch=NULL,
               j=list(length_new_interval= pmin(xend2,yend2)-pmax(xstart2,ystart2)+1L)
               ]
  })





  #sanity check: compare results
   #weirdness non-equi join nameswitch
  setnames(out2a, c("xend","xstart"),c("start","end"))
  setnames(out2b, c("xend","xstart"),c("start","end"))
  setcolorder(out2a, c("id","start","end","avg_value"))
  setcolorder(out2b, c("id","start","end","avg_value"))
  setkey(out2a,id,start,end)
  setkey(out2b,id,start,end)


  expect_equal(out_package_approach[, list(id,start,end,avg_value=value)], out0)
  expect_equal(out0,out1)
  expect_equal(out1,out2a)
  expect_equal(out2a,out2b)



})


