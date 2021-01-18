#' @useDynLib intervalaverage, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL


create_unused_name <- function(x,reserved_cols){
  for(i in 1:length(x)){
    while(x[i] %in% reserved_cols){
      x[i] <- paste0("i.",x[i])
    }
  }
  x
}



#' grid expand an arbitrary number of data.tables
#'
#' similar to data.table::CJ and base::expand.grid except for rows of data.tables.
#'
#' CJ.dt computes successive cartesian join over rows of each table
#' paying no attention to whatever the tables are keyed on.
#'
#' @param ... data.tables
#' @param groups a character vector corresponding to
#' column names of grouping vars in all of the data.tables
#' @examples
#' #' CJ.dt(data.table(c(1,2,2),c(1,1,1)),data.table(c("a","b"),c("c","d")))
#' #If you want to expand x to unique values of a non-unique columns in y
#' x <- data.table(c(1,2,3),c("a","b","b"))
#' y <- data.table(id=c(1,2,2,1,3),value=c(2,4,1,7,3))
#' z <- CJ.dt(x, y[,list(id=unique(id))])
#' #if you want to merge this back to y
#' y[z,on="id",allow.cartesian=TRUE] #or z[y,on="id",allow.cartesian=TRUE]
#' @import data.table
#' @export
CJ.dt <- function(...,groups=NULL) {
  l = list(...)
  stopifnot( all(sapply(l,is.data.table)) )
  l <- lapply(l,copy) #take a copy of each table
    #copying is probably not going to cause memory or timing issues here
    #because if you don't have the space to copy the individual tables there's
      #no way you'll be able to grid expand them together
     #with that said, the rest of this function is written to try to return each table to its
     #original state so could be modified to remove the copying
     #if this is slowing a particular task down
     #edit out the copy at your own risk and know that if the function returns an error
       #all your tables will likely be modified because they have an extra column
  EVAL <- function(...)eval(parse(text=paste0(...)))
  if(any(sapply(l,nrow)==0)){stop("one or more data.tables have no rows")}
  stopifnot(all(sapply(l,data.table::is.data.table)))


  #while kvar is in names of any of the data.tables, keep prepending i. to it until it's not
  kvar <- create_unused_name("k",unlist(lapply(l,names)))

  invars <- create_unused_name(paste0("in",1:length(l)),
                             unlist(lapply(l,names)))

  for(i in 1:length(l)){
    l[[i]][,(kvar):=1]
    l[[i]][,(invars[i]):=TRUE]
    data.table::setkeyv(l[[i]],c(kvar,groups))
  }

  mymerge = function(x,y) x[y, allow.cartesian = TRUE]
  out <- Reduce(mymerge,l)
  out[,(kvar):=NULL]

  for(i in 1:length(l)){
    l[[i]][,(kvar):=NULL]
    l[[i]][,(invars[i]):=NULL]
  }

  out <- EVAL("out[",paste0(paste0("!is.na(",invars,")"),collapse="&"),"]")
  out[,(invars):=NULL]
  out[]
}



#' Test for self-overlap
#'
#' Test whether a data.table contains intervals which partially or completely overlap
#' with other intervals in different rows, possibly within groups
#'
#' @param x A data.table with two columns defining closed intervals (see also interval_vars).
#' @param interval_vars A length-2 character vector corresponding to column names of x which designate
#' the closed (inclusive) starting and ending intervals. The column name
#' specifying the lower-bound column must be specified first.
#' @param group_vars NULL or a character vector corresponding to column names of x.
#' overlap checks will occur within groups defined by the columns specified here.
#' @param verbose prints additional information, default is FALSE
#' @return length-1 logical vector. TRUE if there are overlaps, FALSE otherwise.
#'
#' @examples
#' x <- data.table(start=c(1,2),end=c(3,4))
#' is.overlapping(x,c("start","end")) #the interval 1,3 overlaps with the interval 2,4
#' y <- data.table(start=c(1,3),end=c(2,4))
#' is.overlapping(y,c("start","end")) #the interval 1,2 doesn't overlap other intervals in y
#' z <- data.table(start=c(1,3,1,2),end=c(2,4,3,4),id=c(1,1,2,2))
#' is.overlapping(z,c("start","end"),"id")
#' @export
is.overlapping <- function(x,interval_vars,group_vars=NULL,verbose=FALSE){
  if("._.irow"%in%names(x)){
    stop("._.irow cannot be in names of x. rename this column. If you didn't expect this
         column to be here, it may because this function previously crashed. You can
         safely delete this column. Please submit a bug report to the github repo.")
  }
  ._.irow <- i.._.irow <-  NULL # due to NSE notes in R CMD check
  stopifnot(is.data.table(x))

  is_not_preferred_keyx <- !identical(key(x), c(group_vars,interval_vars))
  if(is_not_preferred_keyx){
    statex <- savestate(x)
    if(verbose){message("setkeyv(x,c(group_vars,interval_vars)) prior to calling is.overlap is recommended to save unnecessary row reordering")}
    on.exit(setstate(x,statex))
  }

    data.table::setkeyv(x,c(group_vars,interval_vars))
    stopifnot(!"._.irow"%in%names(x))
    stopifnot(!"i_._.irow"%in%names(x))
    x[,._.irow:=1:nrow(x)]
    z <- data.table::foverlaps(x,x)
    out <- z[._.irow!=i.._.irow]

  nrow(out)!=0
}



#' time-weighted average of values measured over intervals
#'
#' \code{intervalaverage} takes values recorded over
#' non-overlapping intervals and averages them to defined intervals, possibly within
#' groups (individuals/monitors/locations/etc).  This function could be used to take averages over long
#' intervals
#' of values measured over short intervals and/or to take short "averages" of values measured over
#' longer intervals (ie, downsample without smoothing). Measurement intervals and averaging intervals need
#' not align. In the event that an averaging interval contains more than one measurement interval,
#' a weighted average is calculated (ie each measurement is weighted on the duration of its interval's
#' overlap with the averaging period interval).
#'
#' All intervals are treated as closed (ie inclusive of the start and end values in interval_vars)
#'
#'
#' x and y are not copied but rather passed by reference to function internals
#' but the order of these data.tables is restored on function completion or error,
#'
#' When required_percentage is less than 100, xminstart and xmaxend may be useful to
#' determine whether an average meets specified coverage requirements in terms of not
#' just percent of missingness but whether values are represented through the range of the y interval
#'
#'
#'
#' @param x a data.table containing values measured over intervals. see
#' `interval_vars` parameter
#' for how to specify interval columns and `value_vars` for how to
#'  specify value columns.
#' intervals in `x` must must be completely non-overlapping within
#' groups defined by group_vars. if `group_vars` is specified (non-`NULL`), `x` must
#' also contain columns specified in `group_vars`.
#' @param y a data.table object containing intervals over which averages of `x` values should be computed.
#' averaging intervals in `y`, unlike measurement intervals in `x`, may be overlapping within groups.
#' if `group_vars` is specified (non-`NULL`),  `y` must contains those `group_vars` column names
#'  (and this would allow different averaging periods for each group)
#' @param interval_vars a length-2 character vector of column names in both `x` and `y`.
#' These column names specify columns in x and y that define
#' closed (inclusive) starting and ending intervals. The column name
#' specifying the lower-bound column must be specified first.
#' these columns in x and y must all be of the same class and either be integer or IDate.
#' The interval_vars character vector cannot be named. This is reserved for future use allowing
#' different interval_vars column names in x and y.
#' @param value_vars a character vector of column names in `x`. This specifies
#' the columns to be averaged.
#' @param group_vars A character vector of column names in both x and y.
#' The interaction of
#' these variables define groups in which averages of `x` values will be taken.
#' specifying subjects/monitors/locations within which to take averages.
#' By default this is `NULL`, in which case averages are taken over the entire `x`
#' dataset for each `y` period.
#' The group_vars character vector cannot be named. This is reserved for future use allowing
#' different interval_vars column names in x and y.
#' @param required_percentage This percentage of the duration of each (possibly group-specific)
#' `y` interval must be observed
#' and nonmissing for a specific `value_var` in `x` in order for the return table to
#' contain a nonmissing average of the `value_var` for that `y` interval.  If the percentage
#' of the nonmissing `value_var` observations is less than `required_percentage` an NA will be be returned
#' for that average.
#'  The default is 100, meaning that if \emph{any} portion of a `y` interval is either not recorded or
#'  missing in `x`, then the corresponding return row will contain a an NA for the average of that
#'  `value_var`.
#' @param skip_overlap_check by default, FALSE. setting this to TRUE will skip
#'  internal checks to make sure x intervals are non-overlapping within
#'   groups defined by group_vars.
#'    intervals in x must be non-overlapping,
#'    but you may want to skip this check if you've  already checked this because
#'    it is computationally intensive for large datasets.
#' @param verbose include printed timing information? by default, FALSE
#' @return returns a data.table object.
#' Rows of the return data.table correspond to intervals from y. i.e, the number
#' of rows of the return will be the number of rows of y.
#' Columns of the returned data.table are as follows: \cr
#' - grouping variables as specified in `group_vars` \cr
#' - interval columns corresponding to intervals in y. These columns are named the
#'   same they were in x and y and as specified in `interval_vars`
#' - value variable columns from x, averaged to periods in y.
#'    named the same as they were in x \cr
#' - \code{yduration}: the length of the interval (ie as a count) specified in y \cr
#' - \code{xduration}: the total length of the intervals (ie as a count)
#'   from x that fall into this interval from y. this will be equal to
#'   yduration if x is comprehensive for (ie, fully covers)  this interval from y. \cr
#' - \code{nobs_<value_vars>}: for each \code{value_var} specified, this is the count of
#'  non-missing values from x that fall into this interval from y. this will be
#'   equal to xduration if the value_var contains no NA values over the y
#'   interval. If there are NAs in value variables, then \code{nobs_<value_vars>}
#'    will be different from \code{xduration} and won't necessarily be all the same
#'     for each value_var.
#' - \code{xminstart}: For each returned interval (ie the intervals from Y) the minimum of the
#' start intervals represented in x.  If the start of the earliest x interval is less than the start
#' of the y interval, the minimum of the y interval is returned. Note, this is the minimum start
#'  time in x matching with the y interval whether or not any value_vars were missing or not for that start time.
#'  If you need non-missing minimum start times, you could remove NA intervals from
#'  x prior to calling intervalaverage (this would need to be done separately for each value_var).
#' - \code{xmaxend}:  similar to xminstart but the maximum of the end intervals represented in x.
#'  Again, this does not pay attention to whether the interval in x had non-missing value_vars.
#' @examples
#'x <- data.table(start=seq(1L,by=7L,length=6),
#'                end=seq(7L,by=7L,length=6),
#'                pm25=c(10,12,8,14,22,18))
#'
#'y <- data.table(start=seq(3L,by=7L,length=6),
#'                end=seq(9L,by=7L,length=6))
#'
#'z <- intervalaverage(x,y,interval_vars=c("start","end"),
#'                     value_vars=c("pm25"))
#'
#' #also see vignette for more extensive examples
#' @export
intervalaverage <- function(x,
                            y,
                            interval_vars,
                            value_vars,
                            group_vars=NULL,
                            required_percentage=100,
                            skip_overlap_check=FALSE,
                            verbose=FALSE
){

  #due to NSE: to avoid notes in R CMD BUILD
  intervalaverage__xstart_copy <- intervalaverage__xend_copy <- intervalaverage__xend_copy <-
    intervalaverage__ystart_copy <-  intervalaverage__yend_copy <-
    yduration  <- xminstart <- xmaxend <- NULL



  stopifnot(data.table::is.data.table(x))
  stopifnot(data.table::is.data.table(y))


  EVAL <- function(...)eval(parse(text=paste0(...)))


  #save states so that added columns are dropped.
  #x also may be reordered due to the overlap check. this will restore the order
  statex <- savestate(x)
  on.exit(setstate(x,statex))
  statey <- savestate(y)
  on.exit(setstate(y,statey),add=TRUE)


  stopifnot(is.null(names(interval_vars)))
  stopifnot(is.null(names(group_vars)))

  if( any(c("yduration","xduration","xminstart","xmaxend")%in% c(interval_vars,value_vars,group_vars))){
    stop(paste0("column(s) named 'yduration', 'xduration', 'xminstart', or 'xmaxend' has been detected in interval_vars,",
                " value_vars, or group_vars.         These column names ('yduration', 'xduration', 'xminstart', or 'xmaxend')",
                " are reserved for the output. please rename this (or these) column(s) in x/y/interval_vars/value_vars/group_vars."))
  }

  if(!is.null(group_vars)){
    if(!all(group_vars %in% names(x))){
      stop("every value in group_vars must be a columname in x")
    }
    if(!all(group_vars %in% names(y))){
      stop("every value in group_vars must be a columname in y")
    }
  }


   if(!all(interval_vars %in% names(x))){
     stop("every value in interval_vars must be a columname in x")
   }
   if(!all(interval_vars %in% names(y))){
     stop("every value in interval_vars must be a columname in y")
   }

  if(!all(value_vars %in% names(x))){
    stop("every value in value_vars must be a columname in x")
  }


  if(x[,any(sapply(.SD,function(m){any(is.na(m))})) ,.SDcols=interval_vars]){
    stop("columns corresponding to interval_vars cannot be missing in x")
  }

  if(y[,any(sapply(.SD,function(m){any(is.na(m))})) ,.SDcols=interval_vars]){
    stop("columns corresponding to interval_vars cannot be missing in y")
  }

  if(x[,!all(sapply(.SD,is.integer)|sapply(.SD,function(x){class(x)%in% c("IDate")})),.SDcols=interval_vars]){
    stop("interval_vars must correspond to columns in x of class integer or IDate")
  }
  if(x[,class(.SD[[1]])[1]!=class(.SD[[2]])[1],.SDcols=interval_vars]){
    stop("interval_vars must correspond to columns in x of the same class")
  }

  if(y[,!all(sapply(.SD,is.integer)|sapply(.SD,function(x){any(class(x)%in%c("IDate"))})),.SDcols=interval_vars]){
    stop("interval_vars must correspond to columns in y of class integer or IDate")
  }
  if(y[,class(.SD[[1]])[1]!=class(.SD[[2]])[1],.SDcols=interval_vars]){
    stop("interval_vars must correspond to columns in y of the same class")
  }

  #stop if there are variables specified in both groups and interval_vars
  if(!length(intersect(interval_vars,group_vars))==0){
    stop("interval_vars and group_vars cannot refer to the same column(s)")
  }

  if(!length(intersect(value_vars,group_vars))==0){
    stop("value_vars and group_vars cannot refer to the same column(s)")
  }

  if(!length(intersect(value_vars,interval_vars))==0){
    stop("value_vars and interval_vars cannot refer to the same column(s)")
  }


  #stop if interval starts are before interval ends
  if(x[, sum(.SD[[2]]-.SD[[1]] <0)!=0,.SDcols=interval_vars]){
    stop("there exist values in x[[interval_vars[1] ]] that are
         less than corresponding values in  x[[interval_vars[2] ]].
         interval_vars must specify columns corresponding to increasing intervals")
  }

  #check for exact overlaps in x
  if(sum(duplicated(x[,c(group_vars,interval_vars),with=FALSE]))!=0){
    stop("sum(duplicated(x[,c(group_vars,interval_vars),with=FALSE]))!=0 is not TRUE.
         there are replicate/duplicate intervals within groups.
         If you wish to average these together, then do this first")
  }


  ydups <- duplicated(y[,c(group_vars,interval_vars),with=FALSE])
  if(sum(ydups)!=0){
    warning("sum(duplicated(y[,c(group_vars,interval_vars),with=FALSE]))!=0 is not TRUE.
         there are replicate/duplicate intervals within groups of y.
         removing these duplicated rows automatically")
  }else{
    ydups <- FALSE
  }

  if(!skip_overlap_check){

    #stop if there are overlapping periods within groups:
    data.table::setkeyv(x,c(group_vars,interval_vars))
    stopifnot(nrow(data.table::foverlaps(x,x))==nrow(x))
    if(verbose){print(paste(Sys.time(),"passed errorcheck: x is non-overlapping."))}
  }else{
    message("skipping errorcheck. if intervals in x are  overlapping, incorrect results may be returned without error.")
  }


  ### merge x and y ####

  xx <- proc.time()


  ##coerce date to numeric. doing it once here is faser than a million coercison in the by-loop
  interval_dates <- "IDate" %in% class(x[[interval_vars[1]]])
  if(interval_dates){
    stopifnot(!any(c("intervalaverage__start__original_class_copy",
                     "intervalaverage__end__original_class_copy")%in% names(x)))
    stopifnot(!any(c("intervalaverage__start__original_class_copy",
                     "intervalaverage__end__original_class_copy")%in% names(y)))

    setnames(x,interval_vars,
             c("intervalaverage__start__original_class_copy",
               "intervalaverage__end__original_class_copy")
    )
    setnames(y,interval_vars,
             c("intervalaverage__start__original_class_copy",
               "intervalaverage__end__original_class_copy")
    )

    x[, (interval_vars):=lapply(.SD,as.integer),
      .SDcols=c("intervalaverage__start__original_class_copy",
                "intervalaverage__end__original_class_copy")]
    y[, (interval_vars):=lapply(.SD,as.integer),
      .SDcols=c("intervalaverage__start__original_class_copy",
                "intervalaverage__end__original_class_copy")]

    ##set the columns back as they were
    on.exit({
      x[, (interval_vars):=NULL]
      y[, (interval_vars):=NULL]

      setnames(x, c("intervalaverage__start__original_class_copy",
                    "intervalaverage__end__original_class_copy"),
               interval_vars
      )
      setnames(y,c("intervalaverage__start__original_class_copy",
                   "intervalaverage__end__original_class_copy"),
               interval_vars,
      )


    }, add=TRUE, after=FALSE)

  }


  #copy the interval columns so they get carried over to the result
   #if you didn't copy them, the resulting joined columns would be *just the values from*
    #but we need the values from intervals in both x and y to calculate durations

  stopifnot(!any(c("intervalaverage__xstart_copy",
              "intervalaverage__xend_copy")%in% c(names(y),names(x))))

  stopifnot(!any(c("intervalaverage__ystart_copy",
                   "intervalaverage__yend_copy")%in% c(names(x),names(y))))

  x[,`:=`(intervalaverage__xstart_copy=.SD[[1]],
          intervalaverage__xend_copy=.SD[[2]]
          ),.SDcols=interval_vars]

  y[,`:=`(intervalaverage__ystart_copy=.SD[[1]],
          intervalaverage__yend_copy=.SD[[2]]
          ),.SDcols=interval_vars]

  #these copies will be deleted on function exit thanks to the on.exit state restore calls

  value_names <- unlist(lapply(value_vars, function(x) paste0(c("","nobs_"),x)))
  nobs_vars_names <- paste0("nobs_",value_vars)
  q <- x[y[!ydups],
    {

      Cintervalaverage(.SD,
                            intervalaverage__xstart_copy,
                            intervalaverage__xend_copy,
                            intervalaverage__ystart_copy[1],
                            intervalaverage__yend_copy[1],
                            value_names
      )



    },
    by=.EACHI,
    on=c(group_vars,
         paste0(interval_vars[2],">=",interval_vars[1]),
         paste0(interval_vars[1],"<=",interval_vars[2])),
    .SDcols=value_vars]


  #the on nonequi join seems confusing but remember the left side of the corresponds to vars in x
  #and the right side corresponds to vars in y.
    #so this just means take rows meeting both of these conditionds:
      #-the x end is greater than the y start AND
      #-the x start is less than or equal to the y end
    #which just translates to a partial overlap join
    #in addition to matching on groups_vars

    #to x and the right sid

  ###confusingly the returned non-equi join columns are the names from x but the values from y
  #because of this inequality statement
    #comparing interval_vars[2] with interval_val[1] and vice versa
  #the names get switched. switch them back:

  setnames(q, interval_vars, rev(interval_vars))

  #count length of each interval,
  #this will be used to count percentage of observed time x has in intervals of y
  q[,yduration:=as.numeric( .SD[[2]]-.SD[[1]] + 1),.SDcols=c(interval_vars)]
  #this should probably be an integer
  #but was historically a numeric so leaving this as is for now


  #fix column type of dates since .Internal(pmin) converts to numeric
  if(interval_dates){
    q[, xminstart:=as.IDate(xminstart,origin="1970-01-01")]
    q[, xmaxend:=as.IDate(xmaxend,origin="1970-01-01")]
    q[, (interval_vars):=lapply(.SD,as.IDate),.SDcols=interval_vars]
  }



  for(i in 1:length(value_vars)){
    #remove averages with too few observations in period
    #e.g. q[100*nobs_value/yduration < required_percentage, nobs_value:=NA]
    EVAL(paste0("q[100*",nobs_vars_names[i],"/yduration < required_percentage, ", value_vars[i],":=as.numeric(NA)]"))

    #turn NaNs into NAs
    #NaNs occur when rows of y were not matched at all to x
    set(q, i=which(is.nan(q[[value_vars[i]]])),j=value_vars[i],value=as.numeric(NA))
  }

  data.table::setcolorder(q, c(group_vars,interval_vars,value_vars,"yduration","xduration",nobs_vars_names,
                   "xminstart","xmaxend"))

  setkey(q,NULL) #temporary way to get around this bug: https://github.com/Rdatatable/data.table/issues/4603
  setkeyv(q, c(group_vars,interval_vars))


  if(verbose){
    print("everything else:")
    print(proc.time()-xx)
  }
  q[]

}



#slower algorithm. used for testing since this is the simpler approach and less likely to have errors.
#Not recommended for large datasets since this expands the data into a row for every increment.
#not exported
interval_weighted_avg_slow_f <- function(x,
                                         y,
                                         interval_vars,
                                         value_vars,
                                         group_vars=NULL,
                                         required_percentage=100,
                                         skip_overlap_check=FALSE,
                                         verbose=FALSE){

  xminstart <- xmaxend <- NULL

  if(x[,!all(sapply(.SD,is.integer)|sapply(.SD,function(x){class(x)%in% c("IDate")})),.SDcols=interval_vars]){
    stop("interval_vars must correspond to columns in x of class integer or IDate")
  }
  if(x[,class(.SD[[1]])[1]!=class(.SD[[2]])[1],.SDcols=interval_vars]){
    stop("interval_vars must correspond to columns in x of the same class")
  }

  if(y[,!all(sapply(.SD,is.integer)|sapply(.SD,function(x){any(class(x)%in%c("IDate"))})),.SDcols=interval_vars]){
    stop("interval_vars must correspond to columns in y of class integer or IDate")
  }
  if(y[,class(.SD[[1]])[1]!=class(.SD[[2]])[1],.SDcols=interval_vars]){
    stop("interval_vars must correspond to columns in y of the same class")
  }

  EVAL <- function(...)eval(parse(text=paste0(...)))


  if(!is.null(group_vars)){
    if(!all(group_vars %in% names(x))){
      stop("every value in group_vars must be a columname in x")
    }
    if(!all(group_vars %in% names(y))){
      stop("every value in group_vars must be a columname in y")
    }
  }

  #check for exact overlaps
  if(sum(duplicated(x[,c(group_vars,interval_vars),with=FALSE]))!=0){
    stop("sum(duplicated(x[,c(group_vars,interval_vars),with=FALSE]))!=0 is not TRUE.
         there are replicate/duplicate intervals within groups.
         If you wish to average these together, then do this first")
  }


  ydups <- duplicated(y[,c(group_vars,interval_vars),with=FALSE])
  if(sum(ydups)!=0){
    warning("sum(duplicated(y[,c(group_vars,interval_vars),with=FALSE]))!=0 is not TRUE.
         there are replicate/duplicate intervals within groups of y.
         removing these duplicated rows automatically")
    y <- y[!ydups]
  }




  #set keys for testing if there are overlaps
  #groups do not need to be in x but they do need to be in y
  if(!skip_overlap_check){

    data.table::setkeyv(x,c(group_vars,interval_vars))

    #stop if there are overlapping periods within groups:
    stopifnot(nrow(data.table::foverlaps(x,x))==nrow(x))
    if(verbose){print(paste(Sys.time(),"passed errorcheck: x is non-overlapping."))}
  }else{
    message("skipping errorcheck. if intervals in x are  overlapping, incorrect results may be returned without error.")
  }



  #create a length-1 character vector that will be a column name in x,y, and z that is not in use already
  t <- create_unused_name("time",c(names(x),names(y)))

  #create a length-1 character vector that will be a column name in x,y, and z that is not in use already
  # measurement occured on this day but it might be measured as missing
  measurement <- create_unused_name("meas",c(names(x),names(y)))



  if("yduration"%in%names(x)|"yduration"%in%names(y)){
    stop("names(x) or names(y) contains the column name 'yduration'. A column named 'yduration' cannot be present in x or y
         because it will be a special column in the return value. please rename this column.")
  }
  ydur <- "yduration"

  if("xduration"%in%names(x)|"xduration"%in%names(y)){
    stop("names(x) or names(y) contains the column name 'xduration'. A column named 'xduration' cannot be present in x or y
         because it will be a special column in the return value. please rename this column.")
  }
  xdur <- "xduration"


  if("xminstart"%in%names(x)|"xminstart"%in%names(y)){
    stop("names(x) or names(y) contains the column name 'xduration'. A column named 'xminstart' cannot be present in x or y
         because it will be a special column in the return value. please rename this column.")
  }

  if("xmaxend"%in%names(x)|"xmaxend"%in%names(y)){
    stop("names(x) or names(y) contains the column name 'xduration'. A column named 'xmaxend' cannot be present in x or y
         because it will be a special column in the return value. please rename this column.")
  }


  #nobs vars will be the count of non-missing obs for each time-period in y
  #(ie, summed from temp nobs)
  nobs_vars <- paste("nobs",value_vars, sep="_")
  for(i in 1:length(nobs_vars)){
    while(nobs_vars[i]%in% names(x)|nobs_vars[i]%in% names(y)){
      ydur <- paste0("i.",nobs_vars[i])
    }
  }


  #expand x
  #in x values_vars represent an average over an interval
  #in expand_x, intervals are expanded such that there is now a separate row for each smallest observable
  #increment in each interval (represented as variable t).
  #values are repeated over rows that correspond to those increments,
  #and interval_vars columns specifying the original increments are retained

  x_expanded <- EVAL("x[,list(",t,"=",interval_vars[1],":",interval_vars[2],"),
                     by=c(group_vars,interval_vars,value_vars)]")

  x_expanded[,(interval_vars):=NULL]
  x_expanded[, (measurement):=1L]

  y_expanded <-  EVAL("y[,list(",t,"=",interval_vars[1],":",interval_vars[2],"),
                      by=c(interval_vars,group_vars)]")



  data.table::setkeyv(x_expanded,c(t, group_vars))
  data.table::setkeyv(y_expanded,c(t, group_vars))

  z <- x_expanded[y_expanded]

  avg_call <- paste0(value_vars,"=mean(",value_vars,",na.rm=TRUE)",collapse=",")
  paste0()
  nobs_call <- paste0(nobs_vars,"=sum(!is.na(",value_vars,"))",collapse=",")
  ydur_call <- paste0(ydur,"=length(unique(",t,"))")
  xdur_call <- paste0(xdur,"=sum(",measurement,",na.rm=TRUE)")

  #get the min start and max end dates for intervals that were actually provided as inputs in x
  #as of 1/12/2020 I think gforce does not work with :=,
    #so this is written to create a separate table that is then merged
  minmaxtable <-
    EVAL(
      paste0(
       "z[!is.na(",measurement,"),list(xminstart=min(time),xmaxend=max(time)),by=c(interval_vars,group_vars)]"
      )
    )

  if(any(class(x[[interval_vars[1]]])%in%c("IDate"))){
    minmaxtable[, xminstart:=as.IDate(xminstart,origin="1970-01-01")]
    minmaxtable[, xmaxend:=as.IDate(xmaxend,origin="1970-01-01")]
  }



  out <- EVAL(
    paste0(
    "z[,","list(",avg_call,",",ydur_call,",",xdur_call,",",nobs_call,")" ,
              ",by=c(interval_vars,group_vars)]"
    )
    )

  stopifnot(EVAL("out[,all(",xdur,"<=",ydur,")]"))

  ##merge in minmaxtable
  data.table::setkeyv(minmaxtable,c(group_vars,interval_vars))
  data.table::setkeyv(out,c(group_vars,interval_vars))

  out <- minmaxtable[out]


  for(i in 1:length(value_vars)){
    EVAL("out[100*",nobs_vars[i],"/",ydur,"< required_percentage,",value_vars[i],":=NA]")
  }

  data.table::setcolorder(out, c(group_vars,interval_vars,value_vars,ydur,xdur,nobs_vars,
                     "xminstart","xmaxend"))

  data.table::setkeyv(out,c(group_vars,interval_vars))
  out[]
}



#' isolate sections of overlapping intervals
#'
#' Given a set of intervals in a table, isolate sections of intervals that are overlapping
#' with other in intervals (optionally, within groups). Returns a data.table that contains
#' intervals which are mutually non-overlapping or exactly overlapping with other intervals
#' (ie there are no partially overlapping intervals) (optionally within groups).
#' Note that this doesn't just return the intersects; the original interval data is conserved
#' such that for each interval/row in x, the return table has one or more
#' non-overlapping intervals that together form the union of that original interval.
#'
#' All intervals are treated as closed (ie inclusive of the start and end values in the columns
#' specified by interval_vars)
#'
#' x is not copied but rather passed by reference to function internals
#' but the order of this data.tables is restored on function completion or error.
#'
#' @param x A data.table containing a set of intervals.
#' @param interval_vars A length-2 character vector denoting column names in x.
#' these columns must be of the same class and be integer or IDate. The column name
#' specifying the lower-bound column must be specified first.
#' @param group_vars NULL, or a character vector denoting column names in x.
#'  These columns serve as grouping variables such that testing for overlaps and subsequent isolation only occur
#'  within categories defined by the combination of the group variables.
#' @param interval_vars_out The desired column names of the interval columns in the return data.table.
#' By default these columns will be generated to be named \code{c("start","end")}.
#' If x contains columns with the same name as the desired output column
#' names specified in `interval_vars_out`, the function will return in error to avoid naming
#' confusion in the return table.
#' @return A data.table with columns \code{interval_vars_out} which denote the start and
#' stop period for each new interval. This return table also contains columns in x
#' (including the original interval columns).
#' @seealso \code{\link{is.overlapping}} To test if a table contains overlapping intervals within
#'  values of \code{group_vars}
#'
#' @examples
#'set.seed(23)
#'x2 <- data.table(addr_id=rep(1:4,each=3),
#'                 exposure_start=rep(c(1L,7L,14L),times=4),
#'exposure_end=rep(c(7L,14L,21L),times=4),
#'exposure_value=c(rnorm(12))
#')
#'x2z <- isolateoverlaps(x2,interval_vars=c("exposure_start","exposure_end"),group_vars=c("addr_id"))
#'x2z
#'#x2b represents x2 when where exposure values in overlapping intervals have been averaged
#'x2b <- x2z[, list(exposure_value=mean(exposure_value)),by=c("addr_id","start","end")]
#'
#' @export
isolateoverlaps <- function(x,interval_vars,group_vars=NULL,interval_vars_out=c("start","end")){
  EVAL <- function(...)eval(parse(text=paste0(...)))

  stopifnot(is.data.table(x))

  if(x[,!all(sapply(.SD,is.integer)|sapply(.SD,function(x){class(x)%in% c("IDate")})),.SDcols=interval_vars]){
    stop("interval_vars must correspond to columns in x of class integer or IDate")
  }
  if(x[,class(.SD[[1]])[1]!=class(.SD[[2]])[1],.SDcols=interval_vars]){
    stop("interval_vars must correspond to columns in x of the same class")
  }


  if(any(interval_vars_out %in% names(x))){ stop("Column names in x detected to have the same values as interval_vars_out.
                                            This would causes naming conflicts in the return. Choose different
                                            output names by specifying the argument interval_vars_out")}
  variable_var <- create_unused_name("variable",names(x))
  i.interval_vars <- paste0("i.", interval_vars)
    if(any(i.interval_vars%in% names(x))){
      stop("names(x) contains columns ",paste0(interval_vars,collapse=" "),
      " and at least one column named ",paste0(i.interval_vars,collapse=" "),". Columns named ",paste0(i.interval_vars,collapse=" "),
      " cannot be present in x because they are reserved for use by data.table::foverlaps")
      }

  value_var <- create_unused_name("value",names(x))
  open_var <- create_unused_name("open",names(x))
  close_var <- create_unused_name("close",names(x))
  is_end_var <- create_unused_name("is_end",names(x))

  end_next_var <- create_unused_name("end_next",names(x))
  value_next_var <- create_unused_name("value_next",names(x))

  xd <- data.table::melt(x[,c(interval_vars,group_vars),with=FALSE],id.vars=group_vars,
             measure.vars=interval_vars,
             value.name=value_var,
             variable.name=variable_var)



  ##create end variable:
  EVAL("xd[",variable_var,"==interval_vars[2],(is_end_var):=TRUE]")
  EVAL("xd[",variable_var,"==interval_vars[1],(is_end_var):=FALSE]")
  data.table::setorderv(xd,c(group_vars, value_var,is_end_var))

  EVAL("xd[,(end_next_var):=shift(",is_end_var,",type='lead'),by=group_vars]")
  EVAL("xd[,(value_next_var):=shift(",value_var,",type='lead'),by=group_vars]")

  #the last row is missing because you can't get the next row when there aren't any more rows!
   #we don't need that row where end_next is missing so exclude it.
  #when data.table veresion 1.12.3 you can use fifelse to avoid coercing dates to numeric
  temp <- EVAL("xd[,.SD[!is.na(",end_next_var,"),list(
   ",interval_vars[1],"=data.table::fifelse(!",is_end_var,",",value_var,",",value_var,"+1L),
    ",interval_vars[2],"=data.table::fifelse(!",end_next_var,",",value_next_var,"-1L,",value_next_var,")
  )],by=group_vars]")

  temp <- EVAL("temp[",interval_vars[2],">=",interval_vars[1],"]")


  data.table::setkeyv(temp,c(group_vars,interval_vars))


  out <- data.table::foverlaps(x,temp,by.x=c(group_vars,interval_vars))
  data.table::setorderv(out, c(group_vars,
                   paste0("i.",interval_vars),
                   interval_vars))
  data.table::setcolorder(out, c(group_vars,interval_vars,paste0("i.",interval_vars)))

  data.table::setnames(out, interval_vars, interval_vars_out)
  data.table::setnames(out, i.interval_vars, interval_vars)


  out
}

