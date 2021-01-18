

#' Intersect intervals within groups
#'
#' Given two tables each containing a set of intervals, find all interval intersections within groups.
#' Returns a data.table containing all columns from both tables.
#' One use of this function is to take a table  containing an address history
#' (a table containing the intervals when study participants lived at past addresses)
#' and join it to an exposure history table (a complete set of exposure predictions
#' for each address, where the exposures are stored as the average value over a
#' set of intervals) returning the set of
#' exposure intervals at addresses clipped to exactly when the participant lived at that address.
#'
#' All intervals are treated as closed (ie inclusive of the start and end values in the columns
#' specified by interval_vars)
#'
#' x and y are not copied but rather passed by reference to function internals
#' but the order of these data.tables is restored on function completion or error,
#'
#' Technically speaking this is just an inner cartesian join where the last two join variables are
#' doing a non-equi join for partial overlaps. Then each interval intersect is calculated using max and min.
#'
#' If there are columns with the same names in both x and y (including interval_vars
#' but excepting group_vars), the return value will still return both columns. The column in y
#' will be names as it was originally and the column in x will be prepended with the letter i followed with a dot:
#' \code{i.} \cr
#'
#' Note that the function returns the same result if you switch x and y
#'  (with the exception of switched column names in the case of column name conflicts as just discussed)
#'
#'
#' @param x A data.table with two columns defining closed intervals (see also interval_vars parameter)
#' @param y A data.table with two columns defining closed intervals (see also interval_vars parameter)
#' @param interval_vars Either a length-2 character vector denoting column names in both x and y or a named
#'  length-2 character vector where the names are column names in x and the values are column names in y.
#' These column names specify columns in x and y that define
#' closed (inclusive) starting and ending intervals. The column name
#' specifying the lower-bound column must be specified first.
#' these columns in x and y must all be of the same class and either be integer or IDate.
#' @param group_vars NULL, or either a character vector denoting the column name(s) in x and y,
#'  or a named character vector where the name is the column name in x and the value is the column name in y.
#'  This/these column(s) serve as an additional keying variable in the join (ie in addition to the interval join)
#'  such that intervals in x will only be joined to overlapping in intervals in y where the group_vars values are the same.
#' @param interval_vars_out The column names of the interval columns in the return data.table.
#' By default the return table will contain columns \code{c("start","end")}.
#' If your input tables already contain these columns,
#' you need to either specify \code{interval_vars_out} to be non-conflicting names with columns in x and y.
#' Or you rename columns in x and y to not contain columns named \code{c("start","end")}.
#' @param verbose Prints additional information about the function processing.
#' @return A data.table with columns \code{interval_vars_out} which denote the start and
#' stop period for each interval. This return table also contains columns in x and y. See details
#' for how naming conflicts are dealt with.
#' @seealso \code{\link{is.overlapping}} To test if a table contains overlapping intervals within values of \code{group_vars}
#' @examples
#'set.seed(42)
#'y <- data.table(addr_id=c(1,2,2,3,5),
#'ppt_id=c(1,1,1,2,2),
#'addr_start=c(1L,10L,12L,1L,1L),
#'addr_end=c(9L,11L,14L,17L,10L))
#'x <- data.table(addr_id=rep(1:4,each=3),
#'exposure_start=rep(c(1L,8L,15L),times=4),
#'exposure_end=rep(c(7L,14L,21L),times=4),
#'exposure_value=c(rnorm(12))
#')
#'intervalintersect(x,y,
#'interval_vars=c(exposure_start="addr_start",exposure_end="addr_end"),
#'"addr_id")
#'y2 <- data.table(addr_id=c(1,2,2,2,3),
#'ppt_id=c(1,1,1,1,2),
#'addr_start=c(1L,2L,3L,4L,1L),
#'addr_end=c(9L,12L,13L,8L,10L))
#'
#'#intervalintersect will still work when there are overlapping intervals within a table:
#'is.overlapping(y2,interval_vars =c("addr_start","addr_end") ,group_vars="addr_id")
#'
#'intervalintersect(x,y2,
#'                  interval_vars=c(exposure_start="addr_start",exposure_end="addr_end"),
#'"addr_id")
#'
#'
#'x2 <- data.table(addr_id=rep(1:4,each=3),
#'exposure_start=rep(c(1L,7L,14L),times=4),
#'                 exposure_end=rep(c(7L,14L,21L),times=4),
#'exposure_value=c(rnorm(12))
#')
#'is.overlapping(x2,interval_vars =c("exposure_start","exposure_end") ,group_vars="addr_id")
#'
#'intervalintersect(x2,y2,
#'interval_vars=c(exposure_start="addr_start",exposure_end="addr_end"),
#'"addr_id")
#' #however, it may be meaningful isolate intervals of partial-overlap within
#' #each table and deal with them
#' #prior to intersecting the tables together:
#'
#'x2z <- isolateoverlaps(x2,interval_vars=c("exposure_start","exposure_end"),group_vars=c("addr_id"),
#'interval_vars_out=c("exposure_start2","exposure_end2"))
# x2b represents x2 when where exposure values in overlapping intervals have been averaged
#'x2b <- x2z[, list(exposure_value=mean(exposure_value)),
#'  by=c("addr_id","exposure_start2","exposure_end2")]
#'data.table::setnames(x2b, c("exposure_start2","exposure_end2"),c("exposure_start","exposure_end"))
#'
#'y2z <- isolateoverlaps(y2,interval_vars=c("addr_start","addr_end"),group_vars=c("addr_id"),
#'interval_vars_out = c("addr_start2","addr_end2"))
#'y2b <- unique(y2z[, list(addr_id, ppt_id,addr_start2,addr_end2)])
#'data.table::setnames(y2b, c("addr_start2","addr_end2"), c("addr_start","addr_end"))
#y2b represents where address history which are overlapping in the same address are deduped:
#'
#'intervalintersect(x2b,y2b,
#'interval_vars=c(exposure_start="addr_start",exposure_end="addr_end"),
#'"addr_id")
#'
#' @export
intervalintersect <- function(x,y, interval_vars, group_vars=NULL, interval_vars_out=c("start","end"),verbose=FALSE){



  x_interval_vars <- if(!is.null(names(interval_vars))){names(interval_vars)}else{interval_vars}
  y_interval_vars <- interval_vars


  #restore remove any added columns
  #using statefunctions isn't needed here because x and y aren't reordered.
   #only columns are added which need to be removed
  orig_colnames_x <-  copy(names(x))
  orig_colnames_y <-  copy(names(y))
  on.exit({
    newcolsx <- setdiff(names(x),orig_colnames_x)
    set(x,j=newcolsx,value=NULL)
    newcolsy <- setdiff(names(y),orig_colnames_y)
    set(y,j=newcolsy,value=NULL)
  })

  if(x[,!all(sapply(.SD,is.integer)|sapply(.SD,function(x){class(x)%in% c("IDate")})),.SDcols=x_interval_vars]){
    stop("interval_vars must correspond to columns in x of class integer or IDate")
  }
  if(x[,class(.SD[[1]])[1]!=class(.SD[[2]])[1],.SDcols=x_interval_vars]){
    stop("interval_vars must correspond to columns in x of the same class")
  }

  if(y[,!all(sapply(.SD,is.integer)|sapply(.SD,function(x){any(class(x)%in%c("IDate"))})),.SDcols=y_interval_vars]){
    stop("interval_vars must correspond to columns in y of class integer or IDate")
  }
  if(y[,class(.SD[[1]])[1]!=class(.SD[[2]])[1],.SDcols=y_interval_vars]){
    stop("interval_vars must correspond to columns in y of the same class")
  }

  if( any(interval_vars_out %in% names(x))|any(interval_vars_out %in% names(y)) ){
    stop("interval_vars_out cannot be names in x or y. choose a different output name")
  }



  x_group_vars <- if(!is.null(names(group_vars))){names(group_vars)}else{group_vars}
  y_group_vars <- group_vars


  ##generate copies of interval_columns so they're not lost in the non-equijoin
  stopifnot(!any(c("intervalaverage__xstart_copy",
                   "intervalaverage__xend_copy")%in% c(names(y),names(x))))

  stopifnot(!any(c("intervalaverage__ystart_copy",
                   "intervalaverage__yend_copy")%in% c(names(x),names(y))))

  x[,`:=`(intervalaverage__xstart_copy=.SD[[1]],
          intervalaverage__xend_copy=.SD[[2]]
  ),.SDcols=x_interval_vars]

  y[,`:=`(intervalaverage__ystart_copy=.SD[[1]],
          intervalaverage__yend_copy=.SD[[2]]
  ),.SDcols=y_interval_vars]
  #these columns will be deleted on exit of the function. see on.exit above


  #construct on statement:
  group_vars_on <- paste0(x_group_vars,"==",y_group_vars)
  #   paste0(interval_vars[2],">=",interval_vars[1]),
  #paste0(interval_vars[1],"<=",interval_vars[2])

  interval_vars_on <- c(
    paste0(x_interval_vars[2],">=",y_interval_vars[1]),
    paste0(x_interval_vars[1],"<=",y_interval_vars[2])

  )

  z <- x[y,
    on=c(group_vars_on,interval_vars_on),nomatch=NULL]


  #thanks to data.table's non-equijoin weirdness, the non-equijoin columns in z
    #are named as they were in x but have the values that were in y
  #since we made copies of the original values, the join columns are not needed:
  z[, (x_interval_vars):=NULL]


  #take the later of the two start dates and the earlier of the two end dates
  data.table::set(z, j=interval_vars_out[1],
                  value=pmax(z[["intervalaverage__xstart_copy"]],z[["intervalaverage__ystart_copy"]] ) )
  data.table::set(z, j=interval_vars_out[2],
                  value=pmin(z[["intervalaverage__xend_copy"]],z[["intervalaverage__yend_copy"]] ) )

  z[,c("intervalaverage__ystart_copy",
       "intervalaverage__yend_copy"):=NULL]
  z[,c("intervalaverage__xstart_copy",
       "intervalaverage__xend_copy"):=NULL]

  other_vars <- setdiff(names(z),c(interval_vars_out,x_group_vars))
  setcolorder(z, c(x_group_vars, interval_vars_out, other_vars))
  setkey(z,NULL)
  setkeyv(z,c(x_group_vars,interval_vars_out))
  z[]
}
