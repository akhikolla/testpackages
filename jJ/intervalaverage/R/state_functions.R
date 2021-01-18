#functions to save and restore the row order of a data.table used as an argument
#these are useful when you don't want to copy the data.table withina  function
#but also need to sort it
#this allows you to make the function side-effet free (doesn't change the data.table)
 #by changing the data.table back on function exit
 #call savestate at the top of the function body and setstate at the bottom
 #use call setstate in an error trycatch to restore state on error

#these functions are not exported

##exclude names for temporary variable choice

savestate <- function(x){
  original_colnames <- copy(names(x))
  rowindex_colname <- "rowindex"
  while(rowindex_colname %in% names(x)) rowindex_colname <- paste0("i.",rowindex_colname)
  x[,eval(rowindex_colname):=1:.N]

  k <- key(x)
  list(key=k,rowindex_colname=rowindex_colname,original_colnames=original_colnames)
}



setstate <- function(x,state){

  required_columns <- c(state$original_colnames, state$rowindex_colname)
  stopifnot(all(required_columns %in% names(x)))

  setorderv(x,state$rowindex_colname)
  setkeyv(x,state$key)
  #explicitly remove index for clarity,
  #although this would happen anyway in the next line:
  x[, eval(state$rowindex_colname):=NULL]
  #remove all new columns:
  newcols <- setdiff(names(x),state$original_colnames)
  if(length(newcols)){
    set(x,j=newcols,value=NULL)
  }
  x
}
