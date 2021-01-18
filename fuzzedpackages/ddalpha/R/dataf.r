dataf <- function(fdatas, labels) 
{
  .split_data_matrix <-function(fdata) {
    funcs = split(fdata$data, row(fdata$data))
    dataf = lapply(funcs, function(func) list(args = fdata$argvals, vals = func))
  }
  
  if (class(fdatas) == "fdata") {
    if (nrow(fdatas$data) != length(labels))
      stop("the length of 'labels' must correspond to the number of functions in 'fdatas'")
    if (!is.null(fdatas$fdata2d) && fdatas$fdata2d)
      stop("fdata2d = TRUE is not supported")
    
    dataf.labels = as.list(labels)
    dataf.dataf = .split_data_matrix(fdatas)
    res = list(
      dataf = dataf.dataf,
      labels = dataf.labels,
      name = fdatas$names$main,
      args = fdatas$names$xlab,
      vals = fdatas$names$ylab
    )
    class(res) <- "functional"
    return (res)
  }
  
  
  if (length(fdatas) != length(labels))
    stop("'fdatas' and 'labels' must be vectors of the same length")

  dataf.dataf = list()
  dataf.labels = c()
  for (i in range(1:length(fdatas))) {
    fdata = fdatas[[i]]
    lab = labels[[i]]
    
    if (class(fdata) != "fdata")
      stop("elements of 'fdatas' must be of the 'fdata' class")
    if (!is.null(fdata$fdata2d) && fdata$fdata2d)
      stop("fdata2d = TRUE is not supported")
  
    dataf.labels = c(dataf.labels, rep(lab, nrow(fdata$data)))
    
    dataf.dataf = c(dataf.dataf, .split_data_matrix(fdata))
  }  
  
  res = list(
    dataf = as.list(dataf.dataf),
    labels = dataf.labels,
    # just take the names from the first data set
    name = fdatas[[1]]$names$main,
    args = fdatas[[1]]$names$xlab,
    vals = fdatas[[1]]$names$ylab
  )
  class(res) <- "functional"
  return (res)
  
}


